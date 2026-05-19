//! pval — FFT-convolution p-value for genomic interval overlap.
//!
//! Null model: file B is a random rigid-body shift of itself on each chromosome
//! (same interval count, sizes, and relative gaps; only the starting offset is
//! uniformly random).  The test statistic is the total base-pair overlap (binned at
//! 100 bp) between A and the shifted B, summed over all chromosomes.
//!
//! Algorithm
//! ---------
//! For each chromosome c:
//!   1. Build the sparse derivative representation of each file's depth coverage:
//!        d[start_bin] += 1,  d[end_bin] -= 1
//!      so that prefix_sum(d) = g = the coverage-depth array.
//!   2. R2C FFT(d) → D (half spectrum, length n/2+1).  Recover the coverage
//!      spectrum via the integration identity:
//!        G[k] = D[k] / (1 − e^{−2πik/N})   for k ≠ 0
//!        G[0] = Σ_i g[i]                    (total coverage, computed directly)
//!   3. Cross-correlate in-place: G_A[k] *= conj(G_B[k])
//!      C2R IFFT → c[s] = base-pair overlap (in 100 bp bins) at circular shift s.
//!   4. Build a per-chromosome PMF: PMF_c[k] = fraction of shifts where c_c(s) = k.
//!
//! Combine chromosomes: shifts are independent, so
//!        PMF_total = PMF_1 ★ PMF_2 ★ … ★ PMF_C   (convolution = polynomial multiplication)
//!
//! P-value = Σ_{k ≥ observed} PMF_total[k]   (right-tailed).
//!
//! Performance notes
//! -----------------
//! * Real-to-Complex FFTs (realfft) halve both the memory and transform work for
//!   the spatial step; only the Hermitian half of each spectrum is materialised.
//! * All spatial arrays are f32; f64 is kept only for the PMF polynomial multiply
//!   where accumulation across ~24 chromosomes makes precision matter.
//! * Cross-correlation is done in-place (G_A reused); C2R output reuses the
//!   original real derivative buffer, eliminating hot-loop allocations.
//! * Per-chromosome work runs in parallel via Rayon; the sequential PMF fold
//!   follows once all chromosomes are done.

use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use realfft::RealFftPlanner;
use rayon::prelude::*;
use rustfft::{FftPlanner, num_complex::Complex};

// ── Constants ────────────────────────────────────────────────────────────────

const BIN_SIZE: u32 = 100; // bp per bin

// ── hg38 chromosome sizes ────────────────────────────────────────────────────

fn hg38_chrom_sizes() -> &'static [(&'static str, u32)] {
    &[
        ("chr1",  248_956_422),
        ("chr2",  242_193_529),
        ("chr3",  198_295_559),
        ("chr4",  190_214_555),
        ("chr5",  181_538_259),
        ("chr6",  170_805_979),
        ("chr7",  159_345_973),
        ("chr8",  145_138_636),
        ("chr9",  138_394_717),
        ("chr10", 133_797_422),
        ("chr11", 135_086_622),
        ("chr12", 133_275_309),
        ("chr13", 114_364_328),
        ("chr14", 107_043_718),
        ("chr15", 101_991_189),
        ("chr16",  90_338_345),
        ("chr17",  83_257_441),
        ("chr18",  80_373_285),
        ("chr19",  58_617_616),
        ("chr20",  64_444_167),
        ("chr21",  46_709_983),
        ("chr22",  50_818_468),
        ("chrX",  156_040_895),
        ("chrY",   57_227_415),
    ]
}

// ── BED parsing ──────────────────────────────────────────────────────────────

type Intervals = Vec<(u32, u32)>;
type BedMap = HashMap<String, Intervals>;

fn parse_bed(path: &str) -> BedMap {
    let f = File::open(path).unwrap_or_else(|e| panic!("Cannot open {path}: {e}"));
    let is_gz = Path::new(path).extension().and_then(|e| e.to_str()) == Some("gz");
    let reader: Box<dyn BufRead> = if is_gz {
        Box::new(BufReader::new(MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::new(f))
    };
    let mut map: BedMap = HashMap::new();
    for line in reader.lines() {
        let line = line.expect("I/O error");
        let t = line.trim();
        if t.is_empty()
            || t.starts_with('#')
            || t.starts_with("track")
            || t.starts_with("browser")
        {
            continue;
        }
        let mut cols = t.split_ascii_whitespace();
        let (Some(chrom), Some(s), Some(e)) = (cols.next(), cols.next(), cols.next()) else {
            continue;
        };
        let (Ok(start), Ok(end)) = (s.parse::<u32>(), e.parse::<u32>()) else {
            continue;
        };
        map.entry(chrom.to_string()).or_default().push((start, end));
    }
    map
}

// ── Core FFT helpers (f32 spatial) ───────────────────────────────────────────

/// Build a real f32 derivative spike array of length `n` and return it alongside
/// Σᵢ g[i] (total area under the depth curve over `chrom_bins` bins).
///
/// Each spike at raw position x is split between bins floor(x/BIN) and floor(x/BIN)+1
/// with weights (1−frac) and frac, where frac = (x % BIN_SIZE) / BIN_SIZE.
/// This eliminates aliasing for intervals shorter than BIN_SIZE.
/// prefix_sum(d) = g = per-bin coverage depth.
fn build_derivative(intervals: &[(u32, u32)], n: usize, chrom_bins: usize) -> (Vec<f32>, f32) {
    let max_bp = (chrom_bins as u32) * BIN_SIZE;
    let mut d = vec![0.0f32; n];
    for &(s, e) in intervals {
        let s = s.min(max_bp);
        let e = e.min(max_bp);
        let sb = (s / BIN_SIZE) as usize;
        let sf = (s % BIN_SIZE) as f32 / BIN_SIZE as f32;
        if sb < n { d[sb] += 1.0 - sf; }
        if sf > 0.0 && sb + 1 < n { d[sb + 1] += sf; }
        let eb = (e / BIN_SIZE) as usize;
        let ef = (e % BIN_SIZE) as f32 / BIN_SIZE as f32;
        if eb < n { d[eb] -= 1.0 - ef; }
        if ef > 0.0 && eb + 1 < n { d[eb + 1] -= ef; }
    }
    // Σ_i g[i] = Σ_j d[j] · (chrom_bins − j)
    // This formula works correctly with fractional d[j] values.
    let total_cov: f32 = d[..chrom_bins.min(n)]
        .iter()
        .enumerate()
        .filter(|(_, v)| **v != 0.0)
        .map(|(i, v)| v * (chrom_bins - i) as f32)
        .sum();
    (d, total_cov)
}

/// In-place: convert derivative R2C spectrum (length n/2+1) → coverage spectrum.
/// `n` is the full real-space length (needed for the twiddle denominator).
///
///   G[k] = D[k] / (1 − e^{−2πik/N})   for k ≠ 0
///   G[0] = total_cov
fn derivative_to_coverage(spec: &mut [Complex<f32>], total_cov: f32, n: usize) {
    spec[0] = Complex::new(total_cov, 0.0);
    for k in 1..spec.len() {
        let omega = -(std::f32::consts::TAU * k as f32) / n as f32;
        spec[k] /= Complex::new(1.0 - omega.cos(), -omega.sin());
    }
}

// ── Per-chromosome overlap ────────────────────────────────────────────────────

/// Compute c[s] = base-pair overlap (in 100 bp bins) between A and B circularly
/// shifted by s bins, for every s in [0, chrom_bins).
///
/// Memory layout:
///   da / db — real spike arrays (length n, f32); da is reused as C2R output.
///   sa / sb — R2C spectra (length n/2+1, Complex<f32>); sa holds the cross-
///             correlation spectrum in-place before the inverse transform.
fn chrom_overlap_vs_shift(
    a_ivs: &[(u32, u32)],
    b_ivs: &[(u32, u32)],
    chrom_size: u32,
) -> Vec<f32> {
    let n_bins = ((chrom_size + BIN_SIZE - 1) / BIN_SIZE) as usize;
    let n = n_bins.next_power_of_two();

    let mut planner = RealFftPlanner::<f32>::new();
    let r2c = planner.plan_fft_forward(n);
    let c2r = planner.plan_fft_inverse(n);

    let (mut da, cov_a) = build_derivative(a_ivs, n, n_bins);
    let (mut db, cov_b) = build_derivative(b_ivs, n, n_bins);

    // Forward R2C: real spikes → half spectra (length n/2+1)
    let mut sa = r2c.make_output_vec();
    let mut sb = r2c.make_output_vec();
    r2c.process(&mut da, &mut sa).unwrap();
    r2c.process(&mut db, &mut sb).unwrap();

    // Integrate: derivative spectrum → coverage spectrum
    derivative_to_coverage(&mut sa, cov_a, n);
    derivative_to_coverage(&mut sb, cov_b, n);

    // Cross-correlation in-place: sa[k] *= conj(sb[k])
    for (a, b) in sa.iter_mut().zip(sb.iter()) {
        *a *= b.conj();
    }

    // realfft requires the DC and Nyquist bins to be exactly real before C2R.
    // After the division chain and complex multiply, f32 leaves tiny residuals.
    sa[0].im = 0.0;
    let nyquist = sa.len() - 1;
    sa[nyquist].im = 0.0;

    // Inverse C2R: reuse `da` as the output buffer (no extra allocation)
    c2r.process(&mut sa, &mut da).unwrap();

    let norm = n as f32;
    da[..n_bins].iter().map(|&x| x / norm).collect()
}

// ── PMF helpers (f64 for precision across the multi-chromosome fold) ──────────

/// Build a normalised PMF from an array of real-valued f32 overlaps.
/// Returns f64 values so the polynomial multiplication stays precise.
fn build_pmf(overlaps: &[f32]) -> Vec<f64> {
    if overlaps.is_empty() {
        return vec![1.0];
    }
    let max_k = overlaps
        .iter()
        .cloned()
        .fold(0.0f32, f32::max)
        .max(0.0)
        .round() as usize;
    let mut counts = vec![0u64; max_k + 1];
    for &v in overlaps {
        let k = (v.max(0.0).round() as usize).min(max_k);
        counts[k] += 1;
    }
    let n = overlaps.len() as f64;
    counts.iter().map(|&c| c as f64 / n).collect()
}

/// Polynomial multiplication of two PMFs via f64 FFT.
fn pmf_convolve(a: &[f64], b: &[f64], planner: &mut FftPlanner<f64>) -> Vec<f64> {
    let out_len = a.len() + b.len() - 1;
    let n = out_len.next_power_of_two();

    let to_complex = |v: &[f64]| -> Vec<Complex<f64>> {
        let mut buf: Vec<Complex<f64>> = v.iter().map(|&x| Complex::new(x, 0.0)).collect();
        buf.resize(n, Complex::new(0.0, 0.0));
        buf
    };

    let mut fa = to_complex(a);
    let mut fb = to_complex(b);

    let fft = planner.plan_fft_forward(n);
    fft.process(&mut fa);
    fft.process(&mut fb);

    for (a, b) in fa.iter_mut().zip(fb.iter()) {
        *a *= b;
    }

    let ifft = planner.plan_fft_inverse(n);
    ifft.process(&mut fa);

    fa[..out_len].iter().map(|x| x.re / n as f64).collect()
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: pval <a.bed> <b.bed>");
        std::process::exit(1);
    }

    let a = parse_bed(&args[1]);
    let b = parse_bed(&args[2]);

    let n_a: usize = a.values().map(|v| v.len()).sum();
    let n_b: usize = b.values().map(|v| v.len()).sum();
    eprintln!("A: {} intervals across {} chroms", n_a, a.len());
    eprintln!("B: {} intervals across {} chroms", n_b, b.len());

    let empty: Intervals = Vec::new();

    // Per-chromosome FFT work in parallel; each task owns its RealFftPlanner.
    // Rayon preserves slice order so the fold below is deterministic.
    let chrom_results: Vec<Option<(f32, Vec<f64>)>> = hg38_chrom_sizes()
        .par_iter()
        .map(|&(chrom, size)| {
            let a_ivs = a.get(chrom).unwrap_or(&empty);
            let b_ivs = b.get(chrom).unwrap_or(&empty);
            if a_ivs.is_empty() || b_ivs.is_empty() {
                return None;
            }
            let c = chrom_overlap_vs_shift(a_ivs, b_ivs, size);
            Some((c[0], build_pmf(&c)))
        })
        .collect();

    // Sequential PMF fold (convolution is associative, order doesn't matter)
    let mut pmf_planner = FftPlanner::<f64>::new();
    let mut combined_pmf: Vec<f64> = vec![1.0];
    let mut observed = 0.0f64;
    let mut n_processed = 0usize;

    for (obs_c, pmf) in chrom_results.into_iter().flatten() {
        observed += obs_c as f64;
        combined_pmf = pmf_convolve(&combined_pmf, &pmf, &mut pmf_planner);
        n_processed += 1;
    }
    eprintln!("processed {} chromosomes", n_processed);

    // Right-tailed p-value: P(overlap >= observed)
    let obs_k = observed.round() as usize;
    let p_value: f64 = if obs_k < combined_pmf.len() {
        combined_pmf[obs_k..].iter().sum()
    } else {
        0.0
    };

    let e_overlap: f64 = combined_pmf
        .iter()
        .enumerate()
        .map(|(k, &p)| k as f64 * p)
        .sum();

    println!("---");
    println!(
        "observed overlap : {:.0} bins  ({:.0} bp)",
        observed,
        observed * BIN_SIZE as f64
    );
    println!(
        "expected overlap : {:.1} bins  ({:.0} bp)",
        e_overlap,
        e_overlap * BIN_SIZE as f64
    );
    println!("p-value (right)  : {:.4e}", p_value);
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identical_intervals_low_pvalue() {
        let ivs: Intervals = vec![(0, 1000)];
        let c = chrom_overlap_vs_shift(&ivs, &ivs, 10_000);
        // 1000 bp / 100 bp/bin = 10 bins covered
        assert!(c[0] > 8.0 && c[0] < 12.0, "c[0] = {}", c[0]);
        let pmf = build_pmf(&c);
        let obs_k = c[0].round() as usize;
        let pv: f64 = pmf[obs_k..].iter().sum();
        assert!(pv < 0.15, "p-value = {pv}");
    }

    #[test]
    fn test_no_overlap_at_zero() {
        let a: Intervals = vec![(0, 500)];
        let b: Intervals = vec![(5000, 5500)];
        let c = chrom_overlap_vs_shift(&a, &b, 10_000);
        assert!(c[0].abs() < 1.0, "expected ~0 overlap, got {}", c[0]);
    }

    #[test]
    fn test_pmf_sums_to_one() {
        let ivs: Intervals = vec![(200, 800), (1500, 2000)];
        let c = chrom_overlap_vs_shift(&ivs, &ivs, 10_000);
        let pmf = build_pmf(&c);
        let total: f64 = pmf.iter().sum();
        assert!((total - 1.0).abs() < 1e-6, "PMF sums to {total}");
    }

    #[test]
    fn test_pmf_convolve_uniform() {
        let u = vec![0.5, 0.5];
        let mut planner = FftPlanner::new();
        let conv = pmf_convolve(&u, &u, &mut planner);
        assert_eq!(conv.len(), 3);
        assert!((conv[0] - 0.25).abs() < 1e-9);
        assert!((conv[1] - 0.50).abs() < 1e-9);
        assert!((conv[2] - 0.25).abs() < 1e-9);
    }
}
