//! Fourier-based p-value for genomic interval overlap.
//!
//! # Index time
//! [`DepthMap::build`] records a per-chromosome sparse Dirac impulse train
//! (the derivative of coverage: +1 at each interval start, −1 at each end).
//! Storage is O(K) for K intervals.  Saved to `{db}/depthmap/{sid}.bin`.
//!
//! # Query time (two-phase)
//! Phase A — build DB coverage spectra once per DB file (cached):
//!   [`build_db_spectra`] converts each chromosome's impulse train to a full
//!   coverage spectrum G_B[k] via `realfft` (real→complex FFT, all N/2+1
//!   coefficients).  Coverage is built O(N) via prefix sum of the impulse
//!   train, then FFT'd in O(N log N).
//!
//! Phase B — per (query, DB) pair:
//!   1. Build query spectrum G_Q[k] (same pipeline).
//!   2. For each shared chromosome c, form the cross-correlation spectrum
//!      C_c[k] = G_Q[k] · conj(G_B[k]).  This is a pure O(N/2) dot product.
//!   3. Extract the first two moments of the null (random-shift) overlap
//!      distribution via Parseval's theorem — no IFFT required:
//!        null mean:  μ = Σ_c C_c[0].re / N_c
//!        null var:   σ² = Σ_c (1/N_c²) [2·Σ_{k=1}^{N_c/2−1} |C_c[k]|²
//!                                         + |C_c[N_c/2]|²]
//!        observed:   c_obs = Σ_c (1/N_c) [C_c[0].re
//!                                          + 2·Σ_{k=1}^{N_c/2−1} Re(C_c[k])
//!                                          + C_c[N_c/2].re]
//!   4. Gaussian right-tailed p-value: P(Z > z) where z = (c_obs − μ) / σ.
//!
//! # Gaussian null distribution
//! Under the rigid-body-shift null model the total overlap is approximately
//! Gaussian.  The moments are exact for the empirical shift distribution:
//! no sampling, no IFFT, no frequency truncation.  All N/2+1 spectral
//! coefficients contribute to σ², so the full spectrum must be retained.
//! The resulting p-value is the exact continuous-Normal right tail.
//!
//! # Null model
//! All chromosomes are shifted by the same sample index s simultaneously.
//! Variances accumulate independently per chromosome (incommensurate periods
//! make cross-chromosome covariance terms vanish asymptotically).

use std::cell::RefCell;
use std::collections::HashMap;
use std::fs;
use std::io::{self, Read, Write};
use std::path::Path;

use rayon::prelude::*;
use realfft::RealFftPlanner;
use rustfft::num_complex::Complex;

thread_local! {
    static REAL_PLANNER: RefCell<RealFftPlanner<f32>> = RefCell::new(RealFftPlanner::new());
}

const BIN_SIZE: u32 = 100;

// ── Public types ─────────────────────────────────────────────────────────────

/// Chromosome name → (start, end) interval pairs.
pub type BedMap = HashMap<String, Vec<(u32, u32)>>;

/// Precomputed coverage spectrum for one chromosome of one DB file.
/// Produced by [`build_db_spectra`]; reused across multiple queries.
pub struct ChromDbSpec {
    pub chrom: String,
    /// Padded length N = next power of two ≥ n_bins.
    pub n: usize,
    /// Number of 100 bp bins on this chromosome.
    pub n_bins: usize,
    /// V = number of impulses (= 2 × interval count on this chromosome).
    pub v: f64,
    /// G_B[k] for k = 0..N/2+1 (full coverage spectrum, all coefficients).
    pub spec: Vec<Complex<f32>>,
}

/// Sparse Dirac impulse train (derivative coverage) for one chromosome.
#[derive(Debug)]
pub struct ChromDepthMap {
    pub chrom: String,
    pub n_bins: usize,
    /// Bin indices with a +1 impulse (interval starts, floor(s/BIN)).
    pub pos_spikes: Vec<u32>,
    /// Bin indices with a -1 impulse (interval ends, ceil(e/BIN)).
    pub neg_spikes: Vec<u32>,
}

impl ChromDepthMap {
    /// Sum of absolute impulse magnitudes.
    pub fn v(&self) -> f64 {
        (self.pos_spikes.len() + self.neg_spikes.len()) as f64
    }
    /// Total coverage: Σ_j d[j]·(n_bins − j).
    pub fn total_cov(&self) -> f64 {
        let pos: f64 = self
            .pos_spikes
            .iter()
            .map(|&j| (self.n_bins - j as usize) as f64)
            .sum();
        let neg: f64 = self
            .neg_spikes
            .iter()
            .map(|&j| (self.n_bins - j as usize) as f64)
            .sum();
        pos - neg
    }
    fn n(&self) -> usize {
        self.n_bins.next_power_of_two()
    }
}

/// All per-chromosome depth maps for one indexed BED file.
#[derive(Debug)]
pub struct DepthMap {
    pub chroms: Vec<ChromDepthMap>,
}

impl DepthMap {
    /// Build from a BedMap, covering all hg38 chromosomes present in `bed`.
    pub fn build(bed: &BedMap) -> Self {
        let chroms = hg38_chrom_sizes()
            .iter()
            .filter_map(|&(chrom, size)| {
                let ivs = bed.get(chrom)?;
                if ivs.is_empty() {
                    return None;
                }
                let n_bins = ((size + BIN_SIZE - 1) / BIN_SIZE) as usize;
                let (pos, neg) = build_spikes(ivs, n_bins);
                Some(ChromDepthMap {
                    chrom: chrom.to_string(),
                    n_bins,
                    pos_spikes: pos,
                    neg_spikes: neg,
                })
            })
            .collect();
        DepthMap { chroms }
    }

    /// Save to `DMAP\x01\x00` binary format.
    pub fn save(&self, path: &Path) -> io::Result<()> {
        let mut buf = Vec::new();
        buf.extend_from_slice(b"DMAP\x01\x00");
        buf.extend_from_slice(&(self.chroms.len() as u32).to_le_bytes());
        for cdm in &self.chroms {
            let name = cdm.chrom.as_bytes();
            buf.extend_from_slice(&(name.len() as u32).to_le_bytes());
            buf.extend_from_slice(name);
            buf.extend_from_slice(&(cdm.n_bins as u64).to_le_bytes());
            write_u32_slice(&mut buf, &cdm.pos_spikes);
            write_u32_slice(&mut buf, &cdm.neg_spikes);
        }
        fs::File::create(path)?.write_all(&buf)
    }

    /// Load from a binary file produced by [`DepthMap::save`].
    pub fn load(path: &Path) -> io::Result<Self> {
        let mut buf = Vec::new();
        fs::File::open(path)?.read_to_end(&mut buf)?;
        let mut pos = 0;
        if buf.get(pos..pos + 6) != Some(b"DMAP\x01\x00") {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid DMAP magic",
            ));
        }
        pos += 6;
        let nc = read_u32(&buf, &mut pos) as usize;
        let mut chroms = Vec::with_capacity(nc);
        for _ in 0..nc {
            let nlen = read_u32(&buf, &mut pos) as usize;
            let chrom = String::from_utf8(buf[pos..pos + nlen].to_vec())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            pos += nlen;
            let n_bins = read_u64(&buf, &mut pos) as usize;
            let pos_spikes = read_u32_slice(&buf, &mut pos);
            let neg_spikes = read_u32_slice(&buf, &mut pos);
            chroms.push(ChromDepthMap {
                chrom,
                n_bins,
                pos_spikes,
                neg_spikes,
            });
        }
        Ok(DepthMap { chroms })
    }

    /// Total V across all chromosomes.
    pub fn total_v(&self) -> f64 {
        self.chroms.iter().map(|c| c.v()).sum()
    }
}

// ── Public helpers ────────────────────────────────────────────────────────────

/// V for a BedMap: 2 × total interval count (one +1 and one −1 spike per interval).
pub fn bed_map_v(bed: &BedMap) -> f64 {
    bed.values().map(|v| 2 * v.len()).sum::<usize>() as f64
}

/// Build per-chromosome coverage spectra for a DB DepthMap.
/// All N/2+1 coefficients are retained for each chromosome — the full
/// spectrum is required for an accurate Parseval variance estimate.
/// Runs in parallel over chromosomes.
pub fn build_db_spectra(dm: &DepthMap) -> Vec<ChromDbSpec> {
    dm.chroms
        .par_iter()
        .map(|cdm| {
            let spec = chrom_coverage_spectrum(cdm);
            ChromDbSpec {
                chrom: cdm.chrom.clone(),
                n: cdm.n(),
                n_bins: cdm.n_bins,
                v: cdm.v(),
                spec,
            }
        })
        .collect()
}

/// Compute the right-tailed p-value for a (query, DB) pair given pre-built
/// spectra for both sides (produced by [`build_db_spectra`]).
///
/// ### Algorithm
/// For each chromosome present in both `q_spectra` and `db_spectra`,
/// form C_c[k] = G_Q[k] · conj(G_B[k]) (O(N/2) multiplications).
/// Then, using Parseval's theorem on the full spectrum:
///
/// - **null mean** μ = Σ_c C_c[0].re / N_c
///   (expected overlap when query is uniformly shifted)
/// - **null variance** σ² = Σ_c (1/N_c²) [2·Σ_{k=1}^{N_c/2−1} |C_c[k]|²
///   + |C_c[N_c/2]|²]
/// - **observed** c_obs = Σ_c (1/N_c) [C_c[0].re
///   + 2·Σ_{k=1}^{N_c/2−1} Re(C_c[k]) + C_c[N_c/2].re]
///
/// The z-score z = (c_obs − μ) / σ is then converted to a right-tailed
/// probability via the standard-normal survival function.
///
/// Returns `(observed_bins, p_value)` or `None` if no chromosomes are shared.
pub fn compute_pvalue_cached(
    q_spectra: &[ChromDbSpec],
    db_spectra: &[ChromDbSpec],
) -> Option<(f64, f64)> {
    let mut null_mean = 0.0f64;
    let mut null_var = 0.0f64;
    let mut observed = 0.0f64;
    let mut has_data = false;

    for db_cs in db_spectra {
        let q_spec = match q_spectra.iter().find(|s| s.chrom == db_cs.chrom) {
            Some(s) => &s.spec,
            None => continue,
        };

        let len = db_cs.spec.len().min(q_spec.len()); // both should be N/2+1
        let n = db_cs.n as f64;
        let n2 = n * n;

        // k=0 (DC): both spectra are real here.
        let c0_re = (q_spec[0] * db_cs.spec[0].conj()).re as f64;
        null_mean += c0_re / n;
        observed += c0_re / n;

        // k=1..len-2 (positive frequencies, each has a conjugate mirror):
        // factor of 2 from Hermitian symmetry for both observed and variance.
        for k in 1..len.saturating_sub(1) {
            let c = q_spec[k] * db_cs.spec[k].conj();
            observed += 2.0 * c.re as f64 / n;
            null_var += 2.0 * (c.re * c.re + c.im * c.im) as f64 / n2;
        }

        // k=N/2 (Nyquist): real-valued, counted once.
        if len >= 2 {
            let c_nyq = q_spec[len - 1] * db_cs.spec[len - 1].conj();
            let sq_nyq = (c_nyq.re * c_nyq.re + c_nyq.im * c_nyq.im) as f64;
            observed += c_nyq.re as f64 / n;
            null_var += sq_nyq / n2;
        }

        has_data = true;
    }

    if !has_data {
        return None;
    }

    let p_value = if null_var == 0.0 {
        // No spread in the null → p=1 (cannot be more extreme than anything)
        1.0
    } else if null_var > null_mean {
        // Overdispersed: fit a negative binomial to the exact Parseval moments and
        // use its continuous right-tail I_p(k, r) for an exact, non-sampled p-value.
        let p = 1.0 - null_mean / null_var;
        let r = null_mean * null_mean / (null_var - null_mean);
        nb_sf(observed, r, p)
    } else {
        // Equi- or under-dispersed: fall back to the Gaussian approximation.
        normal_sf((observed - null_mean) / null_var.sqrt())
    };

    Some((observed, p_value))
}

/// Convenience: compute p-value directly from two BedMaps (builds spectra inline).
/// Use [`build_db_spectra`] + [`compute_pvalue_cached`] when either file is
/// compared against multiple partners (amortises the O(N log N) forward FFT).
pub fn compute_pvalue(query_bed: &BedMap, db_dm: &DepthMap) -> Option<(f64, f64)> {
    let db_spectra = build_db_spectra(db_dm);
    let q_dm = DepthMap::build(query_bed);
    let q_spectra = build_db_spectra(&q_dm);
    compute_pvalue_cached(&q_spectra, &db_spectra)
}

/// Convert a `parse_bed_file` result into a [`BedMap`].
pub fn intervals_to_bed_map(shards: &HashMap<String, Vec<crate::core::Interval>>) -> BedMap {
    shards
        .iter()
        .map(|(k, ivs)| (k.clone(), ivs.iter().map(|iv| (iv.start, iv.end)).collect()))
        .collect()
}

/// Parse a BED file (plain or gzip) into a [`BedMap`].
pub fn parse_bed_as_map(path: &Path) -> Result<BedMap, crate::io::BedParseError> {
    let shards = crate::io::parse_bed_file(path, 0)?;
    Ok(shards
        .into_iter()
        .map(|(k, ivs)| (k, ivs.into_iter().map(|iv| (iv.start, iv.end)).collect()))
        .collect())
}

/// hg38 chromosome sizes.
pub fn hg38_chrom_sizes() -> &'static [(&'static str, u32)] {
    &[
        ("chr1", 248_956_422),
        ("chr2", 242_193_529),
        ("chr3", 198_295_559),
        ("chr4", 190_214_555),
        ("chr5", 181_538_259),
        ("chr6", 170_805_979),
        ("chr7", 159_345_973),
        ("chr8", 145_138_636),
        ("chr9", 138_394_717),
        ("chr10", 133_797_422),
        ("chr11", 135_086_622),
        ("chr12", 133_275_309),
        ("chr13", 114_364_328),
        ("chr14", 107_043_718),
        ("chr15", 101_991_189),
        ("chr16", 90_338_345),
        ("chr17", 83_257_441),
        ("chr18", 80_373_285),
        ("chr19", 58_617_616),
        ("chr20", 64_444_167),
        ("chr21", 46_709_983),
        ("chr22", 50_818_468),
        ("chrX", 156_040_895),
        ("chrY", 57_227_415),
    ]
}

// ── Statistics helpers ────────────────────────────────────────────────────────

/// Right-tailed P(NB(r, p) ≥ k) for real-valued k via the regularized
/// incomplete beta: I_p(k, r).
///
/// This is the primary p-value for overdispersed null distributions.
/// The continuous generalization (real k) gives an "exact-as-continuous"
/// tail probability without any sampling or discretisation error.
pub fn nb_sf(k: f64, r: f64, p: f64) -> f64 {
    if k <= 0.0 {
        return 1.0;
    }
    reg_inc_beta(k, r, p).clamp(0.0, 1.0)
}

/// Regularized lower incomplete beta  I_x(a, b) = B(x;a,b)/B(a,b).
///
/// Accurate to ~1e-14 for a, b > 0 and 0 ≤ x ≤ 1, using the Lentz
/// continued-fraction algorithm (Numerical Recipes §6.4).
fn reg_inc_beta(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    // log of the prefactor  x^a (1−x)^b / (a B(a,b))
    let log_bt = a * x.ln() + b * (1.0 - x).ln() + lgamma(a + b) - lgamma(a) - lgamma(b);
    // Choose the CF that converges faster via the symmetry I_x(a,b)=1−I_{1-x}(b,a)
    if x < (a + 1.0) / (a + b + 2.0) {
        log_bt.exp() * beta_cf(a, b, x) / a
    } else {
        1.0 - log_bt.exp() * beta_cf(b, a, 1.0 - x) / b
    }
}

/// Continued-fraction kernel for [`reg_inc_beta`] (Lentz's method).
fn beta_cf(a: f64, b: f64, x: f64) -> f64 {
    const FPMIN: f64 = f64::MIN_POSITIVE;
    const EPS: f64 = 3e-15;

    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;

    let mut c = 1.0f64;
    let init = 1.0 - qab * x / qap;
    let mut d = if init.abs() < FPMIN { FPMIN } else { init };
    d = 1.0 / d;
    let mut h = d;

    for m in 1_u32..=200 {
        let mf = m as f64;
        let m2 = 2.0 * mf;

        // Even step
        let aa = mf * (b - mf) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        h *= d * c;

        // Odd step
        let aa = -(a + mf) * (qab + mf) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() <= EPS {
            break;
        }
    }
    h
}

/// Natural log of the gamma function, accurate to ~15 significant figures.
///
/// Uses the Lanczos approximation (g=7, 9 coefficients from Godfrey 2001)
/// with the reflection formula for x < 0.5.
fn lgamma(x: f64) -> f64 {
    if x < 0.5 {
        // Γ(x)Γ(1−x) = π/sin(πx)  ⟹  ln Γ(x) = ln π − ln sin(πx) − ln Γ(1−x)
        std::f64::consts::PI.ln()
            - (std::f64::consts::PI * x).sin().ln()
            - lgamma(1.0 - x)
    } else {
        const G: f64 = 7.0;
        #[rustfmt::skip]
        const C: [f64; 9] = [
            0.999_999_999_999_809_93,
            676.520_368_121_885_10,
           -1259.139_216_722_402_8,
            771.323_428_777_653_13,
           -176.615_029_162_140_59,
             12.507_343_278_686_905,
             -0.138_571_095_265_720_12,
              9.984_369_578_019_571_6e-6,
              1.505_632_735_149_311_6e-7,
        ];
        let z = x - 1.0;
        let mut s = C[0];
        for (i, &ci) in C[1..].iter().enumerate() {
            s += ci / (z + (i + 1) as f64);
        }
        let t = z + G + 0.5;
        0.5 * std::f64::consts::TAU.ln() + (z + 0.5) * t.ln() - t + s.ln()
    }
}

/// Right-tailed probability P(Z > z) for a standard normal Z.
/// Used as a fallback when the null distribution is equi- or under-dispersed.
pub fn normal_sf(z: f64) -> f64 {
    (erfc(z * std::f64::consts::FRAC_1_SQRT_2) * 0.5).clamp(0.0, 1.0)
}

/// Complementary error function, accurate to ~1 ulp for all finite x.
///
/// Uses a Taylor series for |x| < 1.5 and the Lentz continued-fraction
/// algorithm for the Laplace CF otherwise:
///   erfc(x) = exp(−x²)/√π / (x + (1/2)/(x + 1/(x + (3/2)/(x + …))))
fn erfc(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc(-x);
    }
    if x > 27.0 {
        return 0.0; // exp(−x²) underflows before this point
    }

    if x < 1.5 {
        // erf via Taylor series, then erfc = 1 − erf.
        const FRAC_2_SQRT_PI: f64 = 1.128_379_167_095_512_6;
        let x2 = x * x;
        let mut term = x;
        let mut sum = x;
        for n in 1_u32..=60 {
            term *= -x2 / n as f64;
            let delta = term / (2 * n + 1) as f64;
            sum += delta;
            if delta.abs() <= sum.abs() * f64::EPSILON {
                break;
            }
        }
        (1.0 - FRAC_2_SQRT_PI * sum).max(0.0)
    } else {
        // Lentz algorithm for h = x + (1/2)/(x + 1/(x + (3/2)/(x + …))),
        // where erfc(x) = exp(−x²) / (√π · h).
        const INV_SQRT_PI: f64 = 0.564_189_583_547_756_3;
        let tiny = f64::MIN_POSITIVE;
        let mut f = x;
        let mut c = x;
        let mut d = 0.0f64;
        for i in 1_u32..=120 {
            let a = i as f64 * 0.5;
            d = x + a * d;
            if d == 0.0 {
                d = tiny;
            }
            c = x + a / c;
            if c == 0.0 {
                c = tiny;
            }
            d = 1.0 / d;
            let delta = c * d;
            f *= delta;
            if (delta - 1.0).abs() <= f64::EPSILON {
                break;
            }
        }
        (INV_SQRT_PI * (-x * x).exp() / f).max(0.0)
    }
}

// ── Private helpers ───────────────────────────────────────────────────────────

/// Convert intervals to (+1, −1) spike index vectors.
fn build_spikes(intervals: &[(u32, u32)], n_bins: usize) -> (Vec<u32>, Vec<u32>) {
    let mut pos = Vec::with_capacity(intervals.len());
    let mut neg = Vec::with_capacity(intervals.len());
    for &(s, e) in intervals {
        let sb = ((s / BIN_SIZE) as usize).min(n_bins) as u32;
        let eb = (((e + BIN_SIZE - 1) / BIN_SIZE) as usize).min(n_bins) as u32;
        pos.push(sb);
        neg.push(eb);
    }
    (pos, neg)
}

/// Build a real-valued coverage array of length `n = cdm.n()` via prefix sum
/// of the impulse train, then forward-FFT it with `realfft`, returning all
/// N/2+1 complex coefficients.
///
/// Cost: O(N) build + O(N log N) FFT.
fn chrom_coverage_spectrum(cdm: &ChromDepthMap) -> Vec<Complex<f32>> {
    let n = cdm.n();

    let fft = REAL_PLANNER.with(|p| p.borrow_mut().plan_fft_forward(n));
    let mut indata = fft.make_input_vec(); // vec![0f32; n]

    // Build impulse train (derivative of coverage) in-place.
    for &j in &cdm.pos_spikes {
        let idx = (j as usize).min(n - 1);
        indata[idx] += 1.0;
    }
    for &j in &cdm.neg_spikes {
        let idx = (j as usize).min(n - 1);
        indata[idx] -= 1.0;
    }
    // Prefix sum: impulse train → coverage g[j].
    for i in 1..n {
        indata[i] += indata[i - 1];
    }

    let mut spectrum = fft.make_output_vec(); // complex, length n/2+1
    fft.process(&mut indata, &mut spectrum).unwrap();

    spectrum
}

// ── Binary I/O helpers ────────────────────────────────────────────────────────

fn write_u32_slice(buf: &mut Vec<u8>, v: &[u32]) {
    buf.extend_from_slice(&(v.len() as u64).to_le_bytes());
    for &x in v {
        buf.extend_from_slice(&x.to_le_bytes());
    }
}

fn read_u32(buf: &[u8], pos: &mut usize) -> u32 {
    let v = u32::from_le_bytes(buf[*pos..*pos + 4].try_into().unwrap());
    *pos += 4;
    v
}

fn read_u64(buf: &[u8], pos: &mut usize) -> u64 {
    let v = u64::from_le_bytes(buf[*pos..*pos + 8].try_into().unwrap());
    *pos += 8;
    v
}

fn read_u32_slice(buf: &[u8], pos: &mut usize) -> Vec<u32> {
    let len = read_u64(buf, pos) as usize;
    let mut v = Vec::with_capacity(len);
    for _ in 0..len {
        v.push(read_u32(buf, pos));
    }
    v
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn small_bed() -> BedMap {
        let mut m = BedMap::new();
        m.insert("chr22".to_string(), vec![(10_000_000, 10_001_000)]);
        m
    }

    #[test]
    fn test_depth_map_roundtrip() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bin");

        let bed = small_bed();
        let dm = DepthMap::build(&bed);
        assert_eq!(dm.chroms.len(), 1);
        assert_eq!(dm.chroms[0].chrom, "chr22");

        dm.save(&path).unwrap();
        let loaded = DepthMap::load(&path).unwrap();

        assert_eq!(loaded.chroms.len(), 1);
        let orig = &dm.chroms[0];
        let back = &loaded.chroms[0];
        assert_eq!(back.chrom, orig.chrom);
        assert_eq!(back.n_bins, orig.n_bins);
        assert_eq!(back.pos_spikes, orig.pos_spikes);
        assert_eq!(back.neg_spikes, orig.neg_spikes);
    }

    #[test]
    fn test_compute_pvalue_identical() {
        let mut bed = BedMap::new();
        // 10 Mb interval = 100_000 bins on chr22 (n=524288).
        // Querying with the same interval as the DB: observed overlap = 100_000 bins.
        // Under random shifts, overlap is at its maximum only near shift=0, so the
        // Gaussian p-value should be very small.
        bed.insert("chr22".to_string(), vec![(10_000_000, 20_000_000)]);

        let dm = DepthMap::build(&bed);
        let (observed, p_value) = compute_pvalue(&bed, &dm).unwrap();
        assert!(observed > 50_000.0, "observed={observed}");
        assert!(p_value < 0.1, "p_value={p_value}");
    }

    #[test]
    fn test_compute_pvalue_no_shared_chrom() {
        let mut query = BedMap::new();
        query.insert("chr21".to_string(), vec![(10_000_000, 10_001_000)]);

        let mut db_bed = BedMap::new();
        db_bed.insert("chr22".to_string(), vec![(10_000_000, 10_001_000)]);
        let dm = DepthMap::build(&db_bed);

        assert!(compute_pvalue(&query, &dm).is_none());
    }

    #[test]
    fn test_lgamma_values() {
        // Γ(1) = Γ(2) = 1  →  ln Γ = 0
        assert!(lgamma(1.0).abs() < 1e-13);
        assert!(lgamma(2.0).abs() < 1e-13);
        // Γ(0.5) = √π  →  ln Γ = 0.5 ln π
        let expected = 0.5 * std::f64::consts::PI.ln();
        assert!((lgamma(0.5) - expected).abs() < 1e-13, "lgamma(0.5) = {}", lgamma(0.5));
        // Γ(5) = 4! = 24
        assert!((lgamma(5.0) - 24.0_f64.ln()).abs() < 1e-13);
    }

    #[test]
    fn test_reg_inc_beta_known() {
        // I_{0.5}(1, 3) = 0.875   [verified by hand above]
        let v = reg_inc_beta(1.0, 3.0, 0.5);
        assert!((v - 0.875).abs() < 1e-12, "I_0.5(1,3) = {v}");
        // I_{0.5}(3, 3) = 0.5     [symmetric beta, x=0.5]
        let v = reg_inc_beta(3.0, 3.0, 0.5);
        assert!((v - 0.5).abs() < 1e-12, "I_0.5(3,3) = {v}");
        // Boundary values
        assert!(reg_inc_beta(2.0, 3.0, 0.0) == 0.0);
        assert!(reg_inc_beta(2.0, 3.0, 1.0) == 1.0);
    }

    #[test]
    fn test_nb_sf_known() {
        // NB(r=3, p=0.5): P(X ≥ 0) = 1, P(X ≥ 3) = 0.5
        assert!((nb_sf(0.0, 3.0, 0.5) - 1.0).abs() < 1e-14);
        assert!((nb_sf(3.0, 3.0, 0.5) - 0.5).abs() < 1e-12, "nb_sf(3,3,0.5) = {}", nb_sf(3.0, 3.0, 0.5));
        // NB(r=3, p=0.5): P(X ≥ 1) = 1 - P(X=0) = 1 - (0.5)^3 = 0.875
        assert!((nb_sf(1.0, 3.0, 0.5) - 0.875).abs() < 1e-12);
        // Monotone: more extreme k → smaller p
        let p1 = nb_sf(5.0, 2.0, 0.3);
        let p2 = nb_sf(10.0, 2.0, 0.3);
        assert!(p1 > p2, "nb_sf should decrease with k");
    }

    #[test]
    fn test_normal_sf_symmetry() {
        // P(Z > 0) = 0.5
        assert!((normal_sf(0.0) - 0.5).abs() < 1e-15);
        // P(Z > z) + P(Z > -z) = 1
        for z in [1.0, 2.5, 5.0, 10.0] {
            let sum = normal_sf(z) + normal_sf(-z);
            assert!((sum - 1.0).abs() < 1e-12, "z={z}: sum={sum}");
        }
    }

    #[test]
    fn test_cached_matches_direct() {
        let mut bed = BedMap::new();
        bed.insert(
            "chr22".to_string(),
            vec![(10_000_000, 10_001_000), (20_000_000, 20_002_000)],
        );
        let dm = DepthMap::build(&bed);
        let db_spectra = build_db_spectra(&dm);
        let q_dm = DepthMap::build(&bed);
        let q_spectra = build_db_spectra(&q_dm);

        let (obs_cached, pv_cached) = compute_pvalue_cached(&q_spectra, &db_spectra).unwrap();
        let (obs_direct, pv_direct) = compute_pvalue(&bed, &dm).unwrap();

        assert!((obs_cached - obs_direct).abs() < 1.0, "obs mismatch");
        assert!((pv_cached - pv_direct).abs() < 1e-10, "pv mismatch");
    }
}
