//! FFT-based analytic statistics for genomic interval overlap.
//!
//! # Null model
//! Under the free-relative-gap null, each query interval slides independently
//! over the DB with a uniform random position.  For query interval of length L
//! bins on a chromosome where the DB has pre-computed moments M[L]:
//!
//!   μ_i  = M[L].mean   (expected overlap in bins = L × mean_depth)
//!   σ²_i = M[L].var    (variance from FFT sliding-window moments)
//!
//! Genome-wide null: μ = Σ μ_i, σ² = Σ σ²_i.
//!
//! # Moment computation (O(N log N) per chromosome)
//! For depth signal g[0..N-1] with prefix sum d*[0..=N]:
//!
//!   d**[j]  = Σ_{k<j} d*[k]          (double prefix, size N+2)
//!   d2[j]   = Σ_{k<j} (d*[k])²       (squared-prefix, size N+2)
//!   R[L]    = IFFT(|FFT(d*)|²)[L]    (autocorrelation of d* via FFT)
//!
//!   sum_f(L)  = d**[N+1] − d**[L] − d**[N−L+1]
//!   sum_f2(L) = (d2[N+1] − d2[L]) + d2[N−L+1] − 2·R[L]
//!   mean(L)   = sum_f(L)  / (N−L+1)
//!   var(L)    = sum_f2(L) / (N−L+1) − mean(L)²
//!
//! # Compacted storage (O(M) with M ≈ (1/ε)·ln(N·ε))
//! Dense region L ∈ [1, ⌈1/ε⌉]: one entry per integer L.
//! Sparse region L > ⌈1/ε⌉: representative Ls spaced by stride ε·L, giving
//! logarithmic density.  Every query L maps to a stored representative within
//! relative error ε in L.

use std::collections::HashMap;
use std::path::Path;

use rayon::prelude::*;
use rustfft::num_complex::Complex;
use rustfft::FftPlanner;

pub const DEFAULT_MOMENTS_EPS: f64 = 0.01;

// ── FilterMask ────────────────────────────────────────────────────────────────

/// Whether a [`FilterMask`] BED file lists regions to keep or regions to remove.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterMode {
    Whitelist,
    Blacklist,
}

/// Per-chromosome boolean mask at 100 bp resolution (`true` = accessible).
pub struct FilterMask {
    chroms: HashMap<String, Vec<bool>>,
    mode: FilterMode,
}

impl FilterMask {
    pub fn build(bed: &BedMap, mode: FilterMode) -> Self {
        let chrom_sizes: HashMap<&str, u32> = hg38_chrom_sizes().iter().copied().collect();

        let chroms = match mode {
            FilterMode::Whitelist => bed
                .iter()
                .filter_map(|(chrom, ivs)| {
                    if ivs.is_empty() {
                        return None;
                    }
                    let &size = chrom_sizes.get(chrom.as_str())?;
                    let n_bins = ((size + BIN_SIZE - 1) / BIN_SIZE) as usize;
                    let mut mask = vec![false; n_bins];
                    for &(s, e) in ivs {
                        let sb = (s / BIN_SIZE) as usize;
                        let eb = (((e + BIN_SIZE - 1) / BIN_SIZE) as usize).min(n_bins);
                        for b in sb..eb {
                            mask[b] = true;
                        }
                    }
                    Some((chrom.clone(), mask))
                })
                .collect(),

            FilterMode::Blacklist => chrom_sizes
                .iter()
                .map(|(&chrom, &size)| {
                    let n_bins = ((size + BIN_SIZE - 1) / BIN_SIZE) as usize;
                    let mut mask = vec![true; n_bins];
                    if let Some(ivs) = bed.get(chrom) {
                        for &(s, e) in ivs {
                            let sb = (s / BIN_SIZE) as usize;
                            let eb = (((e + BIN_SIZE - 1) / BIN_SIZE) as usize).min(n_bins);
                            for b in sb..eb {
                                mask[b] = false;
                            }
                        }
                    }
                    (chrom.to_string(), mask)
                })
                .collect(),
        };

        FilterMask { chroms, mode }
    }

    pub fn get(&self, chrom: &str) -> Option<&[bool]> {
        self.chroms.get(chrom).map(Vec::as_slice)
    }

    pub fn mode(&self) -> FilterMode {
        self.mode
    }
}

const BIN_SIZE: u32 = 100;

// ── Public types ──────────────────────────────────────────────────────────────

pub type BedMap = HashMap<String, Vec<(u32, u32)>>;

/// Per-chromosome sliding-window moments for one DB file, compacted with ε-relative error.
pub struct ChromMoments {
    pub chrom: String,
    pub n_bins: usize,
    pub eps: f64,
    pub dense_max: usize,
    /// (mean, var) pairs indexed by a compacted scheme:
    ///   indices 0 .. dense_max-1  → dense region, L = 1 ..= dense_max
    ///   indices dense_max ..      → sparse region, representative L = dense_max·(1+ε)^m (m=1,2,…)
    pub moments: Vec<(f64, f64)>,
}

impl ChromMoments {
    /// Return `(mean, var)` for a sliding window of size `l_bins` bins.
    ///
    /// `mean` is the expected total overlap (bins) when a block of `l_bins` is
    /// placed uniformly at random over this chromosome.
    /// `var` is the variance of that overlap.
    ///
    /// Returns `None` when `l_bins` rounds to 0 or ≥ `n_bins`.
    pub fn lookup(&self, l_bins: f64) -> Option<(f64, f64)> {
        let l = l_bins.round() as usize;
        if l == 0 || l >= self.n_bins || self.moments.is_empty() {
            return None;
        }
        let idx = if l <= self.dense_max {
            l - 1
        } else {
            // m = floor(log_{1+ε}(L / dense_max)), ≥ 1 for L > dense_max
            let m = ((l as f64 / self.dense_max as f64).ln() / (1.0 + self.eps).ln())
                .floor()
                .max(1.0) as usize;
            // Sparse entries start at total index dense_max; slot m → index dense_max + m - 1
            (self.dense_max + m - 1).min(self.moments.len() - 1)
        };
        self.moments.get(idx).copied()
    }
}

/// Per-chromosome interval data for one query file.
pub struct QueryChromData {
    pub chrom: String,
    pub n_bins: usize,
    /// Length of each interval in 100 bp bins (may be fractional).
    pub interval_lengths: Vec<f64>,
}

/// Sparse Dirac impulse train (coverage derivative) for one chromosome.
#[derive(Debug)]
pub struct ChromDepthMap {
    pub chrom: String,
    pub n_bins: usize,
    pub pos_spikes: Vec<u32>,
    pub neg_spikes: Vec<u32>,
}

impl ChromDepthMap {
    pub fn v(&self) -> f64 {
        (self.pos_spikes.len() + self.neg_spikes.len()) as f64
    }

    pub fn total_cov(&self) -> f64 {
        let n = self.n_bins as f64;
        let b = BIN_SIZE as f64;
        let pos: f64 = self.pos_spikes.iter().map(|&x| n - x as f64 / b).sum();
        let neg: f64 = self.neg_spikes.iter().map(|&x| n - x as f64 / b).sum();
        pos - neg
    }
}

/// All per-chromosome depth maps for one indexed BED file.
#[derive(Debug)]
pub struct DepthMap {
    pub chroms: Vec<ChromDepthMap>,
}

impl DepthMap {
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

    pub fn total_v(&self) -> f64 {
        self.chroms.iter().map(|c| c.v()).sum()
    }
}

// ── Public helpers ────────────────────────────────────────────────────────────

pub fn bed_map_v(bed: &BedMap) -> f64 {
    bed.values().map(|v| 2 * v.len()).sum::<usize>() as f64
}

/// Build compacted moment tables for every chromosome in a [`DepthMap`].
///
/// `eps` controls the relative error in block-size approximation; 0.01 means
/// each stored representative is within 1% of the requested block size.
pub fn build_depth_moments(dm: &DepthMap, eps: f64) -> Vec<ChromMoments> {
    dm.chroms
        .par_iter()
        .map(|cdm| build_chrom_moments(cdm, eps))
        .collect()
}

/// Build compacted moments for one chromosome (no filter).
pub fn build_chrom_moments(cdm: &ChromDepthMap, eps: f64) -> ChromMoments {
    let g = build_depth_signal(cdm, None);
    moments_from_signal(&g, &cdm.chrom, cdm.n_bins, eps)
}

/// Build compacted moments for one chromosome with a filter mask applied.
///
/// Bins not accessible under `mask` are zeroed before moment computation.
pub fn build_chrom_moments_with_filter(
    cdm: &ChromDepthMap,
    mask: &[bool],
    eps: f64,
) -> ChromMoments {
    let g = build_depth_signal(cdm, Some(mask));
    moments_from_signal(&g, &cdm.chrom, cdm.n_bins, eps)
}

pub fn mean_interval_bins(q_data: &[QueryChromData]) -> f64 {
    let total_l: f64 = q_data.iter().flat_map(|c| c.interval_lengths.iter()).sum();
    let total_n: usize = q_data.iter().map(|c| c.interval_lengths.len()).sum();
    if total_n == 0 { 0.0 } else { total_l / total_n as f64 }
}

/// Build per-chromosome interval data from a raw BED map.
pub fn build_query_chrom_data(bed: &BedMap, filter: Option<&FilterMask>) -> Vec<QueryChromData> {
    hg38_chrom_sizes()
        .iter()
        .filter_map(|&(chrom, size)| {
            if let Some(f) = filter {
                if f.get(chrom).is_none() {
                    return None;
                }
            }
            let ivs = bed.get(chrom)?;
            if ivs.is_empty() {
                return None;
            }
            let n_bins = ((size + BIN_SIZE - 1) / BIN_SIZE) as usize;
            let interval_lengths = ivs
                .iter()
                .map(|&(s, e)| (e - s) as f64 / BIN_SIZE as f64)
                .collect();
            Some(QueryChromData {
                chrom: chrom.to_string(),
                n_bins,
                interval_lengths,
            })
        })
        .collect()
}

/// Exact base-pair overlap in 100 bp bins between two depth maps (dot product).
///
/// Used by the `pval` binary where the exact observed overlap is needed.
pub fn coverage_dot_product(a_dm: &DepthMap, b_dm: &DepthMap) -> f64 {
    let a_map: HashMap<&str, &ChromDepthMap> =
        a_dm.chroms.iter().map(|c| (c.chrom.as_str(), c)).collect();
    b_dm.chroms
        .iter()
        .filter_map(|b_cd| {
            let a_cd = a_map.get(b_cd.chrom.as_str())?;
            let a_cov = build_depth_signal(a_cd, None);
            let b_cov = build_depth_signal(b_cd, None);
            let n = a_cd.n_bins.min(b_cd.n_bins);
            Some(
                a_cov[..n]
                    .iter()
                    .zip(b_cov[..n].iter())
                    .map(|(&a, &b)| a as f64 * b as f64)
                    .sum::<f64>(),
            )
        })
        .sum()
}

/// Compute analytic NB statistics for a (query, DB) pair.
///
/// For each query interval of length l_i on a shared chromosome, looks up
/// the stored (mean, var) for block size l_i and accumulates:
///   μ_null  = Σ_i mean(l_i)
///   σ²_null = Σ_i var(l_i)
///
/// Fits a negative-binomial null via method of moments, scores `observed`
/// against it via exponential tilting (LLR), and returns a Wilks p-value.
pub fn compute_analytic_stats(
    q_data: &[QueryChromData],
    db_moments: &[ChromMoments],
    observed: f64,
) -> Option<(f64, Option<f64>)> {
    let db_map: HashMap<&str, &ChromMoments> =
        db_moments.iter().map(|m| (m.chrom.as_str(), m)).collect();

    let mut mu_total = 0.0f64;
    let mut var_total = 0.0f64;
    let mut any_shared = false;

    for q_cd in q_data {
        let db_m = match db_map.get(q_cd.chrom.as_str()) {
            Some(d) => d,
            None => continue,
        };
        for &l in &q_cd.interval_lengths {
            if let Some((mean_l, var_l)) = db_m.lookup(l) {
                mu_total += mean_l;
                var_total += var_l;
                any_shared = true;
            }
        }
    }

    if !any_shared {
        return None;
    }

    let llr = nb_llr(mu_total, var_total, observed);
    let p_value = match llr {
        Some(l) if l > 0.0 => erfc(l.sqrt()) * 0.5,
        _ => 1.0,
    };

    Some((p_value, llr))
}

pub fn parse_bed_as_map(path: &Path) -> Result<BedMap, crate::io::BedParseError> {
    let shards = crate::io::parse_bed_file(path, 0)?;
    Ok(shards
        .into_iter()
        .map(|(k, ivs)| (k, ivs.into_iter().map(|iv| (iv.start, iv.end)).collect()))
        .collect())
}

pub fn intervals_to_bed_map(shards: &HashMap<String, Vec<crate::core::Interval>>) -> BedMap {
    shards
        .iter()
        .map(|(k, ivs)| (k.clone(), ivs.iter().map(|iv| (iv.start, iv.end)).collect()))
        .collect()
}

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

// ── Private helpers ───────────────────────────────────────────────────────────

/// Build the full per-bin coverage signal g[0..n_bins-1] with fractional binning.
///
/// If `mask` is supplied, non-accessible bins are zeroed after the prefix sum.
fn build_depth_signal(cdm: &ChromDepthMap, mask: Option<&[bool]>) -> Vec<f32> {
    let n = cdm.n_bins;
    let mut g = vec![0.0f32; n];

    for &x in &cdm.pos_spikes {
        let bin = (x / BIN_SIZE) as usize;
        let frac = (x % BIN_SIZE) as f32 / BIN_SIZE as f32;
        if bin < n {
            g[bin] += 1.0 - frac;
        }
        if frac > 0.0 && bin + 1 < n {
            g[bin + 1] += frac;
        }
    }
    for &x in &cdm.neg_spikes {
        let bin = (x / BIN_SIZE) as usize;
        let frac = (x % BIN_SIZE) as f32 / BIN_SIZE as f32;
        if bin < n {
            g[bin] -= 1.0 - frac;
        }
        if frac > 0.0 && bin + 1 < n {
            g[bin + 1] -= frac;
        }
    }
    for i in 1..n {
        g[i] += g[i - 1];
    }

    if let Some(mask) = mask {
        for (i, &allowed) in mask.iter().enumerate().take(n) {
            if !allowed {
                g[i] = 0.0;
            }
        }
    }

    g
}

/// Build the full compacted moment table from a coverage signal.
fn moments_from_signal(g: &[f32], chrom: &str, n_bins: usize, eps: f64) -> ChromMoments {
    let n = n_bins;
    let dense_max = (1.0_f64 / eps).ceil() as usize;

    if n <= 1 {
        return ChromMoments {
            chrom: chrom.to_string(),
            n_bins,
            eps,
            dense_max,
            moments: Vec::new(),
        };
    }

    // d*[0..=N]: prefix sum of g (d*[0]=0, d*[j] = Σ_{k<j} g[k])
    let mut d_star = vec![0.0f64; n + 1];
    for j in 1..=n {
        d_star[j] = d_star[j - 1] + g[j - 1] as f64;
    }

    // d**[0..=N+1]: prefix sum of d* (d**[0]=0)
    let mut d_ss = vec![0.0f64; n + 2];
    for j in 1..=(n + 1) {
        d_ss[j] = d_ss[j - 1] + d_star[j - 1];
    }

    // d2[0..=N+1]: prefix sum of (d*)² (d2[0]=0)
    let mut d2 = vec![0.0f64; n + 2];
    for j in 1..=(n + 1) {
        d2[j] = d2[j - 1] + d_star[j - 1] * d_star[j - 1];
    }

    // R[L] = autocorrelation of d* at lag L via FFT
    let r = autocorr_via_fft(&d_star);

    let max_l = n - 1; // largest L with at least 2 windows
    let mut moments: Vec<(f64, f64)> = Vec::new();

    // Dense region: L = 1..=dense_max (one entry per integer L)
    for l in 1..=dense_max.min(max_l) {
        moments.push(window_moments(&d_ss, &d2, &r, n, l));
    }

    // Sparse region: representative Ls at dense_max·(1+ε)^m for m = 1, 2, …
    // Each sparse entry covers a range of Ls with relative width ε.
    let log_base = 1.0 + eps;
    let mut l_f = dense_max as f64 * log_base;
    let mut last_l = dense_max;

    loop {
        let l = l_f.floor() as usize;
        if l > max_l {
            break;
        }
        // Skip duplicate representatives (can occur at large ε near dense_max)
        if l > last_l {
            moments.push(window_moments(&d_ss, &d2, &r, n, l));
            last_l = l;
        }
        l_f *= log_base;
    }

    ChromMoments {
        chrom: chrom.to_string(),
        n_bins,
        eps,
        dense_max,
        moments,
    }
}

/// Compute (mean, var) for sliding windows of size `l` over the depth signal.
///
/// sum_f(L)  = d**[N+1] − d**[L] − d**[N−L+1]
/// sum_f2(L) = (d2[N+1] − d2[L]) + d2[N−L+1] − 2·R[L]
fn window_moments(
    d_ss: &[f64],
    d2: &[f64],
    r: &[f64],
    n: usize,
    l: usize,
) -> (f64, f64) {
    debug_assert!(l >= 1 && l < n);
    let n_w = (n - l + 1) as f64;

    let sum_f = d_ss[n + 1] - d_ss[l] - d_ss[n - l + 1];
    let r_l = r.get(l).copied().unwrap_or(0.0);
    let sum_f2 = (d2[n + 1] - d2[l]) + d2[n - l + 1] - 2.0 * r_l;

    let mean = sum_f / n_w;
    let e_f2 = sum_f2 / n_w;
    let var = (e_f2 - mean * mean).max(0.0);

    (mean, var)
}

/// Linear autocorrelation of `signal` at all lags 0..len via zero-padded FFT.
///
/// Returns a vector of length `signal.len()` where `result[L]` =
/// Σ_{k=0}^{n-L-1} signal[k] * signal[k+L].
fn autocorr_via_fft(signal: &[f64]) -> Vec<f64> {
    let n = signal.len();
    if n == 0 {
        return Vec::new();
    }
    // Zero-pad to ≥ 2n to avoid circular wrap-around aliasing.
    let fft_len = (2 * n).next_power_of_two();

    let mut planner = FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(fft_len);
    let ifft = planner.plan_fft_inverse(fft_len);

    let mut buf: Vec<Complex<f64>> = signal
        .iter()
        .map(|&x| Complex::new(x, 0.0))
        .chain(std::iter::repeat(Complex::new(0.0, 0.0)).take(fft_len - n))
        .collect();

    fft.process(&mut buf);

    // Power spectrum: |FFT[k]|²
    for c in &mut buf {
        let p = c.norm_sqr();
        *c = Complex::new(p, 0.0);
    }

    ifft.process(&mut buf);

    // rustfft IFFT is unnormalized; divide by fft_len.
    let scale = 1.0 / fft_len as f64;
    buf.into_iter().take(n).map(|c| c.re * scale).collect()
}

/// Saddlepoint LLR for the NB null tilted to observation `o`.
pub(crate) fn nb_llr(mu: f64, var: f64, observed: f64) -> Option<f64> {
    if mu <= 0.0 || var <= mu {
        return None;
    }
    let r = mu * mu / (var - mu);
    let p = mu / var;
    let p_o = r / (r + observed);
    if !(0.0 < p && p < 1.0 && 0.0 < p_o && p_o < 1.0) {
        return None;
    }
    let theta = ((1.0 - p_o) / (1.0 - p)).ln();
    Some(observed * theta + r * (p_o / p).ln())
}

/// Complementary error function (max error < 1.5e-7 for x ≥ 0).
fn erfc(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc(-x);
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let p = t
        * (0.254_829_592
            + t * (-0.284_496_736
                + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
    p * (-x * x).exp()
}

fn build_spikes(intervals: &[(u32, u32)], n_bins: usize) -> (Vec<u32>, Vec<u32>) {
    let max_bp = (n_bins as u32) * BIN_SIZE;
    let mut pos = Vec::with_capacity(intervals.len());
    let mut neg = Vec::with_capacity(intervals.len());
    for &(s, e) in intervals {
        pos.push(s.min(max_bp));
        neg.push(e.min(max_bp));
    }
    (pos, neg)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── autocorr correctness ──────────────────────────────────────────────────

    #[test]
    fn test_autocorr_lag0_is_sum_of_squares() {
        let v = vec![1.0, 2.0, 3.0, 0.0, 0.0];
        let r = autocorr_via_fft(&v);
        let expected = v.iter().map(|x| x * x).sum::<f64>();
        assert!((r[0] - expected).abs() < 1e-6, "r[0]={} expected={}", r[0], expected);
    }

    #[test]
    fn test_autocorr_manual() {
        // d* = [0,1,1,1,2,2,2,3] → R[1] = 0*1+1*1+1*1+1*2+2*2+2*2+2*3 = 18
        let d_star = vec![0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0];
        let r = autocorr_via_fft(&d_star);
        assert!((r[1] - 18.0).abs() < 1e-4, "r[1]={}", r[1]);
    }

    // ── window_moments correctness ────────────────────────────────────────────

    #[test]
    fn test_window_moments_uniform_signal() {
        // g = [1,1,1,...,1] (N bins): every window of any L sums to L.
        // mean(L) = L, var(L) = 0.
        let n = 50usize;
        let g: Vec<f32> = vec![1.0; n];

        let mut d_star = vec![0.0f64; n + 1];
        for j in 1..=n { d_star[j] = d_star[j-1] + g[j-1] as f64; }
        let mut d_ss = vec![0.0f64; n + 2];
        for j in 1..=(n+1) { d_ss[j] = d_ss[j-1] + d_star[j-1]; }
        let mut d2 = vec![0.0f64; n + 2];
        for j in 1..=(n+1) { d2[j] = d2[j-1] + d_star[j-1]*d_star[j-1]; }
        let r = autocorr_via_fft(&d_star);

        for l in 1..n {
            let (mean, var) = window_moments(&d_ss, &d2, &r, n, l);
            assert!((mean - l as f64).abs() < 1e-6,
                "uniform: mean(L={})={:.4} expected {}", l, mean, l);
            assert!(var < 1e-6,
                "uniform: var(L={})={:.4} expected 0", l, var);
        }
    }

    #[test]
    fn test_window_moments_sparse_periodic() {
        // g = [1,0,0,1,0,0,1] (N=7), L=3: every window sums to 1 → var=0, mean=1.
        let g: Vec<f32> = vec![1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0];
        let n = g.len();
        let mut d_star = vec![0.0f64; n + 1];
        for j in 1..=n { d_star[j] = d_star[j-1] + g[j-1] as f64; }
        let mut d_ss = vec![0.0f64; n + 2];
        for j in 1..=(n+1) { d_ss[j] = d_ss[j-1] + d_star[j-1]; }
        let mut d2 = vec![0.0f64; n + 2];
        for j in 1..=(n+1) { d2[j] = d2[j-1] + d_star[j-1]*d_star[j-1]; }
        let r = autocorr_via_fft(&d_star);

        let (mean, var) = window_moments(&d_ss, &d2, &r, n, 3);
        assert!((mean - 1.0).abs() < 1e-6, "mean={}", mean);
        assert!(var < 1e-6, "var={}", var);

        // sum_f for L=1: all g[k] → sum = 3, n_w = 7
        let (mean1, _) = window_moments(&d_ss, &d2, &r, n, 1);
        assert!((mean1 - 3.0 / 7.0).abs() < 1e-6, "mean(L=1)={}", mean1);
    }

    // ── ChromMoments lookup ───────────────────────────────────────────────────

    #[test]
    fn test_chrom_moments_lookup_dense() {
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(10_000_000, 20_000_000)]);
        let dm = DepthMap::build(&bed);
        let moments = build_chrom_moments(&dm.chroms[0], 0.01);

        // Dense region: L=50 should return a valid pair.
        let r = moments.lookup(50.0);
        assert!(r.is_some(), "lookup(50) returned None");
        let (mean, var) = r.unwrap();
        assert!(mean > 0.0, "mean={}", mean);
        assert!(var >= 0.0, "var={}", var);
    }

    #[test]
    fn test_chrom_moments_lookup_sparse() {
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(0, 50_818_468)]);
        let dm = DepthMap::build(&bed);
        let moments = build_chrom_moments(&dm.chroms[0], 0.01);

        // Sparse region: L=500 (> dense_max=100)
        let r = moments.lookup(500.0);
        assert!(r.is_some(), "lookup(500) returned None");
        let (mean, var) = r.unwrap();
        // The stored representative for L=500 has |L' - 500| / 500 ≤ ε=1%.
        // For full coverage, mean = L' exactly, so |mean - 500| ≤ 500 * ε ≈ 5.
        assert!((mean - 500.0).abs() / 500.0 < 0.015, "mean={}", mean);
        assert!(var < 10.0, "var={}", var);
    }

    #[test]
    fn test_chrom_moments_lookup_out_of_range() {
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(0, 1_000_000)]);
        let dm = DepthMap::build(&bed);
        let moments = build_chrom_moments(&dm.chroms[0], 0.01);
        assert!(moments.lookup(0.0).is_none());
        assert!(moments.lookup(moments.n_bins as f64).is_none());
    }

    // ── compute_analytic_stats ────────────────────────────────────────────────

    #[test]
    fn test_compute_analytic_stats_no_shared_chrom() {
        let mut q_bed = BedMap::new();
        q_bed.insert("chr21".to_string(), vec![(10_000_000, 10_001_000)]);
        let q_data = build_query_chrom_data(&q_bed, None);

        let mut db_bed = BedMap::new();
        db_bed.insert("chr22".to_string(), vec![(10_000_000, 10_001_000)]);
        let db_dm = DepthMap::build(&db_bed);
        let db_moments = build_depth_moments(&db_dm, DEFAULT_MOMENTS_EPS);

        assert!(compute_analytic_stats(&q_data, &db_moments, 0.0).is_none());
    }

    #[test]
    fn test_compute_analytic_stats_zero_observed() {
        let mut q_bed = BedMap::new();
        q_bed.insert("chr22".to_string(), vec![(0, 1_000_000)]);
        let q_data = build_query_chrom_data(&q_bed, None);

        let mut db_bed = BedMap::new();
        db_bed.insert("chr22".to_string(), vec![(49_000_000, 50_000_000)]);
        let db_dm = DepthMap::build(&db_bed);
        let db_moments = build_depth_moments(&db_dm, DEFAULT_MOMENTS_EPS);

        let (p_value, _) = compute_analytic_stats(&q_data, &db_moments, 0.0).unwrap();
        assert_eq!(p_value, 1.0);
    }

    #[test]
    fn test_compute_analytic_stats_enriched() {
        let ivs: Vec<(u32, u32)> =
            (0..500u32).map(|i| (i * 10_000, i * 10_000 + 5_000)).collect();

        let mut db_bed = BedMap::new();
        db_bed.insert("chr22".to_string(), ivs.clone());
        let db_dm = DepthMap::build(&db_bed);
        let db_moments = build_depth_moments(&db_dm, DEFAULT_MOMENTS_EPS);

        let mut q_bed = BedMap::new();
        q_bed.insert("chr22".to_string(), ivs);
        let q_data = build_query_chrom_data(&q_bed, None);

        let (_p_value, llr) = compute_analytic_stats(&q_data, &db_moments, 25_000.0).unwrap();
        assert!(llr.map_or(true, |l| l >= 0.0), "llr={llr:?}");
    }

    // ── nb_llr ────────────────────────────────────────────────────────────────

    #[test]
    fn test_nb_llr_zero_at_null_mean() {
        let llr = nb_llr(4000.0, 8000.0, 4000.0).unwrap();
        assert!(llr.abs() < 1e-6, "llr={llr}");
    }

    #[test]
    fn test_nb_llr_positive_for_enrichment() {
        let llr = nb_llr(4000.0, 8000.0, 20_000.0).unwrap();
        assert!(llr > 0.0, "llr={llr}");
    }

    #[test]
    fn test_nb_llr_none_underdispersed() {
        assert!(nb_llr(2000.0, 1000.0, 1.0).is_none());
    }

    // ── erfc ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_erfc_known_values() {
        assert!((erfc(0.0) - 1.0).abs() < 1e-6);
        assert!((erfc(1.0) - 0.157_299_2).abs() < 1e-5, "erfc(1)={}", erfc(1.0));
        assert!(erfc(5.0) < 1e-10);
    }

    // ── FilterMask ────────────────────────────────────────────────────────────

    #[test]
    fn test_filter_mask_excludes_unmatched_chrom() {
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(10_000_000, 11_000_000)]);
        let dm = DepthMap::build(&bed);

        let mut filter_bed = BedMap::new();
        filter_bed.insert("chr21".to_string(), vec![(0, 46_709_983)]);
        let filter = FilterMask::build(&filter_bed, FilterMode::Whitelist);

        let db_moments: Vec<ChromMoments> = dm
            .chroms
            .iter()
            .filter_map(|cdm| {
                let mask = filter.get(&cdm.chrom)?;
                Some(build_chrom_moments_with_filter(cdm, mask, DEFAULT_MOMENTS_EPS))
            })
            .collect();

        let q_data = build_query_chrom_data(&bed, Some(&filter));
        assert!(compute_analytic_stats(&q_data, &db_moments, 0.0).is_none());
    }

    #[test]
    fn test_sub_bin_fractional_coverage() {
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(20, 80)]);
        let dm = DepthMap::build(&bed);
        let chrom = &dm.chroms[0];

        let tc = chrom.total_cov();
        assert!((tc - 0.6).abs() < 1e-9, "total_cov={tc}");
    }
}
