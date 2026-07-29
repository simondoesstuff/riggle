//! FFT-based sliding-window moment computation.
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

use rayon::prelude::*;
use rustfft::num_complex::Complex;
use rustfft::FftPlanner;

use super::depth::{ChromDepthMap, DepthMap, build_depth_signal};
use super::filter::FilterMask;
use super::genome::{BedMap, BIN_SIZE, hg38_chrom_sizes};

pub const DEFAULT_MOMENTS_EPS: f64 = 0.01;

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

pub fn mean_interval_bins(q_data: &[QueryChromData]) -> f64 {
    let total_l: f64 = q_data.iter().flat_map(|c| c.interval_lengths.iter()).sum();
    let total_n: usize = q_data.iter().map(|c| c.interval_lengths.len()).sum();
    if total_n == 0 { 0.0 } else { total_l / total_n as f64 }
}

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

    let r = autocorr_via_fft(&d_star);

    let max_l = n - 1;
    let mut moments: Vec<(f64, f64)> = Vec::new();

    // Dense region: L = 1..=dense_max (one entry per integer L)
    for l in 1..=dense_max.min(max_l) {
        moments.push(window_moments(&d_ss, &d2, &r, n, l));
    }

    // Sparse region: representative Ls at dense_max·(1+ε)^m for m = 1, 2, …
    let log_base = 1.0 + eps;
    let mut l_f = dense_max as f64 * log_base;
    let mut last_l = dense_max;

    loop {
        let l = l_f.floor() as usize;
        if l > max_l {
            break;
        }
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

        let r = moments.lookup(500.0);
        assert!(r.is_some(), "lookup(500) returned None");
        let (mean, var) = r.unwrap();
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
}
