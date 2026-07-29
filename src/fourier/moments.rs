//! FFT-based sliding-window moment computation.
//!
//! # Per-chromosome pipeline (O(N log N))
//! For depth signal g[0..N-1] (built from the sparse impulse train):
//!
//!   d*[j]  = Σ_{k<j} g[k]               (prefix sum of g, size N+1)
//!   d**[j] = Σ_{k<j} d*[k]              (double prefix sum, size N+2)
//!   d2[j]  = Σ_{k<j} (d*[k])²           (SIMD-accelerated squared prefix, size N+2)
//!   R[L]   = IFFT(|FFT(d*)|²)[L]        (autocorrelation via zero-padded FFT)
//!
//! For each block size L = 1..N-1:
//!   n_w      = N - L + 1
//!   sum_f    = d**[N+1] − d**[L] − d**[N−L+1]
//!   sum_f2   = (d2[N+1] − d2[L]) + d2[N−L+1] − 2·R[L]
//!   mean(L)  = sum_f / n_w
//!   var(L)   = sum_f2 / n_w − mean(L)²
//!
//! # Dense storage
//! All N-1 (mean, var) pairs are stored as f32 in `moments[L-1]`.
//! Query lookup for block size L is a direct O(1) array access with no approximation.

use rayon::prelude::*;
use rustfft::num_complex::Complex;
use rustfft::FftPlanner;
use wide::f64x4;

use super::depth::{ChromDepthMap, DepthMap, build_depth_signal};
use super::genome::{BedMap, BIN_SIZE, hg38_chrom_sizes};

/// Dense sliding-window moments for one DB chromosome.
///
/// `moments[L-1]` = `(mean, var)` for a uniformly placed block of L bins.
pub struct ChromMoments {
    pub chrom: String,
    pub n_bins: usize,
    /// (mean, var) pairs for L = 1..n_bins-1, stored as f32.
    pub moments: Vec<(f32, f32)>,
}

impl ChromMoments {
    /// Return `(mean, var)` for a sliding window of `l_bins` bins.
    ///
    /// Returns `None` when `l_bins` rounds to 0 or ≥ `n_bins`.
    pub fn lookup(&self, l_bins: f64) -> Option<(f64, f64)> {
        let l = l_bins.round() as usize;
        if l == 0 || l >= self.n_bins {
            return None;
        }
        let (m, v) = self.moments.get(l - 1)?;
        Some((*m as f64, *v as f64))
    }
}

/// Per-chromosome interval data for one query file.
pub struct QueryChromData {
    pub chrom: String,
    pub n_bins: usize,
    /// Length of each interval in 100 bp bins (may be fractional).
    pub interval_lengths: Vec<f64>,
}

/// Build dense moment tables for every chromosome in a [`DepthMap`].
pub fn build_depth_moments(dm: &DepthMap) -> Vec<ChromMoments> {
    dm.chroms
        .par_iter()
        .map(build_chrom_moments)
        .collect()
}

/// Build dense moments for one chromosome.
pub fn build_chrom_moments(cdm: &ChromDepthMap) -> ChromMoments {
    let g = build_depth_signal(cdm);
    moments_from_signal(&g, &cdm.chrom, cdm.n_bins)
}

/// Build per-chromosome interval data from a raw BED map.
pub fn build_query_chrom_data(bed: &BedMap) -> Vec<QueryChromData> {
    hg38_chrom_sizes()
        .iter()
        .filter_map(|&(chrom, size)| {
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

fn moments_from_signal(g: &[f32], chrom: &str, n_bins: usize) -> ChromMoments {
    let n = n_bins;

    if n <= 1 {
        return ChromMoments {
            chrom: chrom.to_string(),
            n_bins,
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

    // d2[0..=N+1]: prefix sum of (d*)² using SIMD squaring for throughput
    let d2 = d2_prefix_sum(&d_star);

    // Autocorrelation R[L] = IFFT(|FFT(d*)|²)[L]
    let r = autocorr_via_fft(&d_star);

    // Build all (mean, var) pairs for L = 1..n_bins-1
    let max_l = n - 1;
    let mut moments = Vec::with_capacity(max_l);
    for l in 1..=max_l {
        let (mean, var) = window_moments(&d_ss, &d2, &r, n, l);
        moments.push((mean as f32, var as f32));
    }

    ChromMoments {
        chrom: chrom.to_string(),
        n_bins,
        moments,
    }
}

/// Prefix sum of squares of `d_star`, returned as a vec of length `d_star.len() + 1`.
///
/// Uses f64x4 SIMD to compute four squares per iteration before the sequential
/// prefix-sum step, improving throughput on the squaring bottleneck.
fn d2_prefix_sum(d_star: &[f64]) -> Vec<f64> {
    let m = d_star.len(); // = n+1
    let mut d2 = vec![0.0f64; m + 1]; // d2[0]=0, d2[j] = Σ_{k<j} d_star[k]²

    const W: usize = 4;
    let chunks = m / W;

    for c in 0..chunks {
        let b = c * W;
        let v = f64x4::from([d_star[b], d_star[b + 1], d_star[b + 2], d_star[b + 3]]);
        let sq = (v * v).to_array();
        // Sequential prefix additions using the SIMD-computed squares.
        for i in 0..W {
            d2[b + i + 1] = d2[b + i] + sq[i];
        }
    }
    for i in chunks * W..m {
        d2[i + 1] = d2[i] + d_star[i] * d_star[i];
    }

    d2
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

    // ── d2_prefix_sum ─────────────────────────────────────────────────────────

    #[test]
    fn test_d2_prefix_sum_correctness() {
        let d_star: Vec<f64> = (0..20).map(|i| i as f64 * 0.7).collect();

        let got = d2_prefix_sum(&d_star);

        // Brute-force reference
        let mut expected = vec![0.0f64; d_star.len() + 1];
        for j in 1..=d_star.len() {
            expected[j] = expected[j - 1] + d_star[j - 1] * d_star[j - 1];
        }

        for (i, (e, g)) in expected.iter().zip(got.iter()).enumerate() {
            assert!((e - g).abs() < 1e-9, "index {i}: expected {e}, got {g}");
        }
    }

    // ── window_moments correctness ────────────────────────────────────────────

    #[test]
    fn test_window_moments_uniform_signal() {
        // g = [1,1,1,...,1] (N bins): every window of any L sums to L.
        // mean(L) = L, var(L) = 0.
        let n = 50usize;
        let g: Vec<f32> = vec![1.0; n];

        let mut d_star = vec![0.0f64; n + 1];
        for j in 1..=n {
            d_star[j] = d_star[j - 1] + g[j - 1] as f64;
        }
        let mut d_ss = vec![0.0f64; n + 2];
        for j in 1..=(n + 1) {
            d_ss[j] = d_ss[j - 1] + d_star[j - 1];
        }
        let d2 = d2_prefix_sum(&d_star);
        let r = autocorr_via_fft(&d_star);

        for l in 1..n {
            let (mean, var) = window_moments(&d_ss, &d2, &r, n, l);
            assert!(
                (mean - l as f64).abs() < 1e-6,
                "uniform: mean(L={})={:.4} expected {}",
                l, mean, l
            );
            assert!(var < 1e-6, "uniform: var(L={})={:.4} expected 0", l, var);
        }
    }

    #[test]
    fn test_window_moments_sparse_periodic() {
        // g = [1,0,0,1,0,0,1] (N=7), L=3: every window sums to 1 → var=0, mean=1.
        let g: Vec<f32> = vec![1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0];
        let n = g.len();
        let mut d_star = vec![0.0f64; n + 1];
        for j in 1..=n {
            d_star[j] = d_star[j - 1] + g[j - 1] as f64;
        }
        let mut d_ss = vec![0.0f64; n + 2];
        for j in 1..=(n + 1) {
            d_ss[j] = d_ss[j - 1] + d_star[j - 1];
        }
        let d2 = d2_prefix_sum(&d_star);
        let r = autocorr_via_fft(&d_star);

        let (mean, var) = window_moments(&d_ss, &d2, &r, n, 3);
        assert!((mean - 1.0).abs() < 1e-6, "mean={}", mean);
        assert!(var < 1e-6, "var={}", var);

        let (mean1, _) = window_moments(&d_ss, &d2, &r, n, 1);
        assert!((mean1 - 3.0 / 7.0).abs() < 1e-6, "mean(L=1)={}", mean1);
    }

    // ── ChromMoments lookup ───────────────────────────────────────────────────

    #[test]
    fn test_chrom_moments_lookup() {
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(10_000_000, 20_000_000)]);
        let dm = DepthMap::build(&bed);
        let moments = build_chrom_moments(&dm.chroms[0]);

        let r = moments.lookup(50.0);
        assert!(r.is_some(), "lookup(50) returned None");
        let (mean, var) = r.unwrap();
        assert!(mean > 0.0, "mean={}", mean);
        assert!(var >= 0.0, "var={}", var);
    }

    #[test]
    fn test_chrom_moments_lookup_large_l() {
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(0, 50_818_468)]);
        let dm = DepthMap::build(&bed);
        let moments = build_chrom_moments(&dm.chroms[0]);

        let r = moments.lookup(500.0);
        assert!(r.is_some(), "lookup(500) returned None");
        let (mean, _) = r.unwrap();
        // For a fully covered chromosome, mean(L) ≈ L
        assert!((mean - 500.0).abs() / 500.0 < 0.02, "mean={}", mean);
    }

    #[test]
    fn test_chrom_moments_lookup_out_of_range() {
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(0, 1_000_000)]);
        let dm = DepthMap::build(&bed);
        let moments = build_chrom_moments(&dm.chroms[0]);
        assert!(moments.lookup(0.0).is_none());
        assert!(moments.lookup(moments.n_bins as f64).is_none());
    }
}
