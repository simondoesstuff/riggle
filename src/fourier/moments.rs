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
//! # Compact storage
//! Only a subset of L values are stored, tolerating a 1% relative error in L:
//!   - L = 1..T_THRESHOLD (T=100): every block size stored, no approximation.
//!   - L > T: only the first L in each log_{1+eps}-spaced slot is stored.
//!     The number of stored entries is T + ⌊log_{1.01}(N/T)⌋ ≈ 1100 for chr1.
//!
//! Use `compact_index(l)` to map a block size to its 0-indexed position in the store.

use std::cell::RefCell;
use std::collections::HashMap;

use rayon::prelude::*;
use realfft::RealFftPlanner;
use rustfft::num_complex::Complex;
use wide::f64x4;

// One real-FFT planner per rayon worker thread.  Plans for a given size are
// computed once and reused across all chromosomes and files on that thread.
thread_local! {
    static R2C_PLANNER: RefCell<RealFftPlanner<f64>> = RefCell::new(RealFftPlanner::new());
}

use super::depth::{ChromDepthMap, DepthMap, build_depth_signal};
use super::genome::{BedMap, BIN_SIZE, hg38_chrom_sizes};

/// Relative error tolerance for block-size approximation (1%).
pub const EPS: f64 = 0.01;
/// Block sizes L <= T_THRESHOLD are stored exactly (dense); above this,
/// only geometrically spaced samples are kept.
pub const T_THRESHOLD: usize = 100; // floor(1.0 / EPS)

/// Map block size L to its 0-indexed position in the compact moment store.
///
/// For L <= T_THRESHOLD: index = L - 1 (exact, no approximation).
/// For L > T_THRESHOLD: index = T + floor(log_{1+EPS}(L/T)) - 1.
///
/// Two queries L and L' that share an index differ by at most EPS in relative terms.
pub fn compact_index(l: usize) -> usize {
    debug_assert!(l >= 1, "compact_index: l must be >= 1");
    if l <= T_THRESHOLD {
        l - 1
    } else {
        let k = ((l as f64 / T_THRESHOLD as f64).ln() / (1.0_f64 + EPS).ln()) as usize;
        T_THRESHOLD + k - 1
    }
}

/// Compact sliding-window moments for one DB chromosome.
///
/// Stores one `(mean, var)` pair per log-slot: dense for L ≤ T_THRESHOLD,
/// then geometrically spaced above. Use `compact_index(l)` to look up.
pub struct ChromMoments {
    pub chrom: String,
    pub n_bins: usize,
    /// Compact (mean, var) pairs as f32; position given by `compact_index(L)`.
    pub moments: Vec<(f32, f32)>,
}

impl ChromMoments {
    /// Return `(mean, var)` for a sliding window of `l_bins` bins.
    ///
    /// Lookup is O(1). The returned moments may have been computed at the
    /// nearest sampled L (within EPS relative error). Returns `None` when
    /// `l_bins` rounds to 0 or ≥ `n_bins`.
    pub fn lookup(&self, l_bins: f64) -> Option<(f64, f64)> {
        let l = l_bins.ceil() as usize;
        if l == 0 || l >= self.n_bins {
            return None;
        }
        let (m, v) = self.moments.get(compact_index(l))?;
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
///
/// Iterates only the chromosomes present in `bed` (not the full hg38 list).
/// `interval_lengths` are sorted so the stats accumulator can group by unique
/// block size with a single linear pass instead of one lookup per interval.
pub fn build_query_chrom_data(bed: &BedMap) -> Vec<QueryChromData> {
    let size_map: HashMap<&str, u32> =
        hg38_chrom_sizes().iter().map(|&(c, s)| (c, s)).collect();
    bed.iter()
        .filter_map(|(chrom, ivs)| {
            if ivs.is_empty() {
                return None;
            }
            let &size = size_map.get(chrom.as_str())?;
            let n_bins = ((size + BIN_SIZE - 1) / BIN_SIZE) as usize;
            let mut interval_lengths: Vec<f64> = ivs
                .iter()
                .map(|&(s, e)| (e - s) as f64 / BIN_SIZE as f64)
                .collect();
            interval_lengths.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
            Some(QueryChromData {
                chrom: chrom.clone(),
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

    // d*[0..=N]: prefix sum of g  (d*[0]=0, d*[j] = Σ_{k<j} g[k])
    let d_star = prefix_sum_f32_to_f64(g);   // two-pass SIMD

    // d**[0..=N+1]: prefix sum of d*  (d**[0]=0)
    let d_ss = prefix_sum_f64(&d_star);       // two-pass SIMD

    // d2[0..=N+1]: prefix sum of (d*)²  with SIMD squaring + two-pass carry
    let d2 = d2_prefix_sum(&d_star);

    // Autocorrelation R[L] = IFFT(|FFT(d*)|²)[L]
    let r = autocorr_via_fft(&d_star);

    let max_l = n - 1;

    // Estimate compact store size: T linear + log-spaced geometric entries.
    let geom_est = if n > T_THRESHOLD {
        ((max_l as f64 / T_THRESHOLD as f64).ln() / (1.0_f64 + EPS).ln()).ceil() as usize + 1
    } else {
        0
    };
    let mut moments = Vec::with_capacity(T_THRESHOLD.min(max_l) + geom_est);

    // Phase 1: dense linear — every L = 1..=T (no approximation).
    for l in 1..=T_THRESHOLD.min(max_l) {
        let (mean, var) = window_moments(&d_ss, &d2, &r, n, l);
        moments.push((mean as f32, var as f32));
    }

    // Phase 2: geometric spacing for L > T.
    // Canonical L for each log-slot: first integer >= T * (1+EPS)^j.
    // Relative error in L is at most EPS for any query in that slot.
    let mut l_f = T_THRESHOLD as f64;
    let mut prev_l = T_THRESHOLD;
    loop {
        l_f *= 1.0 + EPS;
        let l = l_f.ceil() as usize;
        if l > max_l {
            break;
        }
        if l > prev_l {
            let (mean, var) = window_moments(&d_ss, &d2, &r, n, l);
            moments.push((mean as f32, var as f32));
            prev_l = l;
        }
    }

    ChromMoments {
        chrom: chrom.to_string(),
        n_bins,
        moments,
    }
}

/// Two-pass prefix sum: output[0]=0, output[j] = Σ_{k<j} src[k] (f32→f64).
///
/// Pass 1: independent local prefix sums within each chunk of 4 (no cross-chunk
///         carry, so chunks are processed without sequential dependency).
/// Pass 2: sequential prefix of chunk totals → per-chunk carries.
/// Pass 3: SIMD broadcast-add of carry to each chunk.
fn prefix_sum_f32_to_f64(src: &[f32]) -> Vec<f64> {
    let n = src.len();
    let mut d = vec![0.0f64; n + 1];

    const W: usize = 4;
    let chunks = n / W;

    for c in 0..chunks {
        let b = c * W;
        let g0 = src[b] as f64;
        let g1 = src[b + 1] as f64;
        let g2 = src[b + 2] as f64;
        let g3 = src[b + 3] as f64;
        d[b + 1] = g0;
        d[b + 2] = g0 + g1;
        d[b + 3] = g0 + g1 + g2;
        d[b + 4] = g0 + g1 + g2 + g3;
    }

    let mut carry = 0.0f64;
    let mut carries = vec![0.0f64; chunks];
    for c in 0..chunks {
        carries[c] = carry;
        carry += d[c * W + W];
    }

    for c in 1..chunks {
        let b = c * W;
        let cv = f64x4::splat(carries[c]);
        let out = (f64x4::from([d[b + 1], d[b + 2], d[b + 3], d[b + 4]]) + cv).to_array();
        d[b + 1] = out[0];
        d[b + 2] = out[1];
        d[b + 3] = out[2];
        d[b + 4] = out[3];
    }

    // Tail: d[chunks*W] is correct after pass 3 (or 0 if chunks==0).
    for i in chunks * W..n {
        d[i + 1] = d[i] + src[i] as f64;
    }
    d
}

/// Two-pass prefix sum: output[0]=0, output[j] = Σ_{k<j} src[k] (f64→f64).
fn prefix_sum_f64(src: &[f64]) -> Vec<f64> {
    let n = src.len();
    let mut d = vec![0.0f64; n + 1];

    const W: usize = 4;
    let chunks = n / W;

    for c in 0..chunks {
        let b = c * W;
        d[b + 1] = src[b];
        d[b + 2] = src[b] + src[b + 1];
        d[b + 3] = src[b] + src[b + 1] + src[b + 2];
        d[b + 4] = src[b] + src[b + 1] + src[b + 2] + src[b + 3];
    }

    let mut carry = 0.0f64;
    let mut carries = vec![0.0f64; chunks];
    for c in 0..chunks {
        carries[c] = carry;
        carry += d[c * W + W];
    }

    for c in 1..chunks {
        let b = c * W;
        let cv = f64x4::splat(carries[c]);
        let out = (f64x4::from([d[b + 1], d[b + 2], d[b + 3], d[b + 4]]) + cv).to_array();
        d[b + 1] = out[0];
        d[b + 2] = out[1];
        d[b + 3] = out[2];
        d[b + 4] = out[3];
    }

    for i in chunks * W..n {
        d[i + 1] = d[i] + src[i];
    }
    d
}

/// Two-pass prefix sum of squares of `src`, length `src.len() + 1`.
///
/// SIMD f64x4 squaring in pass 1; two-pass carry propagation eliminates the
/// sequential cross-chunk dependency of the naive approach.
fn d2_prefix_sum(src: &[f64]) -> Vec<f64> {
    let m = src.len();
    let mut d2 = vec![0.0f64; m + 1];

    const W: usize = 4;
    let chunks = m / W;

    for c in 0..chunks {
        let b = c * W;
        let v = f64x4::from([src[b], src[b + 1], src[b + 2], src[b + 3]]);
        let sq = (v * v).to_array();
        d2[b + 1] = sq[0];
        d2[b + 2] = sq[0] + sq[1];
        d2[b + 3] = sq[0] + sq[1] + sq[2];
        d2[b + 4] = sq[0] + sq[1] + sq[2] + sq[3];
    }

    let mut carry = 0.0f64;
    let mut carries = vec![0.0f64; chunks];
    for c in 0..chunks {
        carries[c] = carry;
        carry += d2[c * W + W];
    }

    for c in 1..chunks {
        let b = c * W;
        let cv = f64x4::splat(carries[c]);
        let out = (f64x4::from([d2[b + 1], d2[b + 2], d2[b + 3], d2[b + 4]]) + cv).to_array();
        d2[b + 1] = out[0];
        d2[b + 2] = out[1];
        d2[b + 3] = out[2];
        d2[b + 4] = out[3];
    }

    for i in chunks * W..m {
        d2[i + 1] = d2[i] + src[i] * src[i];
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

/// Linear autocorrelation of `signal` at all lags 0..len via zero-padded real FFT.
///
/// Uses a real-to-complex forward FFT and complex-to-real inverse FFT (roughly
/// half the work and memory of a complex-to-complex FFT on a real input).
/// Plans are cached per rayon worker thread via `R2C_PLANNER`.
///
/// Returns a vector of length `signal.len()` where `result[L]` =
/// Σ_{k=0}^{n-L-1} signal[k] * signal[k+L].
fn autocorr_via_fft(signal: &[f64]) -> Vec<f64> {
    let n = signal.len();
    if n == 0 {
        return Vec::new();
    }
    let fft_len = (2 * n).next_power_of_two();

    // Fetch (or compute) plans from the thread-local planner cache.
    let (r2c, c2r) = R2C_PLANNER.with(|p| {
        let mut p = p.borrow_mut();
        (p.plan_fft_forward(fft_len), p.plan_fft_inverse(fft_len))
    });

    // Real buffer: signal zero-padded to fft_len.
    let mut real_buf: Vec<f64> = Vec::with_capacity(fft_len);
    real_buf.extend_from_slice(signal);
    real_buf.resize(fft_len, 0.0);

    // Complex spectrum: only fft_len/2 + 1 elements needed (half of complex FFT).
    let mut spectrum = vec![Complex::new(0.0, 0.0); fft_len / 2 + 1];

    r2c.process(&mut real_buf, &mut spectrum).expect("r2c FFT failed");

    // Power spectrum in-place: |X[k]|² (real, so imaginary set to 0).
    for c in &mut spectrum {
        let p = c.norm_sqr();
        *c = Complex::new(p, 0.0);
    }

    // Inverse real FFT: writes back into real_buf (fft_len f64).
    c2r.process(&mut spectrum, &mut real_buf).expect("c2r IFFT failed");

    // realfft IFFT is unnormalized; divide by fft_len.
    let scale = 1.0 / fft_len as f64;
    real_buf.truncate(n);
    for x in &mut real_buf {
        *x *= scale;
    }
    real_buf
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
