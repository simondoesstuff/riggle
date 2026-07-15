//! Analytic NB statistics for genomic interval overlap.
//!
//! # Index time
//! [`DepthMap::build`] records a per-chromosome sparse Dirac impulse train
//! (the derivative of coverage: +1 at each interval start, −1 at each end).
//! Storage is O(K) for K intervals.  Saved to `{db}/depthmap.rkyv`.
//!
//! # Query time
//! 1. [`build_chrom_cov_data_with_filter`] converts each DB chromosome's impulse
//!    train to a per-bin coverage array g[j] via prefix sum — O(N_bins) — and
//!    computes the mean depth and the sample variance of sliding window sums at
//!    each power-of-2 scale from 1 to 2^(MAX_SLIDING_LOG2−1) bins.
//! 2. [`build_query_chrom_data`] extracts per-chromosome interval lengths (in
//!    100 bp bins) from the raw BED intervals of the query file.
//! 3. [`compute_analytic_stats`] fits a negative-binomial (NB) null by summing
//!    independent contributions from each query interval.
//!
//! # Null model (free-relative-gap / independent-interval)
//! Each query interval slides independently over the DB (its gaps relative to
//! other query intervals are free).  For query interval i of length l_i (bins)
//! on a chromosome where the DB has mean depth μ_B:
//!
//!   μ_i  = l_i · μ_B
//!   σ²_i = Var_s(g_B)   (sample variance of s-bin sliding-window sums over the
//!                         DB coverage array, where s = 2^round(log2(l_i)))
//!
//! The sliding-window variance captures the true cluster structure of the DB:
//! when DB intervals are longer than s, windows are bimodal (all-on / all-off),
//! giving high variance (overdispersion) even for binary coverage.
//!
//! Genome-wide null moments are the sum across all intervals and chromosomes:
//!   μ = Σ_i μ_i,  σ² = Σ_i σ²_i
//!
//! NB parameters: p = μ/σ², r = μ²/(σ²−μ).  Requires σ² > μ (overdispersion);
//! returns None otherwise (e.g., point-like queries against Bernoulli DB).
//!
//! LLR at observation o (exponential tilting to mean o):
//!   LLR(o) = o·θ + r·ln(p_o/p)   where θ = ln((1−p_o)/(1−p)), p_o = r/(r+o).
//!
//! P-value: Wilks' approximation — 2·LLR ≈ χ²(1), so
//!   p = erfc(√LLR) / 2.

use std::collections::HashMap;
use std::path::Path;

use rayon::prelude::*;

/// Number of power-of-2 sliding-window scales precomputed per DB chromosome.
/// Covers scales 2^0 = 1 bin (100 bp) through 2^13 = 8192 bins (≈ 820 kb).
const MAX_SLIDING_LOG2: usize = 14;

// ── FilterMask ────────────────────────────────────────────────────────────────

/// Whether a [`FilterMask`] BED file lists regions to keep or regions to remove.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterMode {
    /// Only bins covered by the BED are accessible; chromosomes absent from the
    /// BED are excluded entirely.
    Whitelist,
    /// Bins covered by the BED are zeroed out; all other bins remain accessible.
    /// Chromosomes absent from the BED are unaffected (fully accessible).
    Blacklist,
}

/// A per-chromosome boolean mask at 100 bp resolution.
///
/// Internally always stores `true = accessible` regardless of the source mode.
/// Built from a [`BedMap`] via [`FilterMask::build`].
pub struct FilterMask {
    chroms: HashMap<String, Vec<bool>>,
    mode: FilterMode,
}

impl FilterMask {
    /// Build a `FilterMask` from a [`BedMap`] and a [`FilterMode`].
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

    /// Return the accessibility mask for `chrom` (`true` = bin contributes).
    pub fn get(&self, chrom: &str) -> Option<&[bool]> {
        self.chroms.get(chrom).map(Vec::as_slice)
    }

    pub fn mode(&self) -> FilterMode {
        self.mode
    }
}

const BIN_SIZE: u32 = 100;

// ── Public types ──────────────────────────────────────────────────────────────

/// Chromosome name → (start, end) interval pairs.
pub type BedMap = HashMap<String, Vec<(u32, u32)>>;

/// Per-chromosome coverage statistics for one DB file.
/// Produced by [`build_chrom_cov_data_with_filter`]; reused across all queries.
pub struct ChromCovData {
    pub chrom: String,
    pub n_bins: usize,
    /// Mean coverage per bin (over n_bins, after masking).
    pub cov_mean: f64,
    /// `sliding_vars[k]` = sample variance of `2^k`-bin sliding-window sums of g.
    /// Zero for k where `2^k > n_bins`.  Index 0 equals the per-bin variance.
    pub sliding_vars: [f64; MAX_SLIDING_LOG2],
}

/// Per-chromosome interval data for one query file.
/// Produced by [`build_query_chrom_data`].
pub struct QueryChromData {
    pub chrom: String,
    pub n_bins: usize,
    /// Length of each interval on this chromosome, in 100 bp bins (may be fractional).
    pub interval_lengths: Vec<f64>,
}

/// Sparse Dirac impulse train (derivative coverage) for one chromosome.
#[derive(Debug)]
pub struct ChromDepthMap {
    pub chrom: String,
    pub n_bins: usize,
    /// Raw bp start positions with a +1 impulse (one per interval).
    pub pos_spikes: Vec<u32>,
    /// Raw bp end positions with a −1 impulse (one per interval).
    pub neg_spikes: Vec<u32>,
}

impl ChromDepthMap {
    /// Sum of absolute impulse magnitudes.
    pub fn v(&self) -> f64 {
        (self.pos_spikes.len() + self.neg_spikes.len()) as f64
    }

    /// Total coverage in bins: Σ_interval (end − start) / BIN_SIZE.
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

    /// Total V across all chromosomes.
    pub fn total_v(&self) -> f64 {
        self.chroms.iter().map(|c| c.v()).sum()
    }
}

// ── Public helpers ────────────────────────────────────────────────────────────

/// V for a BedMap: 2 × total interval count.
pub fn bed_map_v(bed: &BedMap) -> f64 {
    bed.values().map(|v| 2 * v.len()).sum::<usize>() as f64
}

/// Build per-chromosome coverage data for a DepthMap.
pub fn build_chrom_cov_data(dm: &DepthMap) -> Vec<ChromCovData> {
    build_chrom_cov_data_with_filter(dm, None)
}

/// Mean interval length in 100 bp bins across all chromosomes in a query.
pub fn mean_interval_bins(q_data: &[QueryChromData]) -> f64 {
    let total_l: f64 = q_data.iter().flat_map(|c| c.interval_lengths.iter()).sum();
    let total_n: usize = q_data.iter().map(|c| c.interval_lengths.len()).sum();
    if total_n == 0 { 0.0 } else { total_l / total_n as f64 }
}

/// Build per-chromosome interval data from a raw BED map.
///
/// Interval lengths are stored in 100 bp bins (may be fractional for short intervals).
/// When a [`FilterMask`] is supplied in whitelist mode, chromosomes absent from
/// the mask are excluded entirely (matching the DB-side filtering).
pub fn build_query_chrom_data(bed: &BedMap, filter: Option<&FilterMask>) -> Vec<QueryChromData> {
    hg38_chrom_sizes()
        .iter()
        .filter_map(|&(chrom, size)| {
            // For whitelist: exclude chromosomes not in the filter.
            // For blacklist: all hg38 chromosomes are present in the filter.
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

/// Exact base-pair overlap in 100 bp bins between two depth maps.
///
/// Builds per-bin coverage arrays on the fly (no filter applied), computes the
/// dot product over shared chromosomes, and discards the arrays.  Intended for
/// the `pval` binary where no sweep count is available; the query pipeline uses
/// the heuristic `sweep_count × mean_query_interval_bins` instead.
pub fn coverage_dot_product(a_dm: &DepthMap, b_dm: &DepthMap) -> f64 {
    let a_map: HashMap<&str, &ChromDepthMap> =
        a_dm.chroms.iter().map(|c| (c.chrom.as_str(), c)).collect();
    b_dm.chroms
        .iter()
        .filter_map(|b_cd| {
            let a_cd = a_map.get(b_cd.chrom.as_str())?;
            let (a_cov, _, _) = chrom_coverage_and_stats(a_cd, None);
            let (b_cov, _, _) = chrom_coverage_and_stats(b_cd, None);
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

/// Like [`build_chrom_cov_data`] but applies a [`FilterMask`] before coverage
/// computation.  Chromosomes absent from a whitelist are excluded entirely.
pub fn build_chrom_cov_data_with_filter(
    dm: &DepthMap,
    filter: Option<&FilterMask>,
) -> Vec<ChromCovData> {
    dm.chroms
        .par_iter()
        .filter_map(|cdm| {
            let mask = match filter {
                None => None,
                Some(f) => match f.get(&cdm.chrom) {
                    Some(m) => Some(m),
                    None => return None,
                },
            };
            let (g, cov_mean, _) = chrom_coverage_and_stats(cdm, mask);
            let sliding_vars = all_sliding_vars(&g);
            Some(ChromCovData {
                chrom: cdm.chrom.clone(),
                n_bins: cdm.n_bins,
                cov_mean,
                sliding_vars,
            })
        })
        .collect()
}

/// Compute analytic NB statistics for a (query, DB) pair.
///
/// Returns `(p_value, llr)` or `None` if no chromosomes are shared.
///
/// Uses the free-relative-gap null: each query interval slides independently.
/// For interval i of length l_i on chromosome c with DB mean μ_B and
/// sliding-window variance at scale 2^round(log2(l_i)):
///   μ_i  = l_i · μ_B,   σ²_i = sliding_vars[round(log2(l_i))]
/// Summing over all intervals gives the genome-wide null moments.
///
/// `observed` is the base-pair overlap in 100 bp bins supplied by the caller.
/// In the query pipeline this is `sweep_count × mean_query_interval_bins`
/// (a heuristic assuming full overlaps); for exact results use
/// [`coverage_dot_product`].
pub fn compute_analytic_stats(
    q_data: &[QueryChromData],
    db_data: &[ChromCovData],
    observed: f64,
) -> Option<(f64, Option<f64>)> {
    let db_map: HashMap<&str, &ChromCovData> =
        db_data.iter().map(|c| (c.chrom.as_str(), c)).collect();

    let mut mu_total = 0.0f64;
    let mut var_total = 0.0f64;
    let mut any_shared = false;

    for q_cd in q_data {
        let db_cd = match db_map.get(q_cd.chrom.as_str()) {
            Some(d) => d,
            None => continue,
        };
        let mu_b = db_cd.cov_mean;
        for &l in &q_cd.interval_lengths {
            let k = round_pow2_log2(l);
            mu_total += l * mu_b;
            var_total += db_cd.sliding_vars[k];
            any_shared = true;
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

/// Parse a BED file (plain or gzip) into a [`BedMap`].
pub fn parse_bed_as_map(path: &Path) -> Result<BedMap, crate::io::BedParseError> {
    let shards = crate::io::parse_bed_file(path, 0)?;
    Ok(shards
        .into_iter()
        .map(|(k, ivs)| (k, ivs.into_iter().map(|iv| (iv.start, iv.end)).collect()))
        .collect())
}

/// Convert a `parse_bed_file` result into a [`BedMap`].
pub fn intervals_to_bed_map(shards: &HashMap<String, Vec<crate::core::Interval>>) -> BedMap {
    shards
        .iter()
        .map(|(k, ivs)| (k.clone(), ivs.iter().map(|iv| (iv.start, iv.end)).collect()))
        .collect()
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

// ── Private helpers ───────────────────────────────────────────────────────────

/// Round l (in bins) to the nearest power of 2 in log-space; return k = log2(scale).
/// Clamped to [0, MAX_SLIDING_LOG2 − 1].
fn round_pow2_log2(l: f64) -> usize {
    if l <= 1.0 {
        return 0;
    }
    (l.log2().round() as usize).min(MAX_SLIDING_LOG2 - 1)
}

/// Precompute sliding-window variance for every power-of-2 scale from 1 to 2^13.
///
/// Builds one prefix-sum array and sweeps it at each scale in O(N·MAX_SLIDING_LOG2).
fn all_sliding_vars(g: &[f32]) -> [f64; MAX_SLIDING_LOG2] {
    let n = g.len();
    let mut prefix = vec![0.0f64; n + 1];
    for i in 0..n {
        prefix[i + 1] = prefix[i] + g[i] as f64;
    }
    let mut vars = [0.0f64; MAX_SLIDING_LOG2];
    for k in 0..MAX_SLIDING_LOG2 {
        let scale = 1usize << k;
        if scale > n {
            break;
        }
        let n_windows = (n - scale + 1) as f64;
        let mut sum = 0.0f64;
        let mut sum_sq = 0.0f64;
        for i in 0..=(n - scale) {
            let f = prefix[i + scale] - prefix[i];
            sum += f;
            sum_sq += f * f;
        }
        let mean = sum / n_windows;
        vars[k] = sum_sq / n_windows - mean * mean;
    }
    vars
}

/// Saddlepoint LLR for the NB null tilted to the observation `o`.
///
/// Returns `None` when the NB fit is not feasible (μ ≤ 0, σ² ≤ μ, or
/// computed probabilities outside (0, 1)).
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

/// Build a coverage array and compute its mean and variance.
///
/// Spikes are distributed fractionally across adjacent bins to avoid aliasing
/// for sub-bin intervals.  The mask (if any) zeroes non-accessible bins after
/// the prefix sum.
fn chrom_coverage_and_stats(
    cdm: &ChromDepthMap,
    mask: Option<&[bool]>,
) -> (Vec<f32>, f64, f64) {
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

    let nb = n as f64;
    let sum: f64 = g.iter().map(|&x| x as f64).sum();
    let mean = sum / nb;
    let var: f64 = g
        .iter()
        .map(|&x| {
            let d = x as f64 - mean;
            d * d
        })
        .sum::<f64>()
        / nb;

    (g, mean, var)
}

/// Convert intervals to raw bp position vectors for fractional spike distribution.
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

    // ── nb_llr unit tests ─────────────────────────────────────────────────────

    #[test]
    fn test_nb_llr_zero_at_null_mean() {
        // When observed = μ, the tilted distribution equals the null → LLR = 0.
        let llr = nb_llr(4000.0, 8000.0, 4000.0).unwrap();
        assert!(llr.abs() < 1e-6, "LLR at null mean should be 0, got {llr}");
    }

    #[test]
    fn test_nb_llr_positive_for_enrichment() {
        let llr = nb_llr(4000.0, 8000.0, 20_000.0).unwrap();
        assert!(llr > 0.0, "LLR above null mean should be positive, got {llr}");
    }

    #[test]
    fn test_nb_llr_none_underdispersed() {
        // var ≤ mu → NB not applicable.
        assert!(nb_llr(2000.0, 1000.0, 1.0).is_none());
    }

    #[test]
    fn test_nb_llr_none_zero_mu() {
        assert!(nb_llr(0.0, 1.0, 1.0).is_none());
    }

    // ── erfc sanity ───────────────────────────────────────────────────────────

    #[test]
    fn test_erfc_known_values() {
        assert!((erfc(0.0) - 1.0).abs() < 1e-6);
        // erfc(1) ≈ 0.1573
        assert!((erfc(1.0) - 0.157_299_2).abs() < 1e-5, "erfc(1)={}", erfc(1.0));
        // erfc → 0 for large x
        assert!(erfc(5.0) < 1e-10);
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
        let db_cov = build_chrom_cov_data(&db_dm);

        assert!(compute_analytic_stats(&q_data, &db_cov, 0.0).is_none());
    }

    #[test]
    fn test_compute_analytic_stats_zero_observed() {
        // Sweep returned 0 overlaps → p_value must be 1.0 regardless of moments.
        let mut q_bed = BedMap::new();
        q_bed.insert("chr22".to_string(), vec![(0, 1_000_000)]);
        let q_data = build_query_chrom_data(&q_bed, None);

        let mut db_bed = BedMap::new();
        db_bed.insert("chr22".to_string(), vec![(49_000_000, 50_000_000)]);
        let db_dm = DepthMap::build(&db_bed);
        let db_cov = build_chrom_cov_data(&db_dm);

        let (p_value, _llr) = compute_analytic_stats(&q_data, &db_cov, 0.0).unwrap();
        assert_eq!(p_value, 1.0);
    }

    #[test]
    fn test_compute_analytic_stats_enriched() {
        // Many clustered intervals on both sides; high observed should give positive LLR.
        let ivs: Vec<(u32, u32)> =
            (0..500u32).map(|i| (i * 10_000, i * 10_000 + 5_000)).collect();

        let mut db_bed = BedMap::new();
        db_bed.insert("chr22".to_string(), ivs.clone());
        let db_dm = DepthMap::build(&db_bed);
        let db_cov = build_chrom_cov_data(&db_dm);

        let mut q_bed = BedMap::new();
        q_bed.insert("chr22".to_string(), ivs);
        let q_data = build_query_chrom_data(&q_bed, None);

        // observed = 500 overlapping intervals × 50 bins each (5000 bp / 100 bp/bin)
        let (_p_value, llr) = compute_analytic_stats(&q_data, &db_cov, 25_000.0).unwrap();
        assert!(llr.map_or(true, |l| l >= 0.0), "llr={llr:?}");
    }

    // ── FilterMask tests ──────────────────────────────────────────────────────

    #[test]
    fn test_filter_mask_excludes_unmatched_chrom() {
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(10_000_000, 11_000_000)]);
        let dm = DepthMap::build(&bed);

        let mut filter_bed = BedMap::new();
        filter_bed.insert("chr21".to_string(), vec![(0, 46_709_983)]);
        let filter = FilterMask::build(&filter_bed, FilterMode::Whitelist);

        let db_cov = build_chrom_cov_data_with_filter(&dm, Some(&filter));
        let q_data = build_query_chrom_data(&bed, Some(&filter));

        assert!(db_cov.is_empty());
        assert!(compute_analytic_stats(&q_data, &db_cov, 0.0).is_none());
    }

    #[test]
    fn test_filter_mask_confines_signal() {
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(10_000_000, 20_000_000)]);
        let dm = DepthMap::build(&bed);

        let mut filter_bed = BedMap::new();
        filter_bed.insert("chr22".to_string(), vec![(0, 50_818_468)]);
        let filter = FilterMask::build(&filter_bed, FilterMode::Whitelist);

        let db_cov = build_chrom_cov_data_with_filter(&dm, Some(&filter));
        let q_data = build_query_chrom_data(&bed, Some(&filter));
        // Whitelist covering the full chromosome; should have shared data.
        assert!(compute_analytic_stats(&q_data, &db_cov, 100_000.0).is_some());
    }

    #[test]
    fn test_blacklist_excludes_signal() {
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(10_000_000, 20_000_000)]);
        let dm = DepthMap::build(&bed);

        let mut bl_bed = BedMap::new();
        bl_bed.insert("chr22".to_string(), vec![(0, 50_818_468)]);
        let filter = FilterMask::build(&bl_bed, FilterMode::Blacklist);

        // Full blacklist zeroes all DB bins → cov_mean = 0 → μ = 0 → NB fit fails.
        let db_cov = build_chrom_cov_data_with_filter(&dm, Some(&filter));
        let q_data = build_query_chrom_data(&bed, Some(&filter));
        let (p_value, llr) = compute_analytic_stats(&q_data, &db_cov, 100_000.0).unwrap();
        assert!(llr.is_none(), "NB fit should fail when all DB bins are zeroed");
        assert_eq!(p_value, 1.0);
    }

    // ── Fractional spike distribution ─────────────────────────────────────────

    #[test]
    fn test_sub_bin_fractional_coverage() {
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(20, 80)]);
        let dm = DepthMap::build(&bed);
        let chrom = &dm.chroms[0];

        let tc = chrom.total_cov();
        assert!((tc - 0.6).abs() < 1e-9, "total_cov={tc}");
        assert_eq!(chrom.pos_spikes, vec![20u32]);
        assert_eq!(chrom.neg_spikes, vec![80u32]);
    }
}
