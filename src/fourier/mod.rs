//! Fourier-based p-value for genomic interval overlap.
//!
//! # Index time
//! [`DepthMap::build`] records a per-chromosome sparse Dirac impulse train
//! (the derivative of coverage: +1 at each interval start, −1 at each end).
//! Storage is O(K) for K intervals.  Saved to `{db}/depthmap/{sid}.bin`.
//!
//! # Query time (two-phase)
//! Phase A — build DB coverage spectra once per DB file (cached):
//!   [`build_db_spectra`] converts each chromosome's impulse train to a
//!   truncated coverage spectrum G_B[k] via `realfft` (real→complex FFT).
//!   Coverage is built O(N) via prefix sum of the impulse train, then FFT'd
//!   in O(N log N).  Only the first M coefficients are stored.
//!
//! Phase B — per (query, DB) pair:
//!   1. Build query spectrum G_Q[k] (same real FFT pipeline).
//!   2. Accumulate cross-correlation: C_total[k] += G_Q[k] · conj(G_B[k])
//!      over all shared chromosomes.
//!   3. Single M-length complex IFFT → M samples c[s] of total overlap vs
//!      a circular shift of s·(N/M) bins.
//!   4. Empirical right-tailed p-value: P(overlap ≥ observed) =
//!      fraction of the M shift samples with c[s] ≥ c[0].
//!
//! # Frequency truncation (variance-adaptive)
//! All N/2+1 cross-correlation coefficients are accumulated first.  Then the
//! cumulative power Σ_{k=0}^{M-1} |C[k]|² is scanned in ascending-frequency
//! order until it reaches `variance_threshold × total_power`.  Only those M
//! coefficients are passed to the IFFT, giving M shift samples at coarser
//! resolution.  Lower thresholds (e.g. 0.90) act as a smoothing regulariser;
//! higher thresholds (e.g. 0.9999) approach the full-spectrum result.
//!
//! # Null model
//! All chromosomes are shifted by the same sample index s simultaneously.
//! This is equivalent to adding cross-correlation spectra before the IFFT
//! (linearity of IFFT) and gives a single M-sample null distribution for
//! the combined multi-chromosome overlap.

use std::cell::RefCell;
use std::collections::HashMap;
use std::path::Path;

use rayon::prelude::*;
use realfft::RealFftPlanner;
use rustfft::{FftPlanner, num_complex::Complex};

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
    ///
    /// **Whitelist**: only bins covered by `bed` are marked accessible; chromosomes
    /// absent from `bed` are excluded entirely.
    ///
    /// **Blacklist**: all bins start accessible; bins covered by `bed` are marked
    /// inaccessible.  Chromosomes absent from `bed` remain fully accessible.
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

    /// Return the accessibility mask for `chrom` (`true` = bin contributes to FFT).
    ///
    /// For a whitelist, returns `None` if the chromosome was absent from the BED
    /// (fully inaccessible).  For a blacklist, always returns `Some` since every
    /// chromosome is represented.
    pub fn get(&self, chrom: &str) -> Option<&[bool]> {
        self.chroms.get(chrom).map(Vec::as_slice)
    }

    pub fn mode(&self) -> FilterMode {
        self.mode
    }
}

thread_local! {
    static REAL_PLANNER: RefCell<RealFftPlanner<f32>> = RefCell::new(RealFftPlanner::new());
    static COMPLEX_PLANNER: RefCell<FftPlanner<f32>> = RefCell::new(FftPlanner::new());
}

const BIN_SIZE: u32 = 100;

/// Default fraction of cross-correlation power retained before the IFFT.
/// Near 1.0, the truncation is nearly invisible; lower values (e.g. 0.90)
/// act as a smoothing regulariser and can sharpen biological signal.
pub const DEFAULT_VARIANCE_THRESHOLD: f64 = 1.0;

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
    /// G_B[k] for k = 0..m (coverage spectrum, truncated to m frequencies).
    pub spec: Vec<Complex<f32>>,
}

/// Sparse Dirac impulse train (derivative coverage) for one chromosome.
///
/// Spikes are stored as raw bp positions.  At query time each spike is
/// fractionally distributed across two adjacent bins (linear interpolation
/// between `floor(x/BIN)` and `floor(x/BIN)+1`) to eliminate aliasing for
/// intervals smaller than `BIN_SIZE`.
#[derive(Debug)]
pub struct ChromDepthMap {
    pub chrom: String,
    pub n_bins: usize,
    /// Raw bp start positions with a +1 impulse (one per interval).
    pub pos_spikes: Vec<u32>,
    /// Raw bp end positions with a -1 impulse (one per interval).
    pub neg_spikes: Vec<u32>,
}

impl ChromDepthMap {
    /// Sum of absolute impulse magnitudes.
    pub fn v(&self) -> f64 {
        (self.pos_spikes.len() + self.neg_spikes.len()) as f64
    }
    /// Total coverage in bins: Σ_interval (end − start) / BIN_SIZE.
    ///
    /// With fractional spike distribution, this equals
    ///   Σ_pos (n_bins − x/BIN_SIZE) − Σ_neg (n_bins − x/BIN_SIZE).
    pub fn total_cov(&self) -> f64 {
        let n = self.n_bins as f64;
        let b = BIN_SIZE as f64;
        let pos: f64 = self.pos_spikes.iter().map(|&x| n - x as f64 / b).sum();
        let neg: f64 = self.neg_spikes.iter().map(|&x| n - x as f64 / b).sum();
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

/// Build per-chromosome coverage spectra for a DepthMap.
/// Stores all N/2+1 coefficients from the real forward FFT.
/// Runs in parallel over chromosomes.
pub fn build_db_spectra(dm: &DepthMap) -> Vec<ChromDbSpec> {
    build_db_spectra_with_filter(dm, None)
}

/// Like [`build_db_spectra`] but applies a [`FilterMask`] before the FFT.
///
/// For each chromosome:
/// - If `filter` is `None`, the full coverage vector is used (same as
///   [`build_db_spectra`]).
/// - If `filter` is `Some(f)` and `f` has an entry for the chromosome, the
///   coverage vector is multiplied element-wise by the boolean mask (non-
///   whitelisted bins are zeroed out) before the forward FFT.
/// - If `filter` is `Some(f)` but `f` has **no** entry for the chromosome, that
///   chromosome is excluded entirely (no `ChromDbSpec` is emitted for it).
///
/// This confines both the observed overlap and the shift-null distribution to
/// the same accessible genomic space, implementing a positional prior.
pub fn build_db_spectra_with_filter(dm: &DepthMap, filter: Option<&FilterMask>) -> Vec<ChromDbSpec> {
    dm.chroms
        .par_iter()
        .filter_map(|cdm| {
            let mask = match filter {
                None => None,
                Some(f) => match f.get(&cdm.chrom) {
                    Some(m) => Some(m),
                    None => return None, // chromosome excluded by filter
                },
            };
            let spec = chrom_coverage_spectrum(cdm, mask);
            Some(ChromDbSpec {
                chrom: cdm.chrom.clone(),
                n: cdm.n(),
                n_bins: cdm.n_bins,
                v: cdm.v(),
                spec,
            })
        })
        .collect()
}

/// Compute the right-tailed p-value for a (query, DB) pair given pre-built
/// spectra for both sides (produced by [`build_db_spectra`]).
///
/// ### Algorithm
/// 1. For each chromosome present in both `q_spectra` and `db_spectra`,
///    accumulate the cross-correlation: `C_total[k] += G_Q[k] · conj(G_B[k])`.
///    This is a pure O(N) dot product — no FFTs happen here.
/// 2. Scan `C_total` in ascending-frequency order, accumulating
///    `|C_total[k]|²`.  Stop at the smallest M such that the cumulative power
///    reaches `variance_threshold × total_power`.  Truncate to M coefficients.
/// 3. M-length complex IFFT (via the thread-local planner) → M samples `c[s]`
///    of total overlap vs a circular shift of `s` positions.
/// 4. Empirical p-value = fraction of the M samples with `c[s] ≥ c[0]`.
///    `c[0]` (shift = 0) corresponds to the actual observed overlap.
/// 5. `observed_bins` is the Parseval-normalised sum
///    Σ_chrom Re(Σ_k G_Q[k]·conj(G_B[k])) / N_chrom.
///
/// `variance_threshold` ∈ (0, 1]: fraction of total cross-correlation power
/// to retain.  Use [`DEFAULT_VARIANCE_THRESHOLD`] for a sensible default.
/// Lower values (e.g. 0.90) smooth the null distribution; 1.0 uses all
/// coefficients.
///
/// Returns `(observed_bins, p_value)` or `None` if no chromosomes are shared.
pub fn compute_pvalue_cached(
    q_spectra: &[ChromDbSpec],
    db_spectra: &[ChromDbSpec],
    variance_threshold: f64,
) -> Option<(f64, f64)> {
    // Full spectrum length = max shared-chromosome spectrum length.
    let m_full = db_spectra
        .iter()
        .filter(|db_cs| q_spectra.iter().any(|s| s.chrom == db_cs.chrom))
        .map(|db_cs| db_cs.spec.len())
        .max()?;

    let mut c_total: Vec<Complex<f32>> = vec![Complex::new(0.0, 0.0); m_full];
    let mut observed = 0.0f64;

    for db_cs in db_spectra {
        let q_spec = match q_spectra.iter().find(|s| s.chrom == db_cs.chrom) {
            Some(s) => &s.spec,
            None => continue,
        };

        let m_pair = m_full.min(db_cs.spec.len()).min(q_spec.len());

        // O(N) cross-correlation accumulation — no FFT.
        for k in 0..m_pair {
            c_total[k] += q_spec[k] * db_cs.spec[k].conj();
        }

        // Parseval-normalised observed overlap for this chromosome (in bins).
        let dot: f32 = q_spec[..m_pair]
            .iter()
            .zip(db_cs.spec[..m_pair].iter())
            .map(|(q, b)| (q * b.conj()).re)
            .sum();
        observed += (dot / db_cs.n as f32) as f64;
    }

    // Variance-adaptive frequency cutoff: find the smallest M such that
    // Σ_{k<M} |C_total[k]|² ≥ variance_threshold × total_power, then
    // truncate. This keeps low-frequency structure and discards high-
    // frequency noise, acting as a smoothing regulariser.
    // Fast path: threshold ≥ 1 means "use everything" — skip the scan.
    let m = if variance_threshold >= 1.0 {
        m_full
    } else {
        let total_power: f32 = c_total.iter().map(|c| c.norm_sqr()).sum();
        let target_power = total_power * variance_threshold as f32;
        let mut cumulative_power = 0.0f32;
        let mut cutoff = m_full;
        for (k, c) in c_total.iter().enumerate() {
            cumulative_power += c.norm_sqr();
            if cumulative_power >= target_power {
                cutoff = k + 1;
                break;
            }
        }
        cutoff
    };
    c_total.truncate(m);

    // M-length complex IFFT via the thread-local planner (plan cache hit).
    // c_total[s].re is proportional to the total overlap at shift s across
    // all shared chromosomes. c_total[0].re is the actual observed overlap.
    let ifft = COMPLEX_PLANNER.with(|p| p.borrow_mut().plan_fft_inverse(m));
    ifft.process(&mut c_total);

    let obs_ref = c_total[0].re;
    let count_geq = c_total.iter().filter(|x| x.re >= obs_ref).count();
    let p_value = count_geq as f64 / m as f64;

    Some((observed, p_value))
}

/// Convenience: compute p-value directly from two BedMaps (builds spectra inline).
/// Uses [`DEFAULT_VARIANCE_THRESHOLD`].
/// Use [`build_db_spectra`] + [`compute_pvalue_cached`] when either file is
/// compared against multiple partners (amortises the O(N log N) forward FFT).
pub fn compute_pvalue(query_bed: &BedMap, db_dm: &DepthMap) -> Option<(f64, f64)> {
    let db_spectra = build_db_spectra(db_dm);
    let q_dm = DepthMap::build(query_bed);
    let q_spectra = build_db_spectra(&q_dm);
    compute_pvalue_cached(&q_spectra, &db_spectra, DEFAULT_VARIANCE_THRESHOLD)
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

// ── Private helpers ───────────────────────────────────────────────────────────

/// Convert intervals to raw bp position vectors for fractional spike distribution.
///
/// Positions are clamped to `[0, n_bins × BIN_SIZE]`; the fractional
/// bin-split is applied later in [`chrom_coverage_spectrum`].
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

/// Build a real-valued coverage array of length `n = cdm.n()` via prefix sum
/// of the impulse train, optionally mask it, then forward-FFT it with
/// `realfft`, returning all n/2+1 complex coefficients.
///
/// When `mask` is `Some(m)`, bins whose entry in `m` is `false` are zeroed
/// after the prefix sum but before the FFT.  Bins beyond `m.len()` are already
/// zero-padded by `make_input_vec`.
///
/// Cost: O(N) build + O(N) mask application + O(N log N) FFT.
fn chrom_coverage_spectrum(cdm: &ChromDepthMap, mask: Option<&[bool]>) -> Vec<Complex<f32>> {
    let n = cdm.n();

    let fft = REAL_PLANNER.with(|p| p.borrow_mut().plan_fft_forward(n));
    let mut indata = fft.make_input_vec(); // vec![0f32; n]

    // Build fractionally-weighted impulse train (anti-aliased derivative of coverage).
    // Each spike at raw position x is split between bins floor(x/BIN) and floor(x/BIN)+1
    // with weights (1−frac) and frac, where frac = (x % BIN_SIZE) / BIN_SIZE.
    // For sub-bin intervals this correctly assigns fractional coverage rather than
    // cancelling to zero as the old integer-bin scheme did.
    for &x in &cdm.pos_spikes {
        let bin = (x / BIN_SIZE) as usize;
        let frac = (x % BIN_SIZE) as f32 / BIN_SIZE as f32;
        indata[bin.min(n - 1)] += 1.0 - frac;
        if frac > 0.0 && bin + 1 < n {
            indata[bin + 1] += frac;
        }
    }
    for &x in &cdm.neg_spikes {
        let bin = (x / BIN_SIZE) as usize;
        let frac = (x % BIN_SIZE) as f32 / BIN_SIZE as f32;
        indata[bin.min(n - 1)] -= 1.0 - frac;
        if frac > 0.0 && bin + 1 < n {
            indata[bin + 1] -= frac;
        }
    }
    // Prefix sum: impulse train → coverage g[j].
    for i in 1..n {
        indata[i] += indata[i - 1];
    }

    // Apply positional filter: zero out non-whitelisted bins.
    if let Some(mask) = mask {
        for (i, &allowed) in mask.iter().enumerate() {
            if !allowed {
                indata[i] = 0.0;
            }
        }
        // Bins beyond mask.len() are already 0.0 (zero-padded by make_input_vec).
    }

    let mut spectrum = fft.make_output_vec(); // complex, length n/2+1
    fft.process(&mut indata, &mut spectrum).unwrap();

    spectrum
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_pvalue_identical() {
        let mut bed = BedMap::new();
        // 10 Mb interval = 100_000 bins on chr22 (n=524288).
        // Querying with the same interval as the DB: observed overlap = 100_000 bins.
        // Under random shifts, only shift=0 achieves that maximum overlap, so the
        // empirical p-value should be very small (≈ 1/N).
        bed.insert("chr22".to_string(), vec![(10_000_000, 20_000_000)]);

        let dm = DepthMap::build(&bed);
        let q_dm = DepthMap::build(&bed);
        let db_spectra = build_db_spectra(&dm);
        let q_spectra = build_db_spectra(&q_dm);
        // Use threshold=1.0 (full spectrum) to get enough IFFT samples for a
        // tight empirical p-value — the default threshold is tuned for biology,
        // not for this synthetic single-peak test.
        let (observed, p_value) = compute_pvalue_cached(&q_spectra, &db_spectra, 1.0).unwrap();
        assert!(observed > 50_000.0, "observed={observed}");
        assert!(p_value < 0.1, "p_value={p_value}");
    }

    #[test]
    fn test_filter_mask_excludes_unmatched_chrom() {
        // DB and query both on chr22; filter only covers chr21.
        // With filter active, no chromosomes are shared → None.
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(10_000_000, 11_000_000)]);
        let dm = DepthMap::build(&bed);

        let mut filter_bed = BedMap::new();
        filter_bed.insert("chr21".to_string(), vec![(0, 46_709_983)]);
        let filter = FilterMask::build(&filter_bed, FilterMode::Whitelist);

        let db_spectra = build_db_spectra_with_filter(&dm, Some(&filter));
        let q_spectra = build_db_spectra_with_filter(&dm, Some(&filter));

        // chr22 excluded by filter → no spectra → p-value is None.
        assert!(db_spectra.is_empty());
        let result = compute_pvalue_cached(&q_spectra, &db_spectra, 1.0);
        assert!(result.is_none());
    }

    #[test]
    fn test_filter_mask_confines_signal() {
        // Query and DB identical on chr22; filter covers the whole chromosome.
        // Should still return a low p-value (signal is preserved).
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(10_000_000, 20_000_000)]);
        let dm = DepthMap::build(&bed);

        // Full chr22 whitelist (all bins accessible).
        let mut filter_bed = BedMap::new();
        filter_bed.insert("chr22".to_string(), vec![(0, 50_818_468)]);
        let filter = FilterMask::build(&filter_bed, FilterMode::Whitelist);

        let db_spectra = build_db_spectra_with_filter(&dm, Some(&filter));
        let q_spectra = build_db_spectra_with_filter(&dm, Some(&filter));
        let (observed, p_value) = compute_pvalue_cached(&q_spectra, &db_spectra, 1.0).unwrap();

        assert!(observed > 50_000.0, "observed={observed}");
        assert!(p_value < 0.1, "p_value={p_value}");
    }

    #[test]
    fn test_blacklist_excludes_signal() {
        // Query and DB both on chr22; blacklist covers the entire signal region.
        // After blacklisting, the coverage vectors are all-zero → no cross-
        // correlation power → observed overlap ~0 and p-value should be high.
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(10_000_000, 20_000_000)]);
        let dm = DepthMap::build(&bed);

        // Blacklist covers exactly the signal region (plus some slack).
        let mut bl_bed = BedMap::new();
        bl_bed.insert("chr22".to_string(), vec![(0, 50_818_468)]);
        let filter = FilterMask::build(&bl_bed, FilterMode::Blacklist);

        let db_spectra = build_db_spectra_with_filter(&dm, Some(&filter));
        let q_spectra = build_db_spectra_with_filter(&dm, Some(&filter));
        // chr22 is present in the blacklist mask (all bins blocked), so it is
        // not excluded from spectra — but all bins are zero, so observed ≈ 0.
        let (observed, _p_value) = compute_pvalue_cached(&q_spectra, &db_spectra, 1.0).unwrap();
        assert!(observed.abs() < 1.0, "observed should be ~0 after full blacklist, got {observed}");
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

        let (obs_cached, pv_cached) =
            compute_pvalue_cached(&q_spectra, &db_spectra, DEFAULT_VARIANCE_THRESHOLD).unwrap();
        let (obs_direct, pv_direct) = compute_pvalue(&bed, &dm).unwrap();

        assert!((obs_cached - obs_direct).abs() < 1.0, "obs mismatch");
        assert!((pv_cached - pv_direct).abs() < 0.05, "pv mismatch");
    }

    /// Verify fractional spike distribution for a sub-bin interval.
    ///
    /// Interval [20, 80] fits entirely within bin 0 (100 bp bins).
    /// Expected coverage in fractional bins: (80 − 20) / 100 = 0.60.
    ///
    /// Spike maths:
    ///   start=20 → bin 0 += 0.80, bin 1 += 0.20
    ///   end=80   → bin 0 -= 0.20, bin 1 -= 0.80
    ///   net:       bin 0 = +0.60, bin 1 = −0.60
    ///   prefix sum: g[0] = 0.60, g[1..] = 0.0
    #[test]
    fn test_sub_bin_fractional_coverage() {
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(20, 80)]);
        let dm = DepthMap::build(&bed);
        let chrom = &dm.chroms[0];

        // total_cov should equal (80 − 20) / 100 = 0.6 bins.
        let tc = chrom.total_cov();
        assert!((tc - 0.6).abs() < 1e-9, "total_cov={tc}");

        // Verify by checking raw spike storage (raw bp positions).
        assert_eq!(chrom.pos_spikes, vec![20u32]);
        assert_eq!(chrom.neg_spikes, vec![80u32]);
    }
}
