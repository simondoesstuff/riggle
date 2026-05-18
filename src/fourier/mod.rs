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
//! # Frequency truncation
//! The absolute error from keeping only M frequencies is bounded by
//! C / M  where  C ≤ V_Q · V_B  and  V = Σ|impulse| = 2K.
//! We use M = ceil(√(V_Q · V_B / ε)), rounded to the next power of two,
//! capped at M_MAX = 2^25 (~33 M).  For typical ChIP-seq data M ≪ N,
//! so the real FFT still runs far less than the full N-point transform.
//!
//! # Null model
//! All chromosomes are shifted by the same sample index s simultaneously.
//! This is equivalent to adding cross-correlation spectra before the IFFT
//! (linearity of IFFT) and gives a single M-sample null distribution for
//! the combined multi-chromosome overlap.

use std::cell::RefCell;
use std::collections::HashMap;
use std::fs;
use std::io::{self, Read, Write};
use std::path::Path;

use rayon::prelude::*;
use realfft::RealFftPlanner;
use rustfft::{FftPlanner, num_complex::Complex};

thread_local! {
    static REAL_PLANNER: RefCell<RealFftPlanner<f32>> = RefCell::new(RealFftPlanner::new());
    static COMPLEX_PLANNER: RefCell<FftPlanner<f32>> = RefCell::new(FftPlanner::new());
}

const BIN_SIZE: u32 = 100;

/// Absolute error tolerance in 100 bp overlap bins.
pub const DEFAULT_EPSILON: f64 = 1.0;

/// Hard upper cap on the number of frequencies to retain (must be a power of two).
/// ≈ 33 M — well above the typical genomic working range while bounding memory.
const M_MAX: usize = 1 << 25;
/// Minimum number of frequencies — ensures enough null samples even when V is tiny.
const M_MIN: usize = 64;

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

/// M = number of frequencies to retain; rounded up to the next power of two,
/// capped at M_MAX.
pub fn compute_m(v_query: f64, v_db: f64, epsilon: f64) -> usize {
    if v_query == 0.0 || v_db == 0.0 {
        return M_MIN;
    }
    ((v_query * v_db / epsilon).sqrt().ceil() as usize)
        .max(M_MIN)
        .min(M_MAX)
        .next_power_of_two()
        .min(M_MAX)
}

/// Build per-chromosome coverage spectra for a DB DepthMap.
/// `m` is the number of frequencies to retain; computed via [`compute_m`].
/// Runs in parallel over chromosomes.
pub fn build_db_spectra(dm: &DepthMap, m: usize) -> Vec<ChromDbSpec> {
    dm.chroms
        .par_iter()
        .map(|cdm| {
            let spec = chrom_coverage_spectrum(cdm, m);
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
/// `m` is the number of null samples (IFFT length); must be a power of two.
/// The per-chromosome spectra may have been built at a larger M — only the
/// first `m` coefficients are used.
///
/// ### Algorithm
/// 1. For each chromosome present in both `q_spectra` and `db_spectra`,
///    accumulate the cross-correlation: `C_total[k] += G_Q[k] · conj(G_B[k])`.
///    This is a pure O(M) dot product — no FFTs happen here.
/// 2. Single M-length complex IFFT of `C_total` (via the thread-local planner)
///    → M samples `c[s]` of total overlap vs a circular shift of s·N/M bins.
/// 3. Empirical p-value = fraction of the M samples with `c[s] ≥ c[0]`.
///    `c[0]` (shift = 0) corresponds to the actual observed overlap.
/// 4. `observed_bins` is the Parseval-normalised sum
///    Σ_chrom Re(Σ_k G_Q[k]·conj(G_B[k])) / N_chrom.
///
/// Both `q_spectra` and `db_spectra` are produced by [`build_db_spectra`].
/// Call [`build_db_spectra`] once per file and reuse across multiple comparisons
/// to amortise the O(N log N) forward FFT cost.
///
/// Returns `(observed_bins, p_value)` or `None` if no chromosomes are shared.
pub fn compute_pvalue_cached(
    q_spectra: &[ChromDbSpec],
    db_spectra: &[ChromDbSpec],
    m: usize,
) -> Option<(f64, f64)> {
    let mut c_total: Vec<Complex<f32>> = vec![Complex::new(0.0, 0.0); m];
    let mut observed = 0.0f64;
    let mut has_data = false;

    for db_cs in db_spectra {
        let q_spec = match q_spectra.iter().find(|s| s.chrom == db_cs.chrom) {
            Some(s) => &s.spec,
            None => continue,
        };

        let m_pair = m.min(db_cs.spec.len()).min(q_spec.len());

        // O(M) cross-correlation accumulation — no FFT.
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

        has_data = true;
    }

    if !has_data {
        return None;
    }

    // Single M-length complex IFFT via the thread-local planner (plan cache hit).
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
/// Use [`build_db_spectra`] + [`compute_pvalue_cached`] when either file is
/// compared against multiple partners (amortises the O(N log N) forward FFT).
pub fn compute_pvalue(query_bed: &BedMap, db_dm: &DepthMap, epsilon: f64) -> Option<(f64, f64)> {
    let v_q = bed_map_v(query_bed);
    let v_db = db_dm.total_v();
    let m = compute_m(v_q, v_db, epsilon);
    let db_spectra = build_db_spectra(db_dm, m);
    let q_dm = DepthMap::build(query_bed);
    let q_spectra = build_db_spectra(&q_dm, m);
    compute_pvalue_cached(&q_spectra, &db_spectra, m)
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
/// of the impulse train, then forward-FFT it with `realfft`, returning the
/// first `m` complex coefficients (capped at n/2+1).
///
/// Cost: O(N) build + O(N log N) FFT — fast even for large M because rustfft
/// uses SIMD and the full FFT benefits from the efficient plan cache.
fn chrom_coverage_spectrum(cdm: &ChromDepthMap, m: usize) -> Vec<Complex<f32>> {
    let n = cdm.n();
    let m_cap = m.min(n / 2 + 1);

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

    spectrum[..m_cap].to_vec()
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
        // Under random shifts, only shift=0 achieves that maximum overlap, so the
        // empirical p-value should be very small (≈ 1/M).
        bed.insert("chr22".to_string(), vec![(10_000_000, 20_000_000)]);

        let dm = DepthMap::build(&bed);
        let (observed, p_value) = compute_pvalue(&bed, &dm, DEFAULT_EPSILON).unwrap();
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

        assert!(compute_pvalue(&query, &dm, DEFAULT_EPSILON).is_none());
    }

    #[test]
    fn test_compute_m_capped() {
        // Very large V should not exceed M_MAX
        let m = compute_m(1e8, 1e8, 1.0);
        assert!(m <= M_MAX);
        assert!(m.is_power_of_two());
    }

    #[test]
    fn test_cached_matches_direct() {
        let mut bed = BedMap::new();
        bed.insert(
            "chr22".to_string(),
            vec![(10_000_000, 10_001_000), (20_000_000, 20_002_000)],
        );
        let dm = DepthMap::build(&bed);
        let v_q = bed_map_v(&bed);
        let v_db = dm.total_v();
        let m = compute_m(v_q, v_db, DEFAULT_EPSILON);
        let db_spectra = build_db_spectra(&dm, m);
        let q_dm = DepthMap::build(&bed);
        let q_spectra = build_db_spectra(&q_dm, m);

        let (obs_cached, pv_cached) = compute_pvalue_cached(&q_spectra, &db_spectra, m).unwrap();
        let (obs_direct, pv_direct) = compute_pvalue(&bed, &dm, DEFAULT_EPSILON).unwrap();

        assert!((obs_cached - obs_direct).abs() < 1.0, "obs mismatch");
        assert!((pv_cached - pv_direct).abs() < 0.05, "pv mismatch");
    }
}
