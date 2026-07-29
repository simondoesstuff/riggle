use std::collections::HashMap;
use std::path::Path;

use wide::f32x8;

use super::genome::{BedMap, BIN_SIZE, hg38_chrom_sizes};

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

pub fn bed_map_v(bed: &BedMap) -> f64 {
    bed.values().map(|v| 2 * v.len()).sum::<usize>() as f64
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

/// Exact base-pair overlap in 100 bp bins between two depth maps (dot product).
pub fn coverage_dot_product(a_dm: &DepthMap, b_dm: &DepthMap) -> f64 {
    let a_map: HashMap<&str, &ChromDepthMap> =
        a_dm.chroms.iter().map(|c| (c.chrom.as_str(), c)).collect();
    b_dm.chroms
        .iter()
        .filter_map(|b_cd| {
            let a_cd = a_map.get(b_cd.chrom.as_str())?;
            let a_cov = build_depth_signal(a_cd);
            let b_cov = build_depth_signal(b_cd);
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

/// Build the full per-bin coverage signal g[0..n_bins-1] with fractional binning.
///
/// The impulse train is assembled in-place; the prefix sum (cumsum) is then
/// applied using a SIMD two-pass scan (f32x8) for throughput.
pub(crate) fn build_depth_signal(cdm: &ChromDepthMap) -> Vec<f32> {
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

    simd_cumsum_f32(&mut g);
    g
}

/// Two-pass SIMD prefix sum (cumulative sum) for f32 slices.
///
/// Pass 1: local prefix sums within independent chunks of 8 (f32x8).
/// Pass 2: sequential prefix of chunk totals → per-chunk carry offsets.
/// Pass 3: add carry offsets to each chunk with a SIMD broadcast-add.
fn simd_cumsum_f32(g: &mut [f32]) {
    let n = g.len();
    if n <= 1 {
        return;
    }

    const W: usize = 8;
    let full_chunks = n / W;
    let tail = full_chunks * W;

    // Pass 1: within each chunk, compute local prefix sums (chunks are independent).
    for c in 0..full_chunks {
        let b = c * W;
        for i in 1..W {
            g[b + i] += g[b + i - 1];
        }
    }

    // Tail: sequential prefix sum (< W elements, no inter-chunk dependency yet).
    for i in tail + 1..n {
        g[i] += g[i - 1];
    }

    if full_chunks == 0 {
        return;
    }

    // Pass 2: prefix sum of chunk totals → carry for each chunk.
    // After pass 1, g[c*W + W - 1] = local sum of chunk c.
    let mut carry = 0.0f32;
    let mut carries = vec![0.0f32; full_chunks];
    for c in 0..full_chunks {
        carries[c] = carry;
        carry += g[c * W + W - 1];
    }
    // `carry` now equals the sum of all complete-chunk elements.

    // Pass 3: broadcast-add chunk offsets using f32x8 (skip chunk 0, offset = 0).
    for c in 1..full_chunks {
        let off = carries[c];
        let b = c * W;
        let v = f32x8::from([
            g[b],
            g[b + 1],
            g[b + 2],
            g[b + 3],
            g[b + 4],
            g[b + 5],
            g[b + 6],
            g[b + 7],
        ]);
        let out = (v + f32x8::splat(off)).to_array();
        g[b..b + W].copy_from_slice(&out);
    }

    // Add the complete-chunk carry to each tail element.
    if tail < n && carry != 0.0 {
        for i in tail..n {
            g[i] += carry;
        }
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sub_bin_fractional_coverage() {
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(20, 80)]);
        let dm = DepthMap::build(&bed);
        let chrom = &dm.chroms[0];

        let tc = chrom.total_cov();
        assert!((tc - 0.6).abs() < 1e-9, "total_cov={tc}");
    }

    #[test]
    fn test_simd_cumsum_correctness() {
        // Verify the two-pass SIMD cumsum matches a brute-force sequential cumsum.
        let original: Vec<f32> = (0..100).map(|i| i as f32 * 0.5 - 25.0).collect();

        let mut expected = original.clone();
        for i in 1..expected.len() {
            expected[i] += expected[i - 1];
        }

        let mut got = original.clone();
        simd_cumsum_f32(&mut got);

        for (i, (e, g)) in expected.iter().zip(got.iter()).enumerate() {
            assert!(
                (e - g).abs() < 1e-4,
                "mismatch at index {i}: expected {e}, got {g}"
            );
        }
    }

    #[test]
    fn test_simd_cumsum_non_multiple_of_8() {
        let original: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0];

        let mut expected = original.clone();
        for i in 1..expected.len() {
            expected[i] += expected[i - 1];
        }

        let mut got = original.clone();
        simd_cumsum_f32(&mut got);

        for (i, (e, g)) in expected.iter().zip(got.iter()).enumerate() {
            assert!((e - g).abs() < 1e-4, "index {i}: expected {e}, got {g}");
        }
    }
}
