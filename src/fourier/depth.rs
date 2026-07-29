use std::collections::HashMap;
use std::path::Path;

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

/// Build the full per-bin coverage signal g[0..n_bins-1] with fractional binning.
///
/// If `mask` is supplied, non-accessible bins are zeroed after the prefix sum.
pub(crate) fn build_depth_signal(cdm: &ChromDepthMap, mask: Option<&[bool]>) -> Vec<f32> {
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
}
