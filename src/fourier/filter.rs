use std::collections::HashMap;

use super::genome::{BedMap, BIN_SIZE, hg38_chrom_sizes};

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fourier::{
        ChromMoments, DEFAULT_MOMENTS_EPS, DepthMap, build_chrom_moments_with_filter,
        build_query_chrom_data, compute_analytic_stats,
    };

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
        let lookup = |chrom: &str, l: f64| {
            db_moments.iter().find(|m| m.chrom == chrom)?.lookup(l)
        };
        assert!(compute_analytic_stats(&q_data, lookup, 0.0).is_none());
    }
}
