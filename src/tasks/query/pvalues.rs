use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use rayon::prelude::*;

use crate::fourier::{
    ChromMoments, DEFAULT_MOMENTS_EPS, DepthMap, FilterMask, QueryChromData,
    build_chrom_moments_with_filter, build_depth_moments, build_query_chrom_data,
    compute_analytic_stats, mean_interval_bins, parse_bed_as_map,
};
use crate::io::{MappedDepthStore, MappedMomentStore};
use crate::matrix::SparseMatrix;

use super::PValueResult;

/// For every non-zero (q_sid, d_sid) pair in `counts`, compute analytic NB stats.
///
/// Phase A: Load DB DepthMaps from `depthmap.rkyv` and, when available and no
/// filter is active, pre-stored moments from `momentmap.rkyv`.
///
/// Phase B: Build query interval data once per query file.
///
/// Phase C: For each DB file obtain moments (from store or computed on-the-fly
/// via FFT), then score every overlapping query.
pub(super) fn compute_analytic_pvalues(
    counts: &SparseMatrix,
    query_file_paths: &[PathBuf],
    db_path: &Path,
    filter: Option<&FilterMask>,
) -> Vec<PValueResult> {
    let mut by_db: HashMap<u32, Vec<usize>> = HashMap::new();
    for (q_sid, row) in counts.outer_iterator().enumerate() {
        for (d_sid, &cnt) in row.iter() {
            if cnt > 0 {
                by_db.entry(d_sid as u32).or_default().push(q_sid);
            }
        }
    }
    if by_db.is_empty() {
        return Vec::new();
    }

    let needed_q_sids: HashSet<usize> = by_db.values().flatten().copied().collect();

    // Phase A: open depth-map store (always needed) and moment store (optional).
    let store_path = db_path.join("depthmap.rkyv");
    let store = match MappedDepthStore::open(&store_path) {
        Ok(s) => s,
        Err(_) => return Vec::new(),
    };
    let moment_store_path = db_path.join("momentmap.rkyv");
    let moment_store = MappedMomentStore::open(&moment_store_path).ok();

    let db_depthmaps: HashMap<u32, DepthMap> = by_db
        .par_iter()
        .filter_map(|(&d_sid, _)| Some((d_sid, store.get(d_sid)?)))
        .collect();

    // Phase B: build query interval data once per query file.
    let query_cov_data: HashMap<usize, (Vec<QueryChromData>, f64)> = needed_q_sids
        .par_iter()
        .filter_map(|&q_sid| {
            let path = query_file_paths.get(q_sid)?;
            let bed = parse_bed_as_map(path).ok()?;
            let q_data = build_query_chrom_data(&bed, filter);
            let mean_iv = mean_interval_bins(&q_data);
            Some((q_sid, (q_data, mean_iv)))
        })
        .collect();

    // Phase C: obtain moments once per DB source, then score every query.
    //
    // Hot path: `moment_store.get_sid(d_sid)` returns a zerocopy handle backed
    // by the mmap — no f32→f64 conversion, no Vec allocation.  The chrom→index
    // map inside it is O(chroms) and built once per d_sid.
    //
    // Fallback: when the store is absent or a positional filter is active,
    // moments are computed on-the-fly via FFT (O(N log N) per DB file).  A
    // HashMap is built once per d_sid so chrom lookup is O(1) in the inner loop.
    by_db
        .par_iter()
        .flat_map(|(&d_sid, q_sids)| -> Vec<PValueResult> {
            let dm = match db_depthmaps.get(&d_sid) {
                Some(d) => d,
                None => return Vec::new(),
            };

            // Zerocopy stored moments (no filter, store present).
            let sid_moments = if filter.is_none() {
                moment_store.as_ref().and_then(|ms| ms.get_sid(d_sid))
            } else {
                None
            };

            // Fallback: compute moments once per d_sid when stored data is unavailable.
            let fallback: Option<Vec<ChromMoments>> = if sid_moments.is_none() {
                Some(if let Some(f) = filter {
                    build_moments_with_filter(dm, f)
                } else {
                    build_depth_moments(dm, DEFAULT_MOMENTS_EPS)
                })
            } else {
                None
            };

            // Pre-index fallback by chrom so the inner loop stays O(1) per interval.
            let fallback_map: Option<HashMap<&str, &ChromMoments>> =
                fallback.as_ref().map(|fb| fb.iter().map(|m| (m.chrom.as_str(), m)).collect());

            q_sids
                .iter()
                .filter_map(|&q_sid| {
                    let (q_data, mean_iv) = query_cov_data.get(&q_sid)?;
                    let sweep_count = counts.get(q_sid, d_sid as usize).copied().unwrap_or(0);
                    let observed_bins = sweep_count as f64 * mean_iv;

                    let lookup = |chrom: &str, l: f64| -> Option<(f64, f64)> {
                        if let Some(ref sm) = sid_moments {
                            sm.lookup(chrom, l)
                        } else if let Some(ref fm) = fallback_map {
                            fm.get(chrom)?.lookup(l)
                        } else {
                            None
                        }
                    };

                    let (p_value, llr) =
                        compute_analytic_stats(q_data, lookup, observed_bins)?;
                    Some(PValueResult {
                        query_id: q_sid,
                        db_sid: d_sid,
                        observed_bins,
                        p_value,
                        llr,
                    })
                })
                .collect()
        })
        .collect()
}

fn build_moments_with_filter(dm: &DepthMap, filter: &FilterMask) -> Vec<ChromMoments> {
    dm.chroms
        .iter()
        .filter_map(|cdm| {
            let mask = filter.get(&cdm.chrom)?;
            Some(build_chrom_moments_with_filter(cdm, mask, DEFAULT_MOMENTS_EPS))
        })
        .collect()
}
