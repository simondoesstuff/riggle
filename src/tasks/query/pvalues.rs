use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use rayon::prelude::*;

use crate::fourier::{QueryChromData, build_query_chrom_data, compute_analytic_stats,
    parse_bed_as_map};
use crate::io::MappedMomentStore;
use crate::matrix::SparseMatrix;

use super::PValueResult;

/// For every non-zero (q_sid, d_sid) pair in `counts`, compute analytic NB stats.
///
/// Phase A: Open `momentmap.rkyv`; return empty if absent.
/// Phase B: Build query interval data once per query file.
/// Phase C: For each DB source, look up pre-stored moments (O(1) per interval),
///          accumulate μ_null / σ²_null, fit NB, score, return p-value and LLR.
///
/// `overlap` supplies the exact base-pair overlap in 100 bp bins for each
/// (q_sid, d_sid) pair, computed during the sweep phase.
pub(super) fn compute_analytic_pvalues(
    counts: &SparseMatrix,
    query_file_paths: &[PathBuf],
    db_path: &Path,
    overlap: &HashMap<(usize, usize), f32>,
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

    // Phase A: open moment store; abort if unavailable.
    let moment_store_path = db_path.join("momentmap.rkyv");
    let moment_store = match MappedMomentStore::open(&moment_store_path) {
        Ok(s) => s,
        Err(_) => return Vec::new(),
    };

    let needed_q_sids: HashSet<usize> = by_db.values().flatten().copied().collect();

    // Phase B: build query interval data once per query file.
    let query_cov_data: HashMap<usize, Vec<QueryChromData>> = needed_q_sids
        .par_iter()
        .filter_map(|&q_sid| {
            let path = query_file_paths.get(q_sid)?;
            let bed = parse_bed_as_map(path).ok()?;
            let q_data = build_query_chrom_data(&bed);
            Some((q_sid, q_data))
        })
        .collect();

    // Phase C: score each (d_sid, q_sid) pair using stored moments.
    by_db
        .par_iter()
        .flat_map(|(&d_sid, q_sids)| -> Vec<PValueResult> {
            let sid_moments = match moment_store.get_sid(d_sid) {
                Some(m) => m,
                None => return Vec::new(),
            };

            q_sids
                .iter()
                .filter_map(|&q_sid| {
                    let q_data = query_cov_data.get(&q_sid)?;
                    // Exact overlap in bins from the sweep phase.
                    let observed_bins =
                        overlap.get(&(q_sid, d_sid as usize)).copied().unwrap_or(0.0) as f64;

                    let lookup = |chrom: &str, l: f64| sid_moments.lookup(chrom, l);
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
