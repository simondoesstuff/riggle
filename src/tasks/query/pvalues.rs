use std::collections::HashMap;
use std::path::Path;

use rayon::prelude::*;

use crate::fourier::{QueryChromData, compute_analytic_stats};
use crate::io::MappedMomentStore;
use crate::matrix::SparseMatrix;

use super::PValueResult;

/// For every non-zero (q_sid, d_sid) pair in `counts`, compute analytic NB stats.
///
/// `query_chrom_data[q_sid]` holds the pre-built per-chromosome interval data for
/// each query file (built once during `parse_file_batch`; no re-parse at query time).
///
/// Phase A: Open `momentmap.rkyv`; return empty if absent.
/// Phase B: For each DB source with hits, look up pre-stored moments (O(1) per
///          unique block size), accumulate μ_null / σ²_null, fit NB, return p-value
///          and LLR.  Zero-overlap pairs are not emitted; callers treat absent
///          entries as p = 1.
pub(super) fn compute_analytic_pvalues(
    counts: &SparseMatrix,
    query_chrom_data: &[Vec<QueryChromData>],
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

    // Phase B: score each (d_sid, q_sid) pair using stored moments.
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
                    let q_data = query_chrom_data.get(q_sid)?;
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
