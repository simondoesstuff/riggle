use crate::core::Interval;
use crate::fourier::BIN_SIZE;
use crate::io::MappedJumpTable;
use crate::matrix::{DenseMatrix, OverlapMatrix};

/// Collect every overlapping (query, DB) interval pair and its intersection length.
///
/// Same jump-table positioning as [`query_sweep`]; returns owned pairs rather than
/// accumulating into matrices.  Used by `--intervals` mode.
pub fn query_sweep_pairs(
    db_layer: &[Interval],
    layer_max_size: u32,
    jump_table: &MappedJumpTable,
    query_block: &[Interval],
) -> Vec<(Interval, Interval, u32)> {
    let mut pairs = Vec::new();
    if query_block.is_empty() || db_layer.is_empty() {
        return pairs;
    }
    for q in query_block {
        let lo = q.start.saturating_sub(layer_max_size);
        let approx = jump_table.lookup(lo).min(db_layer.len());
        let start_idx = approx + db_layer[approx..].partition_point(|d| d.start < lo);
        for d in &db_layer[start_idx..] {
            if d.start >= q.end {
                break;
            }
            if d.end > q.start {
                let intersection_bp = q.end.min(d.end) - q.start.max(d.start);
                pairs.push((*q, *d, intersection_bp));
            }
        }
    }
    pairs
}

/// Per-query scan over a flat sorted layer, with O(1) cold-start via jump table.
///
/// For each query interval, uses the jump table to land near the first database
/// interval that could overlap it, then scans forward collecting hits.  Every
/// query is treated independently, as if all query intervals were well separated.
///
/// ### Inputs
///
/// - `db_layer`: flat sorted `&[Interval]` from a `MappedLayer`.
/// - `layer_max_size`: the exclusive upper bound on interval sizes in this
///   layer (`LayerConfig::layer_max_size(K)`). Used to bound the search window.
/// - `jump_table`: O(1) index for this layer.  `lookup(lo)` returns a
///   conservative lower-bound index; a short forward scan completes positioning.
/// - `query_block`: sorted `&[Interval]` where `sid` is the row index in
///   `results` (Q_SID).
/// - `results`: `[Q_sids][D_sids]` dense count accumulator, mutated in place.
/// - `overlap`: parallel f32 accumulator; each hit adds the exact overlap in
///   100 bp bins — `(min(q.end, d.end) − max(q.start, d.start)) / BIN_SIZE`.
pub fn query_sweep(
    db_layer: &[Interval],
    layer_max_size: u32,
    jump_table: &MappedJumpTable,
    query_block: &[Interval],
    results: &mut DenseMatrix,
    overlap: &mut OverlapMatrix,
) {
    if query_block.is_empty() || db_layer.is_empty() {
        return;
    }

    for q in query_block {
        // d overlaps q iff d.start < q.end && d.end > q.start.
        // Since d.end <= d.start + layer_max_size, d.end > q.start implies
        // d.start > q.start - layer_max_size.  Use the jump table for an O(1)
        // cold start to that lower bound, then scan forward within one tile.
        let lo = q.start.saturating_sub(layer_max_size);
        let approx = jump_table.lookup(lo).min(db_layer.len());
        let start_idx = approx + db_layer[approx..].partition_point(|d| d.start < lo);
        let q_sid = q.sid as usize;
        for d in &db_layer[start_idx..] {
            if d.start >= q.end {
                break;
            }
            if d.end > q.start {
                let d_sid = d.sid as usize;
                results.add(q_sid, d_sid, 1);
                let overlap_bins =
                    (q.end.min(d.end) - q.start.max(d.start)) as f32 / BIN_SIZE as f32;
                overlap.add(q_sid, d_sid, overlap_bins);
            }
        }
    }
}

/// Like [`query_sweep`] but writes to a sub-window of Q_SID rows.
///
/// `q_sid_base` is the global Q_SID of row 0 in `results`.  Each query's
/// local row index is `q.sid - q_sid_base`, so `results` need only be
/// `window_rows × num_sources` rather than the full `num_queries × num_sources`.
/// Keeping `results` small (a few hundred KB) lets it stay in L2/L3 cache
/// rather than forcing DRAM traffic on every write.
///
/// The caller must guarantee that all intervals in `query_block` have
/// `q_sid_base <= q.sid < q_sid_base + results.num_rows()`.
pub fn query_sweep_windowed(
    db_layer: &[Interval],
    layer_max_size: u32,
    jump_table: &MappedJumpTable,
    query_block: &[Interval],
    q_sid_base: u32,
    results: &mut DenseMatrix,
    overlap: &mut OverlapMatrix,
) {
    if query_block.is_empty() || db_layer.is_empty() {
        return;
    }
    for q in query_block {
        let lo = q.start.saturating_sub(layer_max_size);
        let approx = jump_table.lookup(lo).min(db_layer.len());
        let start_idx = approx + db_layer[approx..].partition_point(|d| d.start < lo);
        let local_q_sid = (q.sid - q_sid_base) as usize;
        for d in &db_layer[start_idx..] {
            if d.start >= q.end {
                break;
            }
            if d.end > q.start {
                let d_sid = d.sid as usize;
                results.add(local_q_sid, d_sid, 1);
                let overlap_bins =
                    (q.end.min(d.end) - q.start.max(d.start)) as f32 / BIN_SIZE as f32;
                overlap.add(local_q_sid, d_sid, overlap_bins);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::{build_jump_table, write_jump_table};
    use tempfile::NamedTempFile;

    fn iv(start: u32, end: u32, sid: u32) -> Interval {
        Interval::new(start, end, sid)
    }

    fn make_jt(db: &[Interval], layer_max_size: u32) -> (NamedTempFile, MappedJumpTable) {
        let tile_size = (layer_max_size / 2).max(1);
        let f = NamedTempFile::new().unwrap();
        let tbl = build_jump_table(db, tile_size);
        write_jump_table(f.path(), &tbl).unwrap();
        let jt = MappedJumpTable::open(f.path(), tile_size).unwrap();
        (f, jt)
    }

    fn sweep(
        db: &[Interval],
        layer_max_size: u32,
        query: &[Interval],
        nrows: usize,
        ncols: usize,
    ) -> (DenseMatrix, OverlapMatrix) {
        let (_f, jt) = make_jt(db, layer_max_size);
        let mut results = DenseMatrix::new(nrows, ncols);
        let mut overlap = OverlapMatrix::new(nrows, ncols);
        query_sweep(db, layer_max_size, &jt, query, &mut results, &mut overlap);
        (results, overlap)
    }

    #[test]
    fn test_basic_overlap() {
        // Q=[40,160) overlaps D=[50,150) by 100 bp = 1 bin
        let (results, overlap) = sweep(&[iv(50, 150, 0)], 200, &[iv(40, 160, 0)], 1, 1);
        assert_eq!(results.get(0, 0), 1);
        assert!((overlap.get(0, 0) - 1.0).abs() < 1e-5, "expected 1.0 bin, got {}", overlap.get(0, 0));
    }

    #[test]
    fn test_no_overlap() {
        let (results, overlap) = sweep(&[iv(50, 100, 0)], 200, &[iv(200, 300, 0)], 1, 1);
        assert_eq!(results.get(0, 0), 0);
        assert_eq!(overlap.get(0, 0), 0.0);
    }

    #[test]
    fn test_multiple_queries_and_db() {
        // D0=[25,75), D1=[75,125); Q0=[0,50), Q1=[50,100)
        let db = vec![iv(25, 75, 0), iv(75, 125, 1)];
        let query = vec![iv(0, 50, 0), iv(50, 100, 1)];
        let (results, overlap) = sweep(&db, 200, &query, 2, 2);

        assert_eq!(results.get(0, 0), 1);
        assert_eq!(results.get(0, 1), 0);
        assert_eq!(results.get(1, 0), 1);
        assert_eq!(results.get(1, 1), 1);
        // Q0∩D0 = [25,50) = 25 bp = 0.25 bins
        assert!((overlap.get(0, 0) - 0.25).abs() < 1e-5, "Q0∩D0={}", overlap.get(0, 0));
        // Q1∩D0 = [50,75) = 25 bp = 0.25 bins
        assert!((overlap.get(1, 0) - 0.25).abs() < 1e-5, "Q1∩D0={}", overlap.get(1, 0));
        // Q1∩D1 = [75,100) = 25 bp = 0.25 bins
        assert!((overlap.get(1, 1) - 0.25).abs() < 1e-5, "Q1∩D1={}", overlap.get(1, 1));
    }

    #[test]
    fn test_dead_zone_skipped_via_jump_table() {
        let db = vec![iv(0, 50, 0), iv(10_000, 10_050, 1)];
        let query = vec![iv(10_000, 10_100, 0)];
        let (results, overlap) = sweep(&db, 100, &query, 1, 2);

        assert_eq!(results.get(0, 0), 0); // D0 far away, no overlap
        assert_eq!(results.get(0, 1), 1); // D1 overlaps
        assert_eq!(overlap.get(0, 0), 0.0);
        // Q∩D1 = [10000,10050) = 50 bp = 0.5 bins
        assert!((overlap.get(0, 1) - 0.5).abs() < 1e-5, "Q∩D1={}", overlap.get(0, 1));
    }

    #[test]
    fn test_touching_intervals_not_overlapping() {
        let (results, overlap) = sweep(&[iv(0, 50, 0)], 200, &[iv(50, 100, 0)], 1, 1);
        assert_eq!(results.get(0, 0), 0);
        assert_eq!(overlap.get(0, 0), 0.0);
    }

    #[test]
    fn test_empty_db() {
        let (results, overlap) = sweep(&[], 200, &[iv(0, 100, 0)], 1, 1);
        assert_eq!(results.get(0, 0), 0);
        assert_eq!(overlap.get(0, 0), 0.0);
    }

    #[test]
    fn test_empty_query() {
        let (results, overlap) = sweep(&[iv(0, 100, 0)], 200, &[], 1, 1);
        assert_eq!(results.get(0, 0), 0);
        assert_eq!(overlap.get(0, 0), 0.0);
    }
}
