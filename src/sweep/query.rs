use crate::core::Interval;
use crate::fourier::BIN_SIZE;
use crate::io::layer::Position;
use crate::io::MappedJumpTable;
use crate::matrix::{DenseMatrix, OverlapMatrix};

/// Collect every overlapping (query, DB) interval pair and its intersection length.
///
/// Same jump-table positioning as [`query_sweep`]; returns owned pairs rather than
/// accumulating into matrices.  Used by `--intervals` mode.
pub fn query_sweep_pairs(
    db_positions: &[Position],
    db_sids: &[u32],
    layer_max_size: u32,
    jump_table: &MappedJumpTable,
    query_block: &[Interval],
) -> Vec<(Interval, Interval, u32)> {
    let mut pairs = Vec::new();
    if query_block.is_empty() || db_positions.is_empty() {
        return pairs;
    }
    for q in query_block {
        let lo = q.start.saturating_sub(layer_max_size);
        let approx = jump_table.lookup(lo).min(db_positions.len());
        let start_idx = approx + db_positions[approx..].partition_point(|p| p.start < lo);
        let pos_slice = &db_positions[start_idx..];
        let sid_slice = &db_sids[start_idx..];
        for (i, p) in pos_slice.iter().enumerate() {
            if p.start >= q.end {
                break;
            }
            if p.end > q.start {
                let d_sid = unsafe { *sid_slice.get_unchecked(i) };
                let intersection_bp = q.end.min(p.end) - q.start.max(p.start);
                let d = Interval::new(p.start, p.end, d_sid);
                pairs.push((*q, d, intersection_bp));
            }
        }
    }
    pairs
}

/// Per-query scan over a flat sorted layer, with O(1) cold-start via jump table.
pub fn query_sweep(
    db_positions: &[Position],
    db_sids: &[u32],
    layer_max_size: u32,
    jump_table: &MappedJumpTable,
    query_block: &[Interval],
    results: &mut DenseMatrix,
    overlap: &mut OverlapMatrix,
) {
    if query_block.is_empty() || db_positions.is_empty() {
        return;
    }

    for q in query_block {
        let lo = q.start.saturating_sub(layer_max_size);
        let approx = jump_table.lookup(lo).min(db_positions.len());
        let start_idx = approx + db_positions[approx..].partition_point(|p| p.start < lo);
        let q_sid = q.sid as usize;
        let pos_slice = &db_positions[start_idx..];
        let sid_slice = &db_sids[start_idx..];
        for (i, p) in pos_slice.iter().enumerate() {
            if p.start >= q.end {
                break;
            }
            if p.end > q.start {
                let d_sid = unsafe { *sid_slice.get_unchecked(i) } as usize;
                results.add(q_sid, d_sid, 1);
                let overlap_bins =
                    (q.end.min(p.end) - q.start.max(p.start)) as f32 / BIN_SIZE as f32;
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
pub fn query_sweep_windowed(
    db_positions: &[Position],
    db_sids: &[u32],
    layer_max_size: u32,
    jump_table: &MappedJumpTable,
    query_block: &[Interval],
    q_sid_base: u32,
    results: &mut DenseMatrix,
    overlap: &mut OverlapMatrix,
) {
    if query_block.is_empty() || db_positions.is_empty() {
        return;
    }
    for q in query_block {
        let lo = q.start.saturating_sub(layer_max_size);
        let approx = jump_table.lookup(lo).min(db_positions.len());
        let start_idx = approx + db_positions[approx..].partition_point(|p| p.start < lo);
        let local_q_sid = (q.sid - q_sid_base) as usize;
        let pos_slice = &db_positions[start_idx..];
        let sid_slice = &db_sids[start_idx..];
        for (i, p) in pos_slice.iter().enumerate() {
            if p.start >= q.end {
                break;
            }
            if p.end > q.start {
                let d_sid = unsafe { *sid_slice.get_unchecked(i) } as usize;
                results.add(local_q_sid, d_sid, 1);
                let overlap_bins =
                    (q.end.min(p.end) - q.start.max(p.start)) as f32 / BIN_SIZE as f32;
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

    fn split_db(db: &[Interval]) -> (Vec<Position>, Vec<u32>) {
        let positions = db.iter().map(|iv| Position { start: iv.start, end: iv.end }).collect();
        let sids = db.iter().map(|iv| iv.sid).collect();
        (positions, sids)
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
        let (positions, sids) = split_db(db);
        let mut results = DenseMatrix::new(nrows, ncols);
        let mut overlap = OverlapMatrix::new(nrows, ncols);
        query_sweep(&positions, &sids, layer_max_size, &jt, query, &mut results, &mut overlap);
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
        assert!((overlap.get(0, 0) - 0.25).abs() < 1e-5, "Q0∩D0={}", overlap.get(0, 0));
        assert!((overlap.get(1, 0) - 0.25).abs() < 1e-5, "Q1∩D0={}", overlap.get(1, 0));
        assert!((overlap.get(1, 1) - 0.25).abs() < 1e-5, "Q1∩D1={}", overlap.get(1, 1));
    }

    #[test]
    fn test_dead_zone_skipped_via_jump_table() {
        let db = vec![iv(0, 50, 0), iv(10_000, 10_050, 1)];
        let query = vec![iv(10_000, 10_100, 0)];
        let (results, overlap) = sweep(&db, 100, &query, 1, 2);

        assert_eq!(results.get(0, 0), 0);
        assert_eq!(results.get(0, 1), 1);
        assert_eq!(overlap.get(0, 0), 0.0);
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
