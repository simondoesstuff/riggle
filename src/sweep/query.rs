use crate::core::Interval;
use crate::io::MappedJumpTable;
use crate::matrix::DenseMatrix;

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
/// - `results`: `[Q_sids][D_sids]` dense accumulator, mutated in place.
pub fn query_sweep(
    db_layer: &[Interval],
    layer_max_size: u32,
    jump_table: &MappedJumpTable,
    query_block: &[Interval],
    results: &mut DenseMatrix,
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
        for d in &db_layer[start_idx..] {
            if d.start >= q.end {
                break;
            }
            if d.end > q.start {
                results.add(q.sid as usize, d.sid as usize, 1);
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

    #[test]
    fn test_basic_overlap() {
        let db = vec![iv(50, 150, 0)];
        let query = vec![iv(40, 160, 0)];
        let mut results = DenseMatrix::new(1, 1);
        let (_f, jt) = make_jt(&db, 200);
        query_sweep(&db, 200, &jt, &query, &mut results);
        assert_eq!(results.get(0, 0), 1);
    }

    #[test]
    fn test_no_overlap() {
        let db = vec![iv(50, 100, 0)];
        let query = vec![iv(200, 300, 0)];
        let mut results = DenseMatrix::new(1, 1);
        let (_f, jt) = make_jt(&db, 200);
        query_sweep(&db, 200, &jt, &query, &mut results);
        assert_eq!(results.get(0, 0), 0);
    }

    #[test]
    fn test_multiple_queries_and_db() {
        // D0=[25,75), D1=[75,125)
        let db = vec![iv(25, 75, 0), iv(75, 125, 1)];
        // Q0=[0,50), Q1=[50,100)
        let query = vec![iv(0, 50, 0), iv(50, 100, 1)];
        let mut results = DenseMatrix::new(2, 2);
        let (_f, jt) = make_jt(&db, 200);
        query_sweep(&db, 200, &jt, &query, &mut results);

        assert_eq!(results.get(0, 0), 1);
        assert_eq!(results.get(0, 1), 0);
        assert_eq!(results.get(1, 0), 1);
        assert_eq!(results.get(1, 1), 1);
    }

    #[test]
    fn test_dead_zone_skipped_via_jump_table() {
        // DB has intervals at 0-50 and 10000-10050; query only at 10000-10100.
        let db = vec![iv(0, 50, 0), iv(10_000, 10_050, 1)];
        let query = vec![iv(10_000, 10_100, 0)];
        let mut results = DenseMatrix::new(1, 2);
        let (_f, jt) = make_jt(&db, 100);
        query_sweep(&db, 100, &jt, &query, &mut results);

        assert_eq!(results.get(0, 0), 0); // D0 far away, no overlap
        assert_eq!(results.get(0, 1), 1); // D1 overlaps
    }

    #[test]
    fn test_touching_intervals_not_overlapping() {
        let db = vec![iv(0, 50, 0)];
        let query = vec![iv(50, 100, 0)];
        let mut results = DenseMatrix::new(1, 1);
        let (_f, jt) = make_jt(&db, 200);
        query_sweep(&db, 200, &jt, &query, &mut results);
        assert_eq!(results.get(0, 0), 0);
    }

    #[test]
    fn test_empty_db() {
        let db: Vec<Interval> = vec![];
        let query = vec![iv(0, 100, 0)];
        let mut results = DenseMatrix::new(1, 1);
        let (_f, jt) = make_jt(&db, 200);
        query_sweep(&db, 200, &jt, &query, &mut results);
        assert_eq!(results.get(0, 0), 0);
    }

    #[test]
    fn test_empty_query() {
        let db = vec![iv(0, 100, 0)];
        let query: Vec<Interval> = vec![];
        let mut results = DenseMatrix::new(1, 1);
        let (_f, jt) = make_jt(&db, 200);
        query_sweep(&db, 200, &jt, &query, &mut results);
        assert_eq!(results.get(0, 0), 0);
    }
}
