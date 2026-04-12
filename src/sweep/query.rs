use std::cmp::Reverse;
use std::collections::BinaryHeap;

use bitvec::prelude::*;

use crate::core::Interval;
use crate::io::MappedJumpTable;
use crate::matrix::DenseMatrix;

/// Single-pass dual-active-set sweep over a flat sorted layer.
///
/// Counts set-level intersections between a query block and a database layer
/// without allocating any intermediate result intervals.
///
/// ### Inputs
///
/// - `db_layer`: flat sorted `&[Interval]` from a `MappedLayer`. Every
///   interval is present exactly once, so no deduplication is needed.
/// - `layer_max_size`: the exclusive upper bound on interval sizes in this
///   layer (`LayerConfig::layer_max_size(K)`). Used for the fast-forward
///   heuristic.
/// - `jump_table`: O(1) fast-forward index for this layer. The dead-zone skip
///   uses a table lookup + short linear scan instead of binary search.
/// - `query_block`: sorted `&[Interval]` where `sid` is the row index in
///   `results` (Q_SID).
/// - `results`: `[Q_sids][D_sids]` dense accumulator, mutated in place.
///
/// ### Algorithm
///
/// At each step, the sweep advances to time `T = min(next D end, next Q end,
/// next D start, next Q start)` and processes exactly one event:
///
/// 1. **D ends** — pop from `active_d`; decrement `overlaps_frame[d_sid]`;
///    if it hits zero, clear `active_mask[d_sid]`.
/// 2. **Q ends** — pop from `active_q`.
/// 3. **D starts** — push to `active_d`; increment `overlaps_frame[d_sid]`;
///    if it was zero, set `active_mask[d_sid]`;
///    forward-match against every active Q: `results[q_sid][d_sid] += 1`.
/// 4. **Q starts** — push to `active_q`; sparsity-aware catch-up via
///    `active_mask.iter_ones()` (TZCNT): `results[q_sid][d_sid] += overlaps_frame[d_sid]`
///    for each set bit, skipping zero entries entirely.
///
/// ### Fast-forward heuristic
///
/// When `active_q` is empty and the next D interval cannot possibly overlap
/// the next Q interval (i.e. `d.start + layer_max_size < q.start`), the
/// entire D dead zone is skipped.  If a jump table is present, an O(1) tile
/// lookup places the cursor within one tile-width of the target; a short
/// forward scan (≤ tile_size coordinate units) completes the positioning.
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

    let num_d_sids = results.num_cols();
    let mut overlaps_frame = vec![0u32; num_d_sids];
    // Tracks which D_sids have a non-zero entry in overlaps_frame.
    // Enables sparsity-aware catch-up in Event 4 via TZCNT (iter_ones).
    let mut active_mask: BitVec<u64, Lsb0> = bitvec![u64, Lsb0; 0; num_d_sids];

    // Min-heap by end: Reverse((end, sid))
    let mut active_d: BinaryHeap<Reverse<(u32, u32)>> = BinaryHeap::new();
    let mut active_q: BinaryHeap<Reverse<(u32, u32)>> = BinaryHeap::new();

    let mut d_cursor = 0usize;
    let mut q_cursor = 0usize;

    loop {
        // ----------------------------------------------------------------
        // Fast-forward: skip dead zone when no queries are active
        // ----------------------------------------------------------------
        if active_q.is_empty() && q_cursor < query_block.len() && d_cursor < db_layer.len() {
            let next_q_start = query_block[q_cursor].start;
            // Any D interval starting before (next_q_start - layer_max_size) has
            // end <= start + layer_max_size < next_q_start, so it can never overlap.
            let dead_zone_end = next_q_start.saturating_sub(layer_max_size);

            if db_layer[d_cursor].start < dead_zone_end {
                // Clear any stale active_d state (all those intervals are expired).
                active_d.clear();
                overlaps_frame.fill(0);
                active_mask.fill(false);
                // O(1) Jump table lookup to get the guaranteed lower bound index
                let approx = jump_table.lookup(dead_zone_end).min(db_layer.len());

                // O(log D) Binary search on the remaining slice
                // We slice from `approx` to the end, find the relative offset, and add it back.
                let relative_offset =
                    db_layer[approx..].partition_point(|d| d.start < dead_zone_end);
                d_cursor = approx + relative_offset;

                if d_cursor >= db_layer.len() {
                    break; // No more DB intervals can overlap any query
                }
            }
        }

        let next_d_end = active_d.peek().map(|x| x.0.0).unwrap_or(u32::MAX);
        let next_q_end = active_q.peek().map(|x| x.0.0).unwrap_or(u32::MAX);
        let next_d_start = db_layer.get(d_cursor).map(|d| d.start).unwrap_or(u32::MAX);
        let next_q_start = query_block
            .get(q_cursor)
            .map(|q| q.start)
            .unwrap_or(u32::MAX);

        let cursor = next_d_end
            .min(next_q_end)
            .min(next_d_start)
            .min(next_q_start);

        if cursor == u32::MAX {
            break;
        }

        // ----------------------------------------------------------------
        // Event 1: D ends — remove expired DB interval, decrement frame
        // ----------------------------------------------------------------
        if cursor == next_d_end {
            let Reverse((_, d_sid)) = active_d.pop().unwrap();
            overlaps_frame[d_sid as usize] -= 1;
            if overlaps_frame[d_sid as usize] == 0 {
                active_mask.set(d_sid as usize, false);
            }
            continue;
        }

        // ----------------------------------------------------------------
        // Event 2: Q ends — remove expired query interval
        // ----------------------------------------------------------------
        if cursor == next_q_end {
            active_q.pop();
            continue;
        }

        // ----------------------------------------------------------------
        // Event 3: D starts — add to active set, forward-match active queries
        // ----------------------------------------------------------------
        if cursor == next_d_start {
            let d = db_layer[d_cursor];
            active_d.push(Reverse((d.end, d.sid)));
            overlaps_frame[d.sid as usize] += 1;
            if overlaps_frame[d.sid as usize] == 1 {
                active_mask.set(d.sid as usize, true);
            }

            // Forward match: every already-active query overlaps this new D.
            for Reverse((q_end, q_sid)) in active_q.iter() {
                if *q_end > d.start {
                    results.add(*q_sid as usize, d.sid as usize, 1);
                }
            }

            d_cursor += 1;
            continue;
        }

        // ----------------------------------------------------------------
        // Event 4: Q starts — catch up with current overlaps_frame (SIMD)
        // ----------------------------------------------------------------
        if cursor == next_q_start {
            let q = query_block[q_cursor];
            active_q.push(Reverse((q.end, q.sid)));

            // Sparsity-aware catch-up: iterate only non-zero D_sids via TZCNT.
            let row = results.row_mut(q.sid as usize);
            for d_sid in active_mask.iter_ones() {
                row[d_sid] += overlaps_frame[d_sid];
            }

            q_cursor += 1;
            continue;
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

    /// Build a real on-disk jump table from `db` with tile_size = layer_max_size / 2.
    /// Returns the tempfile (must stay alive) and the opened table.
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

        // Q0=[0,50) ∩ D0=[25,75) → overlap at [25,50) → count 1
        assert_eq!(results.get(0, 0), 1);
        // Q0=[0,50) ∩ D1=[75,125) → 50 ≤ 75, no overlap
        assert_eq!(results.get(0, 1), 0);
        // Q1=[50,100) ∩ D0=[25,75) → overlap at [50,75) → count 1
        assert_eq!(results.get(1, 0), 1);
        // Q1=[50,100) ∩ D1=[75,125) → overlap at [75,100) → count 1
        assert_eq!(results.get(1, 1), 1);
    }

    #[test]
    fn test_fast_forward_skips_dead_zone() {
        // DB has intervals at 0-50 and 10000-10050
        // Query is only at 10000-10100; fast-forward must skip the 0-50 interval.
        let db = vec![iv(0, 50, 0), iv(10_000, 10_050, 1)];
        let query = vec![iv(10_000, 10_100, 0)];
        let mut results = DenseMatrix::new(1, 2);
        let (_f, jt) = make_jt(&db, 100);
        query_sweep(&db, 100, &jt, &query, &mut results);

        assert_eq!(results.get(0, 0), 0); // D0 in dead zone
        assert_eq!(results.get(0, 1), 1); // D1 overlaps
    }

    #[test]
    fn test_fast_forward_same_result_as_naive() {
        // Large gap: fast-forward (layer_max=10) and naive (layer_max=u32::MAX) must agree.
        let db: Vec<Interval> = (0..100)
            .map(|i| iv(i * 10, i * 10 + 5, i % 3))
            .chain((0..10).map(|i| iv(100_000 + i * 10, 100_000 + i * 10 + 5, i % 3)))
            .collect();
        let query = vec![iv(100_000, 100_050, 0), iv(100_060, 100_100, 1)];

        let mut fast = DenseMatrix::new(2, 3);
        let mut naive = DenseMatrix::new(2, 3);

        let (_ff, jt_fast) = make_jt(&db, 10);
        let (_fn, jt_naive) = make_jt(&db, u32::MAX);
        query_sweep(&db, 10, &jt_fast, &query, &mut fast);
        query_sweep(&db, u32::MAX, &jt_naive, &query, &mut naive);

        for r in 0..2 {
            for c in 0..3 {
                assert_eq!(fast.get(r, c), naive.get(r, c), "mismatch at ({r},{c})");
            }
        }
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
