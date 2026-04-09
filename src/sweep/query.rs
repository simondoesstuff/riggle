use std::cmp::Reverse;
use std::collections::BinaryHeap;

use crate::core::Interval;
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
/// - `query_block`: sorted `&[Interval]` where `sid` is the row index in
///   `results` (Q_SID).
/// - `results`: `[Q_sids][D_sids]` dense accumulator, mutated in place.
///
/// ### Algorithm
///
/// At each step, the sweep advances to time `T = min(next D end, next Q end,
/// next D start, next Q start)` and processes exactly one event:
///
/// 1. **D ends** — pop from `active_d`; decrement `overlaps_frame[d_sid]`.
/// 2. **Q ends** — pop from `active_q`.
/// 3. **D starts** — push to `active_d`; increment `overlaps_frame[d_sid]`;
///    forward-match against every active Q: `results[q_sid][d_sid] += 1`.
/// 4. **Q starts** — push to `active_q`; catch-up via SIMD-friendly loop:
///    `results[q_sid] += overlaps_frame`.
///
/// ### Fast-forward heuristic
///
/// When `active_q` is empty and the next D interval cannot possibly overlap
/// the next Q interval (i.e. `d.start + layer_max_size < q.start`), the
/// entire D dead zone is skipped via binary search.
pub fn query_sweep(
    db_layer: &[Interval],
    layer_max_size: u32,
    query_block: &[Interval],
    results: &mut DenseMatrix,
) {
    if query_block.is_empty() || db_layer.is_empty() {
        return;
    }

    let num_d_sids = results.num_cols();
    let mut overlaps_frame = vec![0u32; num_d_sids];

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
                // Binary search for the first D whose start is within reach.
                d_cursor = db_layer.partition_point(|d| d.start < dead_zone_end);
                if d_cursor >= db_layer.len() {
                    break; // No more DB intervals can overlap any query
                }
            }
        }

        let next_d_end = active_d.peek().map(|x| x.0.0).unwrap_or(u32::MAX);
        let next_q_end = active_q.peek().map(|x| x.0.0).unwrap_or(u32::MAX);
        let next_d_start = db_layer.get(d_cursor).map(|d| d.start).unwrap_or(u32::MAX);
        let next_q_start = query_block.get(q_cursor).map(|q| q.start).unwrap_or(u32::MAX);

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

            // SIMD-friendly: compiler auto-vectorizes this loop (AVX/SSE).
            let row = results.row_mut(q.sid as usize);
            for (d_sid, &count) in overlaps_frame.iter().enumerate() {
                row[d_sid] += count;
            }

            q_cursor += 1;
            continue;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(start: u32, end: u32, sid: u32) -> Interval {
        Interval::new(start, end, sid)
    }

    #[test]
    fn test_basic_overlap() {
        let db = vec![iv(50, 150, 0)];
        let query = vec![iv(40, 160, 0)];
        let mut results = DenseMatrix::new(1, 1);
        query_sweep(&db, 200, &query, &mut results);
        assert_eq!(results.get(0, 0), 1);
    }

    #[test]
    fn test_no_overlap() {
        let db = vec![iv(50, 100, 0)];
        let query = vec![iv(200, 300, 0)];
        let mut results = DenseMatrix::new(1, 1);
        query_sweep(&db, 200, &query, &mut results);
        assert_eq!(results.get(0, 0), 0);
    }

    #[test]
    fn test_multiple_queries_and_db() {
        // D0=[25,75), D1=[75,125)
        let db = vec![iv(25, 75, 0), iv(75, 125, 1)];
        // Q0=[0,50), Q1=[50,100)
        let query = vec![iv(0, 50, 0), iv(50, 100, 1)];
        let mut results = DenseMatrix::new(2, 2);
        query_sweep(&db, 200, &query, &mut results);

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
        // Query is only at 10000-10100
        // Fast-forward must skip the 0-50 interval
        let db = vec![iv(0, 50, 0), iv(10_000, 10_050, 1)];
        let query = vec![iv(10_000, 10_100, 0)];
        let mut results = DenseMatrix::new(1, 2);
        // layer_max_size = 100 (intervals <= 100 in this layer)
        query_sweep(&db, 100, &query, &mut results);

        // D0=[0,50) is in dead zone for Q=[10000,...) with layer_max=100
        assert_eq!(results.get(0, 0), 0); // no overlap with D0
        assert_eq!(results.get(0, 1), 1); // overlap with D1
    }

    #[test]
    fn test_fast_forward_same_result_as_naive() {
        // Construct a scenario with a large gap; the fast-forward path and the
        // naive path (with a tiny layer_max_size that never triggers) must agree.
        let db: Vec<Interval> = (0..100)
            .map(|i| iv(i * 10, i * 10 + 5, i % 3))
            .chain((0..10).map(|i| iv(100_000 + i * 10, 100_000 + i * 10 + 5, i % 3)))
            .collect();
        let query = vec![iv(100_000, 100_050, 0), iv(100_060, 100_100, 1)];

        let mut fast = DenseMatrix::new(2, 3);
        let mut naive = DenseMatrix::new(2, 3);

        query_sweep(&db, 10, &query, &mut fast);     // triggers fast-forward
        query_sweep(&db, u32::MAX, &query, &mut naive); // no fast-forward

        for r in 0..2 {
            for c in 0..3 {
                assert_eq!(
                    fast.get(r, c),
                    naive.get(r, c),
                    "mismatch at ({r},{c})"
                );
            }
        }
    }

    #[test]
    fn test_touching_intervals_not_overlapping() {
        // Intervals that touch at a boundary: [0,50) and [50,100) do NOT overlap.
        let db = vec![iv(0, 50, 0)];
        let query = vec![iv(50, 100, 0)];
        let mut results = DenseMatrix::new(1, 1);
        query_sweep(&db, 200, &query, &mut results);
        assert_eq!(results.get(0, 0), 0);
    }

    #[test]
    fn test_empty_db() {
        let db: Vec<Interval> = vec![];
        let query = vec![iv(0, 100, 0)];
        let mut results = DenseMatrix::new(1, 1);
        query_sweep(&db, 200, &query, &mut results);
        assert_eq!(results.get(0, 0), 0);
    }

    #[test]
    fn test_empty_query() {
        let db = vec![iv(0, 100, 0)];
        let query: Vec<Interval> = vec![];
        let mut results = DenseMatrix::new(1, 1);
        query_sweep(&db, 200, &query, &mut results);
        assert_eq!(results.get(0, 0), 0);
    }
}
