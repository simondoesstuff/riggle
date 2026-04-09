use std::cmp::Reverse;
use std::collections::BinaryHeap;

use crate::core::TaggedInterval;
use crate::io::MappedChunk;
use crate::matrix::DenseMatrix;

/// Query sweep over a mapped chunk to count intersections.
///
/// Implements a single-pass dual-active-set sweep line:
///
/// - **Active_D** (min-heap by end): tracks active database intervals.
/// - **Active_Q** (min-heap by end): tracks active query intervals.
/// - **overlaps_frame**: 1-D array `[D_sids] u32` — current count of active DB
///   intervals per source ID at the sweep line.
///
/// At each event time T (the minimum of all upcoming starts/ends), four event
/// types are processed in strict order (ends before starts, D before Q at the
/// same coordinate):
///
/// 1. **D ends**: pop expired DB intervals, decrement `overlaps_frame[sid]`.
/// 2. **Q ends**: pop expired query intervals.
/// 3. **D starts**: push DB interval to Active_D, increment `overlaps_frame[sid]`,
///    then forward-match against all currently active queries
///    (`results[q_idx][sid] += 1`).
/// 4. **Q starts**: push query to Active_Q, then catch-up by adding the entire
///    `overlaps_frame` row into `results[q_idx]` (SIMD vector addition).
///
/// **Fast-forward heuristic**: when Active_Q is empty, D intervals that start
/// and end before the next query start are skipped entirely.  D intervals that
/// start before but end after the next query start are added directly to
/// `overlaps_frame` so the catch-up step in Event 4 accounts for them.
///
/// **Deduplication**: each DB interval is extracted from its canonical tile
/// (the first tile it appears in within this chunk), so every interval is seen
/// exactly once regardless of how many tiles it spans.
pub fn query_sweep(
    chunk: &MappedChunk,
    tile_size: u32,
    query_batch: &[(usize, TaggedInterval)], // (query_idx, interval), sorted by start
    results: &mut DenseMatrix,
) {
    // TODO: unclear why query_batch is (usize, TaggedInterval) when the query SID is in the TaggedInterval.
    // I am going to use the TaggedInterval SID instead of the outer tuple SID -- to be adjusted elsewhere
    // in the codebase later.
    let queries = query_batch.iter().map(|(_, iv)| iv).collect::<Vec<_>>();

    if queries.is_empty() || chunk.num_tiles() == 0 {
        return;
    }

    let num_tiles = chunk.num_tiles();
    let num_sources = results.num_cols();

    // Sweep-line state
    let mut overlaps_frame = vec![0u32; num_sources];
    // Min-heap by end: Reverse((end, start, sid))
    let mut active_d: BinaryHeap<Reverse<(u32, u32, u32)>> = BinaryHeap::new();
    // Min-heap by end: Reverse((end, sid))
    let mut active_q: BinaryHeap<Reverse<(u32, u32)>> = BinaryHeap::new();

    // Iteration tracking
    let mut active_tile_idx = 0usize;
    let mut d_idx = 0usize;

    for next_query in queries.iter() {
        // select tile
        // allows jumps due to large gaps between query intervals

        let min_tile_idx = (next_query.iv.start / tile_size) as usize;

        if active_tile_idx < min_tile_idx {
            // move tile, wrap database cursor
            d_idx = 0;
            active_tile_idx = min_tile_idx;
            continue;
        }

        // TODO: confirm this repeated get_tile isn't incurring
        // unnecessary overhead
        let active_tile = &chunk.get_tile(active_tile_idx).unwrap().intervals;

        if d_idx >= active_tile.len() {
            active_tile_idx = (active_tile_idx + 1).min(num_tiles - 1);
            continue;
        }

        // advance the cursor

        let next_d_end = active_d.peek().map(|x| x.0.0).unwrap_or(u32::MAX);
        let next_q_end = active_q.peek().map(|x| x.0.0).unwrap_or(u32::MAX);
        let next_d_start = active_tile
            .get(d_idx)
            .map(|x| x.iv.start.to_native())
            .unwrap_or(u32::MAX);
        let next_q_start = next_query.iv.start;

        let cursor = next_d_end
            .min(next_q_end)
            .min(next_d_start)
            .min(next_q_start);

        if cursor == u32::MAX {
            break; // we've exhausted a stream
        }

        // event-based state machine

        // D ends
        // clean up the frame
        if cursor == next_d_end {
            let Reverse((_d_end, _d_start, d_sid)) = active_d.pop().unwrap();
            overlaps_frame[d_sid as usize] -= 1;
            continue;
        }

        // Q ends
        if cursor == next_q_end {
            active_q.pop().unwrap();
            continue;
        }

        // D starts
        // update active query SIDs
        if cursor == next_d_start {
            let d_start = active_tile[d_idx].iv.start.to_native();
            let d_end = active_tile[d_idx].iv.end.to_native();
            let d_sid = active_tile[d_idx].sid.to_native();

            active_d.push(Reverse((d_end, d_start, d_sid)));
            overlaps_frame[d_sid as usize] += 1;

            // forward-match against active queries
            for Reverse((q_end, q_sid)) in active_q.iter() {
                if *q_end > d_start {
                    results.add(*q_sid as usize, d_sid as usize, 1);
                }
            }

            d_idx += 1;
            continue;
        }

        // Q starts
        // update active database SIDs
        if cursor == next_q_start {
            let q_end = next_query.iv.end;
            let q_sid = next_query.sid;

            active_q.push(Reverse((q_end, q_sid)));

            // catch-up with the frame
            let row = results.row_mut(q_sid as usize);

            // The compiler should auto-vectorize this loop into AVX/SSE block additions
            for (d_sid, &count) in overlaps_frame.iter().enumerate() {
                row[d_sid] += count;
            }

            continue;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Tile;
    use crate::io::{ChunkHeader, write_chunk};
    use tempfile::NamedTempFile;

    fn create_test_chunk(tiles: Vec<Tile>, tile_size: u32) -> (NamedTempFile, u32) {
        let mut header = ChunkHeader::new(0, 0, 0, tiles.len() as u32 * tile_size);

        let tile_bytes: Vec<_> = tiles
            .iter()
            .map(|t| rkyv::to_bytes::<rkyv::rancor::Error>(t).unwrap())
            .collect();

        let mut offset = 0u32;
        for bytes in &tile_bytes {
            header.tile_offsets.push(offset);
            offset += bytes.len() as u32;
        }

        let file = NamedTempFile::new().unwrap();
        write_chunk(file.path(), &header, &tiles).unwrap();
        (file, tile_size)
    }

    #[test]
    fn test_query_sweep_basic() {
        // Interval 50-150 appears in both tiles; should be counted exactly once.
        let mut tile0 = Tile::new(0);
        tile0.intervals.push(TaggedInterval::new(50, 150, 0));
        let mut tile1 = Tile::new(100);
        tile1.intervals.push(TaggedInterval::new(50, 150, 0));

        let (file, tile_size) = create_test_chunk(vec![tile0, tile1], 100);
        let mapped = MappedChunk::open(file.path()).unwrap();

        let mut results = DenseMatrix::new(1, 1);
        let queries = vec![(0, TaggedInterval::new(40, 160, 100))];

        query_sweep(&mapped, tile_size, &queries, &mut results);

        assert_eq!(results.get(0, 0), 1);
    }

    #[test]
    fn test_query_sweep_crossing_tiles() {
        // Interval 50-250 spans three tiles; deduplication must yield count of 1.
        let mut tile0 = Tile::new(0);
        tile0.intervals.push(TaggedInterval::new(50, 250, 0));
        let mut tile1 = Tile::new(100);
        tile1.intervals.push(TaggedInterval::new(50, 250, 0));
        let mut tile2 = Tile::new(200);
        tile2.intervals.push(TaggedInterval::new(50, 250, 0));

        let (file, tile_size) = create_test_chunk(vec![tile0, tile1, tile2], 100);
        let mapped = MappedChunk::open(file.path()).unwrap();

        let mut results = DenseMatrix::new(1, 1);
        let queries = vec![(0, TaggedInterval::new(0, 300, 100))];

        query_sweep(&mapped, tile_size, &queries, &mut results);

        assert_eq!(results.get(0, 0), 1);
    }

    #[test]
    fn test_query_sweep_no_overlap() {
        let mut tile0 = Tile::new(0);
        tile0.intervals.push(TaggedInterval::new(50, 150, 0));

        let (file, tile_size) = create_test_chunk(vec![tile0], 100);
        let mapped = MappedChunk::open(file.path()).unwrap();

        let mut results = DenseMatrix::new(1, 1);
        let queries = vec![(0, TaggedInterval::new(200, 300, 100))];

        query_sweep(&mapped, tile_size, &queries, &mut results);

        assert_eq!(results.get(0, 0), 0);
    }

    #[test]
    fn test_query_sweep_multiple_queries() {
        let mut tile0 = Tile::new(0);
        tile0.intervals.push(TaggedInterval::new(25, 75, 0)); // sid 0
        tile0.intervals.push(TaggedInterval::new(75, 125, 1)); // sid 1

        let (file, tile_size) = create_test_chunk(vec![tile0], 100);
        let mapped = MappedChunk::open(file.path()).unwrap();

        let mut results = DenseMatrix::new(2, 2);

        let queries = vec![
            (0, TaggedInterval::new(0, 50, 100)),   // Q0 = [0, 50)
            (1, TaggedInterval::new(50, 100, 101)), // Q1 = [50, 100)
        ];

        query_sweep(&mapped, tile_size, &queries, &mut results);

        // Q0=[0,50) ∩ D0=[25,75) = [25,50) → count 1
        assert_eq!(results.get(0, 0), 1);
        // Q0=[0,50) ∩ D1=[75,125): 50 ≤ 75, no overlap
        assert_eq!(results.get(0, 1), 0);
        // Q1=[50,100) ∩ D0=[25,75) = [50,75) → count 1
        assert_eq!(results.get(1, 0), 1);
        // Q1=[50,100) ∩ D1=[75,125) = [75,100) → count 1
        assert_eq!(results.get(1, 1), 1);
    }
}
