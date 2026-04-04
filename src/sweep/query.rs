use crate::core::{ArchivedTile, TaggedInterval};
use crate::io::MappedChunk;
use crate::matrix::{BitwiseMask, DenseMatrix};

use super::sweep_tiles;

/// Query sweep over a mapped chunk to count intersections
///
/// # Arguments
/// * `mapped_chunk` - Memory-mapped chunk data
/// * `tile_size` - Size of each tile
/// * `query_batch` - Queries sorted by start coordinate
/// * `results` - Dense matrix to accumulate counts (rows = queries, cols = db sids)
/// * `mask` - Bitmask tracking non-zero entries
///
/// # Note
/// Query indices in `query_batch` are their positions in the sorted order,
/// which correspond to row indices in `results`.
///
/// # Deduplication Strategy
/// Intervals spanning multiple tiles appear in running_counts of each tile they span.
/// To avoid double-counting, we only count running_counts and end_ivs in the
/// "first overlap tile" for each query - determined by `query.start >= tile_start`.
/// start_ivs are always unique to their tile, so no deduplication needed.
pub fn query_sweep(
    mapped_chunk: &MappedChunk,
    tile_size: u32,
    query_batch: &[(usize, TaggedInterval)], // (query_index, interval)
    results: &mut DenseMatrix,
    mask: &mut BitwiseMask,
) {
    if query_batch.is_empty() || mapped_chunk.num_tiles() == 0 {
        return;
    }

    let chunk_start = mapped_chunk.start_coord();

    // Use sweep_tiles to iterate over active queries for each tile
    sweep_tiles(
        mapped_chunk.num_tiles(),
        chunk_start,
        tile_size,
        query_batch,
        |q| q.1.iv,
        |tile_idx, tile_start, _tile_end, active_queries| {
            let tile = match mapped_chunk.get_tile(tile_idx) {
                Ok(t) => t,
                Err(_) => return, // Skip tiles we can't read
            };

            for (query_idx, query) in active_queries {
                // Determine if this is the first tile where query overlaps
                // This is true when query.start falls within [tile_start, tile_end)
                // Used to avoid double-counting running_counts and end_ivs
                let is_first_overlap_tile = query.iv.start >= tile_start;

                // Process running counts only in first overlap tile
                if is_first_overlap_tile {
                    process_running_counts(tile, *query_idx, results, mask);
                }

                // Process start_ivs (intervals starting in this tile) - always unique
                process_start_ivs(tile, tile_start, query, *query_idx, results, mask);

                // Process end_ivs only in first overlap tile
                if is_first_overlap_tile {
                    process_end_ivs(tile, tile_start, query, *query_idx, results, mask);
                }
            }
        },
    );
}

/// Process running counts from a tile for a query
/// Called only for the first overlap tile (deduplication handled by caller)
fn process_running_counts(
    tile: &ArchivedTile,
    query_idx: usize,
    results: &mut DenseMatrix,
    mask: &mut BitwiseMask,
) {
    for entry in tile.running_counts.iter() {
        let sid: u32 = entry.0.into();
        let count: u32 = entry.1.into();

        let sid_usize = sid as usize;
        if sid_usize < results.num_cols() {
            results.add(query_idx, sid_usize, count);
            mask.flag(query_idx, sid_usize);
        }
    }
}

/// Process intervals starting in this tile
/// Uses binary search to find the first relevant interval, then iterates
fn process_start_ivs(
    tile: &ArchivedTile,
    tile_start: u32,
    query: &TaggedInterval,
    query_idx: usize,
    results: &mut DenseMatrix,
    mask: &mut BitwiseMask,
) {
    if tile.start_ivs.is_empty() {
        return;
    }

    // Binary search to find first interval that might overlap with query
    // We want intervals where: iv_start >= query.start AND iv_start < query.end
    // Since iv_start = tile_start + offset, we need offset >= query.start - tile_start
    let query_start_offset = query.iv.start.saturating_sub(tile_start);

    // Find first interval with offset >= query_start_offset
    let start_idx = tile.start_ivs.partition_point(|entry| {
        let offset: u32 = entry.offset.into();
        offset < query_start_offset
    });

    // Iterate from start_idx until we exceed query.end
    for entry in tile.start_ivs.iter().skip(start_idx) {
        let offset: u32 = entry.offset.into();
        let sid: u32 = entry.sid.into();
        let iv_start = tile_start + offset;

        // Stop if we've passed the query end
        if iv_start >= query.iv.end {
            break;
        }

        // Count this overlap
        let sid_usize = sid as usize;
        if sid_usize < results.num_cols() {
            results.add(query_idx, sid_usize, 1);
            mask.flag(query_idx, sid_usize);
        }
    }
}

/// Process intervals ending in this tile
/// Uses binary search to find the first relevant interval, then iterates
/// Called only for the first overlap tile (deduplication handled by caller)
fn process_end_ivs(
    tile: &ArchivedTile,
    tile_start: u32,
    query: &TaggedInterval,
    query_idx: usize,
    results: &mut DenseMatrix,
    mask: &mut BitwiseMask,
) {
    if tile.end_ivs.is_empty() {
        return;
    }

    // Binary search to find first interval that might overlap with query
    // We want intervals where: iv_end > query.start (interval ends after query starts)
    // Since iv_end = tile_start + offset, we need offset > query.start - tile_start
    let query_start_offset = query.iv.start.saturating_sub(tile_start);

    // Find first interval with offset > query_start_offset
    let start_idx = tile.end_ivs.partition_point(|entry| {
        let offset: u32 = entry.offset.into();
        offset <= query_start_offset
    });

    // Iterate from start_idx
    for entry in tile.end_ivs.iter().skip(start_idx) {
        let sid: u32 = entry.sid.into();

        // The interval ends in this tile, started somewhere before
        // We overlap if query.start < iv_end (which we ensured with binary search)
        // Note: We don't need to check upper bound since end_ivs only contains
        // intervals ending IN this tile (before tile_end)

        // Count this overlap
        let sid_usize = sid as usize;
        if sid_usize < results.num_cols() {
            results.add(query_idx, sid_usize, 1);
            mask.flag(query_idx, sid_usize);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{OffsetSid, Tile};
    use crate::io::{ChunkHeader, write_chunk};
    use tempfile::NamedTempFile;

    fn create_test_chunk(tiles: Vec<Tile>, tile_size: u32) -> (NamedTempFile, u32) {
        let mut header = ChunkHeader::new(0, 0, 0, tiles.len() as u32 * tile_size);

        // Calculate tile offsets
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
        // Create a chunk with 2 tiles, tile_size = 100
        let mut tile0 = Tile::new(0);
        tile0.start_ivs.push(OffsetSid::new(50, 0)); // interval starting at 50, sid 0

        let mut tile1 = Tile::new(100);
        tile1.end_ivs.push(OffsetSid::new(50, 0)); // interval ending at 150, sid 0

        let (file, tile_size) = create_test_chunk(vec![tile0, tile1], 100);
        let mapped = MappedChunk::open(file.path()).unwrap();

        let mut results = DenseMatrix::new(1, 1);
        let mut mask = BitwiseMask::new(1, 1);

        // Query that overlaps with the interval (50-150)
        let queries = vec![(0, TaggedInterval::new(40, 160, 100))];

        query_sweep(&mapped, tile_size, &queries, &mut results, &mut mask);

        // Should count the interval once (counted from start_ivs in tile0, not from end_ivs in tile1)
        assert_eq!(results.get(0, 0), 1);
    }

    #[test]
    fn test_query_sweep_running_counts() {
        // Create a chunk with 3 tiles, tile_size = 100
        let mut tile0 = Tile::new(0);
        let mut tile1 = Tile::new(100);
        let mut tile2 = Tile::new(200);

        // Interval spans all tiles
        tile0.running_counts.push((0, 1));
        tile1.running_counts.push((0, 1));
        tile2.running_counts.push((0, 1));

        let (file, tile_size) = create_test_chunk(vec![tile0, tile1, tile2], 100);
        let mapped = MappedChunk::open(file.path()).unwrap();

        let mut results = DenseMatrix::new(1, 1);
        let mut mask = BitwiseMask::new(1, 1);

        // Query that spans all tiles
        let queries = vec![(0, TaggedInterval::new(0, 300, 100))];

        query_sweep(&mapped, tile_size, &queries, &mut results, &mut mask);

        // Should count only once despite appearing in 3 tiles
        assert_eq!(results.get(0, 0), 1);
    }

    #[test]
    fn test_query_sweep_no_overlap() {
        let mut tile0 = Tile::new(0);
        tile0.start_ivs.push(OffsetSid::new(50, 0)); // interval at 50

        let (file, tile_size) = create_test_chunk(vec![tile0], 100);
        let mapped = MappedChunk::open(file.path()).unwrap();

        let mut results = DenseMatrix::new(1, 1);
        let mut mask = BitwiseMask::new(1, 1);

        // Query that doesn't overlap
        let queries = vec![(0, TaggedInterval::new(200, 300, 100))];

        query_sweep(&mapped, tile_size, &queries, &mut results, &mut mask);

        assert_eq!(results.get(0, 0), 0);
    }

    #[test]
    fn test_query_sweep_multiple_queries() {
        let mut tile0 = Tile::new(0);
        tile0.start_ivs.push(OffsetSid::new(25, 0)); // interval at 25, sid 0
        tile0.start_ivs.push(OffsetSid::new(75, 1)); // interval at 75, sid 1

        let (file, tile_size) = create_test_chunk(vec![tile0], 100);
        let mapped = MappedChunk::open(file.path()).unwrap();

        let mut results = DenseMatrix::new(2, 2);
        let mut mask = BitwiseMask::new(2, 2);

        let queries = vec![
            (0, TaggedInterval::new(0, 50, 100)), // overlaps interval at 25
            (1, TaggedInterval::new(50, 100, 101)), // overlaps interval at 75
        ];

        query_sweep(&mapped, tile_size, &queries, &mut results, &mut mask);

        assert_eq!(results.get(0, 0), 1); // query 0 overlaps sid 0
        assert_eq!(results.get(0, 1), 0); // query 0 doesn't overlap sid 1
        assert_eq!(results.get(1, 0), 0); // query 1 doesn't overlap sid 0
        assert_eq!(results.get(1, 1), 1); // query 1 overlaps sid 1
    }
}
