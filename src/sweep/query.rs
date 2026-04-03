use crate::core::{ArchivedTile, TaggedInterval};
use crate::io::MappedChunk;
use crate::matrix::{BitwiseMask, DenseMatrix};

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

    // For each query, track which sids we've already counted from running_counts
    // to avoid double-counting when a query spans multiple tiles
    let mut seen_running: Vec<Vec<u32>> = vec![Vec::new(); query_batch.len()];

    // Head pointer into query_batch
    let mut head = 0;

    // Process each tile
    for tile_id in 0..mapped_chunk.num_tiles() {
        let tile = match mapped_chunk.get_tile(tile_id) {
            Ok(t) => t,
            Err(_) => continue,
        };

        let tile_start = chunk_start + tile_id as u32 * tile_size;
        let tile_end = tile_start + tile_size;

        // Advance head past queries that end before this tile
        while head < query_batch.len() && query_batch[head].1.iv.end <= tile_start {
            head += 1;
        }

        // Process all queries that intersect this tile
        let mut scan = head;
        while scan < query_batch.len() && query_batch[scan].1.iv.start < tile_end {
            let (query_idx, query) = &query_batch[scan];

            // Process running counts (intervals spanning entire tile)
            process_running_counts(
                tile,
                query,
                *query_idx,
                results,
                mask,
                &mut seen_running[scan - head],
            );

            // Process start_ivs (intervals starting in this tile)
            process_start_ivs(tile, tile_start, query, *query_idx, results, mask);

            // Process end_ivs (intervals ending in this tile)
            process_end_ivs(
                tile,
                tile_start,
                query,
                *query_idx,
                results,
                mask,
                &seen_running[scan - head],
            );

            scan += 1;
        }

        // Clear seen_running for queries that no longer intersect future tiles
        for (i, (_, q)) in query_batch.iter().enumerate().skip(head).take(scan - head) {
            if q.iv.end <= tile_end {
                // Query ends in this tile, clear its tracking
                if i >= head && i - head < seen_running.len() {
                    seen_running[i - head].clear();
                }
            }
        }
    }
}

/// Process running counts from a tile for a query
fn process_running_counts(
    tile: &ArchivedTile,
    _query: &TaggedInterval,
    query_idx: usize,
    results: &mut DenseMatrix,
    mask: &mut BitwiseMask,
    seen: &mut Vec<u32>,
) {
    for entry in tile.running_counts.iter() {
        let sid: u32 = entry.0.into();
        let count: u32 = entry.1.into();

        // Check if we've already counted this sid for this query
        // (from a previous tile it also spans)
        if seen.contains(&sid) {
            continue;
        }

        // Mark as seen
        seen.push(sid);

        // Add count
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
    let query_start_offset = query.iv.start.saturating_sub(tile_start) as u16;

    // Find first interval with offset >= query_start_offset
    let start_idx = tile.start_ivs.partition_point(|entry| {
        let offset: u16 = entry.0.into();
        offset < query_start_offset
    });

    // Iterate from start_idx until we exceed query.end
    for entry in tile.start_ivs.iter().skip(start_idx) {
        let offset: u16 = entry.0.into();
        let sid: u32 = entry.1.into();
        let iv_start = tile_start + offset as u32;

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
fn process_end_ivs(
    tile: &ArchivedTile,
    tile_start: u32,
    query: &TaggedInterval,
    query_idx: usize,
    results: &mut DenseMatrix,
    mask: &mut BitwiseMask,
    seen: &[u32],
) {
    if tile.end_ivs.is_empty() {
        return;
    }

    // Binary search to find first interval that might overlap with query
    // We want intervals where: iv_end > query.start (interval ends after query starts)
    // Since iv_end = tile_start + offset, we need offset > query.start - tile_start
    let query_start_offset = query.iv.start.saturating_sub(tile_start) as u16;

    // Find first interval with offset > query_start_offset
    let start_idx = tile.end_ivs.partition_point(|entry| {
        let offset: u16 = entry.0.into();
        offset <= query_start_offset
    });

    // Iterate from start_idx
    for entry in tile.end_ivs.iter().skip(start_idx) {
        let sid: u32 = entry.1.into();

        // The interval ends in this tile, started somewhere before
        // We overlap if query.start < iv_end (which we ensured with binary search)
        // Note: We don't need to check upper bound since end_ivs only contains
        // intervals ending IN this tile (before tile_end)

        // Don't double-count if we already counted this from running_counts
        if seen.contains(&sid) {
            continue;
        }

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
    use crate::core::Tile;
    use crate::io::{write_chunk, ChunkHeader};
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
        tile0.start_ivs.push((50, 0)); // interval starting at 50, sid 0

        let mut tile1 = Tile::new(100);
        tile1.end_ivs.push((50, 0)); // interval ending at 150, sid 0

        let (file, tile_size) = create_test_chunk(vec![tile0, tile1], 100);
        let mapped = MappedChunk::open(file.path()).unwrap();

        let mut results = DenseMatrix::new(1, 1);
        let mut mask = BitwiseMask::new(1, 1);

        // Query that overlaps with the interval (50-150)
        let queries = vec![(0, TaggedInterval::new(40, 160, 100))];

        query_sweep(&mapped, tile_size, &queries, &mut results, &mut mask);

        // Should count the interval once
        assert_eq!(results.get(0, 0), 2); // start_iv + end_iv
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
        tile0.start_ivs.push((50, 0)); // interval at 50

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
        tile0.start_ivs.push((25, 0)); // interval at 25, sid 0
        tile0.start_ivs.push((75, 1)); // interval at 75, sid 1

        let (file, tile_size) = create_test_chunk(vec![tile0], 100);
        let mapped = MappedChunk::open(file.path()).unwrap();

        let mut results = DenseMatrix::new(2, 2);
        let mut mask = BitwiseMask::new(2, 2);

        let queries = vec![
            (0, TaggedInterval::new(0, 50, 100)),  // overlaps interval at 25
            (1, TaggedInterval::new(50, 100, 101)), // overlaps interval at 75
        ];

        query_sweep(&mapped, tile_size, &queries, &mut results, &mut mask);

        assert_eq!(results.get(0, 0), 1); // query 0 overlaps sid 0
        assert_eq!(results.get(0, 1), 0); // query 0 doesn't overlap sid 1
        assert_eq!(results.get(1, 0), 0); // query 1 doesn't overlap sid 0
        assert_eq!(results.get(1, 1), 1); // query 1 overlaps sid 1
    }
}
