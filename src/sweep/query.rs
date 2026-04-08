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
/// The overlap between a query Q and database interval D begins exactly at 
/// `max(Q.start, D.start)`. Since tiles perfectly partition the coordinate space, 
/// this overlap start coordinate falls into exactly one tile. We simply only count 
/// the overlap if it "begins" in the current tile.
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
        |tile_idx, tile_start, tile_end, active_queries| {
            let tile = match mapped_chunk.get_tile(tile_idx) {
                Ok(t) => t,
                Err(_) => return, // Skip tiles we can't read
            };

            for (query_idx, query) in active_queries {
                process_intervals(tile, tile_start, tile_end, query, *query_idx, results, mask);
            }
        },
    );
}

/// Process intervals overlapping this tile
/// Uses binary search to find the first relevant interval, then iterates
fn process_intervals(
    tile: &ArchivedTile,
    tile_start: u32,
    tile_end: u32,
    query: &TaggedInterval,
    query_idx: usize,
    results: &mut DenseMatrix,
    mask: &mut BitwiseMask,
) {
    if tile.intervals.is_empty() {
        return;
    }

    // Iterate over intervals in this tile
    // Wait, the archived intervals are TaggedInterval.
    // They are sorted by start coordinate. We can binary search for intervals
    // that start before query.end. Actually, since we want intervals where
    // iv.end > query.start, we could just scan or binary search.
    // But since they are sorted by iv.start, we can binary search to find
    // the first interval where iv.end > query.start? No, end is not sorted.
    // However, they are sorted by start. So we can skip intervals where
    // iv.start >= query.end because they can't overlap.
    
    // Binary search is NOT easily possible for iv.end > query.start because
    // start and end are decoupled. However, since no interval spans an entire
    // tile, an interval overlapping query in this tile must have its start
    // or end in this tile.
    
    for entry in tile.intervals.iter() {
        let iv_start: u32 = entry.iv.start.into();
        let iv_end: u32 = entry.iv.end.into();
        let sid: u32 = entry.sid.into();

        // Stop if we've passed the query end, since intervals are sorted by start
        if iv_start >= query.iv.end {
            break;
        }

        // Check for overlap
        if iv_end > query.iv.start {
            // Deduplication Check
            let overlap_start = query.iv.start.max(iv_start);
            if overlap_start >= tile_start && overlap_start < tile_end {
                let sid_usize = sid as usize;
                if sid_usize < results.num_cols() {
                    results.add(query_idx, sid_usize, 1);
                    mask.flag(query_idx, sid_usize);
                }
            }
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
        tile0.intervals.push(TaggedInterval::new(50, 150, 0)); // interval 50-150, sid 0
        let mut tile1 = Tile::new(100);
        tile1.intervals.push(TaggedInterval::new(50, 150, 0)); // interval 50-150, sid 0

        let (file, tile_size) = create_test_chunk(vec![tile0, tile1], 100);
        let mapped = MappedChunk::open(file.path()).unwrap();

        let mut results = DenseMatrix::new(1, 1);
        let mut mask = BitwiseMask::new(1, 1);

        // Query that overlaps with the interval (50-150)
        let queries = vec![(0, TaggedInterval::new(40, 160, 100))];

        query_sweep(&mapped, tile_size, &queries, &mut results, &mut mask);

        // Should count the interval once
        assert_eq!(results.get(0, 0), 1);
    }

    #[test]
    fn test_query_sweep_crossing_tiles() {
        // Create a chunk with 3 tiles, tile_size = 100
        // Interval 50-250, sid 0. Overlaps tile 0, 1, 2
        let mut tile0 = Tile::new(0);
        tile0.intervals.push(TaggedInterval::new(50, 250, 0));
        let mut tile1 = Tile::new(100);
        tile1.intervals.push(TaggedInterval::new(50, 250, 0));
        let mut tile2 = Tile::new(200);
        tile2.intervals.push(TaggedInterval::new(50, 250, 0));


        let (file, tile_size) = create_test_chunk(vec![tile0, tile1, tile2], 100);
        let mapped = MappedChunk::open(file.path()).unwrap();

        let mut results = DenseMatrix::new(1, 1);
        let mut mask = BitwiseMask::new(1, 1);

        // Query that spans all tiles
        let queries = vec![(0, TaggedInterval::new(0, 300, 100))];

        query_sweep(&mapped, tile_size, &queries, &mut results, &mut mask);

        // Should count the interval once
        assert_eq!(results.get(0, 0), 1);
    }

    #[test]
    fn test_query_sweep_no_overlap() {
        let mut tile0 = Tile::new(0);
        tile0.intervals.push(TaggedInterval::new(50, 150, 0));

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
        tile0.intervals.push(TaggedInterval::new(25, 75, 0)); // interval 25-75, sid 0
        tile0.intervals.push(TaggedInterval::new(75, 125, 1)); // interval 75-125, sid 1

        let (file, tile_size) = create_test_chunk(vec![tile0], 100);
        let mapped = MappedChunk::open(file.path()).unwrap();

        let mut results = DenseMatrix::new(2, 2);
        let mut mask = BitwiseMask::new(2, 2);

        let queries = vec![
            (0, TaggedInterval::new(0, 50, 100)), // overlaps 25-75, sid 0
            (1, TaggedInterval::new(50, 100, 101)), // overlaps 75-125, sid 1
        ];

        query_sweep(&mapped, tile_size, &queries, &mut results, &mut mask);

        assert_eq!(results.get(0, 0), 1); // query 0 overlaps sid 0
        assert_eq!(results.get(0, 1), 0); // query 0 doesn't overlap sid 1
        assert_eq!(results.get(1, 0), 0); // query 1 doesn't overlap sid 0
        assert_eq!(results.get(1, 1), 1); // query 1 overlaps sid 1
    }
}
