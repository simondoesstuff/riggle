use voracious_radix_sort::RadixSort;

use crate::core::{OffsetSid, TaggedInterval, Tile};

use super::sweep_tiles;

/// Build tiles from a sorted batch of intervals within a chunk
///
/// # Arguments
/// * `chunk_bounds` - The interval covered by this chunk [start, end)
/// * `tile_size` - Size of each tile within the chunk
/// * `active_batch` - Intervals sorted by start coordinate that touch this chunk
///
/// # Returns
/// Vector of tiles covering the chunk, each containing interval information
///
/// # Note
/// Since tile_size is > the max interval size for each layer, no interval can
/// span an entire tile. Each interval either starts in the tile, ends in the tile,
/// or both - never spans it completely.
pub fn index_sweep(
    chunk_bounds: crate::core::Interval,
    tile_size: u32,
    active_batch: &[TaggedInterval],
) -> Vec<Tile> {
    if tile_size == 0 {
        return Vec::new();
    }

    let num_tiles = ((chunk_bounds.end - chunk_bounds.start) + tile_size - 1) / tile_size;
    let mut tiles: Vec<Tile> = (0..num_tiles)
        .map(|i| Tile::new(chunk_bounds.start + i * tile_size))
        .collect();

    if active_batch.is_empty() {
        return tiles;
    }

    // Use sweep_tiles to iterate over active intervals for each tile
    sweep_tiles(
        num_tiles as usize,
        chunk_bounds.start,
        tile_size,
        active_batch,
        |iv| iv.iv,
        |tile_idx, tile_start, tile_end, active_items| {
            let tile = &mut tiles[tile_idx];
            let tile_end = tile_end.min(chunk_bounds.end);

            for iv in active_items {
                // Classify the interval's relationship to this tile
                // Since tile_size > max_interval_size, no interval can span a tile
                let starts_in_tile = iv.iv.start >= tile_start && iv.iv.start < tile_end;
                let ends_in_tile = iv.iv.end > tile_start && iv.iv.end <= tile_end;

                if starts_in_tile {
                    let offset = iv.iv.start - tile_start;
                    tile.start_ivs.push(OffsetSid::new(offset, iv.sid));
                }
                if ends_in_tile {
                    let offset = iv.iv.end - tile_start;
                    tile.end_ivs.push(OffsetSid::new(offset, iv.sid));
                }
            }
        },
    );

    // Sort interval lists by offset for binary search during query
    // Uses radix sort O(n * w) where w is word size, vs O(n log n) for comparison sort
    for tile in &mut tiles {
        tile.start_ivs.voracious_sort();
        tile.end_ivs.voracious_sort();
    }

    tiles
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Interval;

    #[test]
    fn test_index_sweep_empty() {
        let bounds = Interval::new(0, 1000);
        let tiles = index_sweep(bounds, 100, &[]);
        assert_eq!(tiles.len(), 10);
        assert!(tiles.iter().all(|t| t.is_empty()));
    }

    #[test]
    fn test_index_sweep_single_interval() {
        let bounds = Interval::new(0, 1000);
        let intervals = vec![TaggedInterval::new(150, 250, 1)];
        let tiles = index_sweep(bounds, 100, &intervals);

        // Interval 150-250 should:
        // - Start in tile 1 (100-200) at offset 50
        // - End in tile 2 (200-300) at offset 50
        assert_eq!(tiles[1].start_ivs.len(), 1);
        assert_eq!(tiles[1].start_ivs[0], OffsetSid::new(50, 1)); // offset 50 from tile start

        assert_eq!(tiles[2].end_ivs.len(), 1);
        assert_eq!(tiles[2].end_ivs[0], OffsetSid::new(50, 1)); // offset 50 from tile start
    }

    #[test]
    fn test_index_sweep_crossing_interval() {
        let bounds = Interval::new(0, 500);
        // Interval crosses tiles but doesn't span any (with proper tile sizing)
        let intervals = vec![TaggedInterval::new(50, 150, 1)];
        let tiles = index_sweep(bounds, 100, &intervals);

        // Interval 50-150:
        // - Tile 0 (0-100): starts at offset 50
        // - Tile 1 (100-200): ends at offset 50

        assert_eq!(tiles[0].start_ivs.len(), 1); // starts at offset 50
        assert_eq!(tiles[0].start_ivs[0], OffsetSid::new(50, 1));

        assert_eq!(tiles[1].end_ivs.len(), 1); // ends at offset 50
        assert_eq!(tiles[1].end_ivs[0], OffsetSid::new(50, 1));
    }

    #[test]
    fn test_index_sweep_multiple_intervals_same_sid() {
        let bounds = Interval::new(0, 300);
        // Two intervals from same sid, both start and end in same tiles
        let intervals = vec![
            TaggedInterval::new(50, 150, 1),
            TaggedInterval::new(60, 160, 1),
        ];
        let tiles = index_sweep(bounds, 100, &intervals);

        // Tile 0: both intervals start
        assert_eq!(tiles[0].start_ivs.len(), 2);
        // Tile 1: both intervals end
        assert_eq!(tiles[1].end_ivs.len(), 2);
    }

    #[test]
    fn test_index_sweep_different_sids() {
        let bounds = Interval::new(0, 200);
        let intervals = vec![
            TaggedInterval::new(50, 150, 1),
            TaggedInterval::new(50, 150, 2),
        ];
        let tiles = index_sweep(bounds, 100, &intervals);

        // Tile 0: both intervals start
        assert_eq!(tiles[0].start_ivs.len(), 2);
        assert!(tiles[0].start_ivs.iter().any(|e| e.sid == 1));
        assert!(tiles[0].start_ivs.iter().any(|e| e.sid == 2));
    }

    #[test]
    fn test_index_sweep_boundary_interval() {
        let bounds = Interval::new(0, 200);
        // Interval exactly at tile boundary
        let intervals = vec![TaggedInterval::new(100, 200, 1)];
        let tiles = index_sweep(bounds, 100, &intervals);

        // Tile 0 (0-100): no overlap (end is 100, exclusive)
        assert!(tiles[0].is_empty());

        // Tile 1 (100-200): starts and ends in this tile
        assert_eq!(tiles[1].start_ivs.len(), 1);
        assert_eq!(tiles[1].end_ivs.len(), 1);
    }

    #[test]
    fn test_index_sweep_point_interval() {
        // While we don't allow empty intervals, a 1-bp interval is fine
        let bounds = Interval::new(0, 100);
        let intervals = vec![TaggedInterval::new(50, 51, 1)];
        let tiles = index_sweep(bounds, 100, &intervals);

        // Single tile, interval starts and ends in it
        assert_eq!(tiles[0].start_ivs.len(), 1);
        assert_eq!(tiles[0].start_ivs[0], OffsetSid::new(50, 1));
        assert_eq!(tiles[0].end_ivs.len(), 1);
        assert_eq!(tiles[0].end_ivs[0], OffsetSid::new(51, 1));
    }

    #[test]
    fn test_index_sweep_sorted_input() {
        let bounds = Interval::new(0, 500);
        let intervals = vec![
            TaggedInterval::new(50, 100, 1),
            TaggedInterval::new(150, 200, 2),
            TaggedInterval::new(250, 350, 3),
            TaggedInterval::new(400, 450, 4),
        ];
        let tiles = index_sweep(bounds, 100, &intervals);

        // Tile 0 (0-100): interval 1 starts and ends
        assert_eq!(tiles[0].start_ivs.len(), 1);
        assert_eq!(tiles[0].end_ivs.len(), 1);

        // Tile 1 (100-200): interval 2 starts and ends
        assert_eq!(tiles[1].start_ivs.len(), 1);
        assert_eq!(tiles[1].end_ivs.len(), 1);

        // Tile 2 (200-300): interval 3 starts
        assert_eq!(tiles[2].start_ivs.len(), 1);

        // Tile 3 (300-400): interval 3 ends
        assert_eq!(tiles[3].end_ivs.len(), 1);

        // Tile 4 (400-500): interval 4 starts and ends
        assert_eq!(tiles[4].start_ivs.len(), 1);
        assert_eq!(tiles[4].end_ivs.len(), 1);
    }
}
