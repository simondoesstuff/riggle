use crate::core::{Interval, TaggedInterval, Tile};

/// Build tiles from a sorted batch of intervals within a chunk
///
/// # Arguments
/// * `chunk_bounds` - The interval covered by this chunk [start, end)
/// * `tile_size` - Size of each tile within the chunk
/// * `active_batch` - Intervals sorted by start coordinate that touch this chunk
///
/// # Returns
/// Vector of tiles covering the chunk, each containing interval information
pub fn index_sweep(
    chunk_bounds: Interval,
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

    // Head pointer into active_batch
    let mut head = 0;

    // Process each tile
    for (tile_idx, tile) in tiles.iter_mut().enumerate() {
        let tile_start = chunk_bounds.start + tile_idx as u32 * tile_size;
        let tile_end = (tile_start + tile_size).min(chunk_bounds.end);

        // Advance head to first interval that could intersect this tile
        // An interval intersects if interval.end > tile_start
        while head < active_batch.len() && active_batch[head].iv.end <= tile_start {
            head += 1;
        }

        // Process all intervals that intersect this tile
        let mut scan = head;
        while scan < active_batch.len() && active_batch[scan].iv.start < tile_end {
            let iv = &active_batch[scan];

            // Classify the interval's relationship to this tile
            let starts_before = iv.iv.start <= tile_start;
            let ends_after = iv.iv.end >= tile_end;
            let starts_in_tile = iv.iv.start >= tile_start && iv.iv.start < tile_end;
            let ends_in_tile = iv.iv.end > tile_start && iv.iv.end <= tile_end;

            if starts_before && ends_after {
                // Interval spans the entire tile - add to running counts
                // Aggregate by sid
                if let Some(entry) = tile.running_counts.iter_mut().find(|(s, _)| *s == iv.sid) {
                    entry.1 += 1;
                } else {
                    tile.running_counts.push((iv.sid, 1));
                }
            } else {
                // Interval partially overlaps the tile
                if starts_in_tile {
                    let offset = (iv.iv.start - tile_start) as u16;
                    tile.start_ivs.push((offset, iv.sid));
                }
                if ends_in_tile {
                    let offset = (iv.iv.end - tile_start) as u16;
                    tile.end_ivs.push((offset, iv.sid));
                }
            }

            scan += 1;
        }
    }

    // Sort interval lists by offset for binary search during query
    for tile in &mut tiles {
        tile.start_ivs.sort_by_key(|(offset, _)| *offset);
        tile.end_ivs.sort_by_key(|(offset, _)| *offset);
    }

    tiles
}

#[cfg(test)]
mod tests {
    use super::*;

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
        assert_eq!(tiles[1].start_ivs[0], (50, 1)); // offset 50 from tile start

        assert_eq!(tiles[2].end_ivs.len(), 1);
        assert_eq!(tiles[2].end_ivs[0], (50, 1)); // offset 50 from tile start
    }

    #[test]
    fn test_index_sweep_spanning_interval() {
        let bounds = Interval::new(0, 500);
        // Interval spans tiles 1 and 2 completely (100-300)
        let intervals = vec![TaggedInterval::new(50, 350, 1)];
        let tiles = index_sweep(bounds, 100, &intervals);

        // Tile 0 (0-100): starts before, ends after -> but interval starts at 50
        // Actually interval 50-350:
        // - Tile 0: starts in tile at 50, doesn't span entire tile (end is at 350 > 100)
        // - Tile 1 (100-200): spans entire tile (50 < 100, 350 > 200)
        // - Tile 2 (200-300): spans entire tile (50 < 200, 350 > 300)
        // - Tile 3 (300-400): starts before (50 < 300), ends in tile at 350

        assert_eq!(tiles[0].start_ivs.len(), 1); // starts at offset 50
        assert!(tiles[0].running_counts.is_empty());

        assert_eq!(tiles[1].running_counts.len(), 1); // spans entire tile
        assert_eq!(tiles[1].running_counts[0], (1, 1));

        assert_eq!(tiles[2].running_counts.len(), 1); // spans entire tile
        assert_eq!(tiles[2].running_counts[0], (1, 1));

        assert_eq!(tiles[3].end_ivs.len(), 1); // ends at offset 50
        assert!(tiles[3].running_counts.is_empty());
    }

    #[test]
    fn test_index_sweep_multiple_intervals_same_sid() {
        let bounds = Interval::new(0, 300);
        let intervals = vec![
            TaggedInterval::new(0, 300, 1), // spans all tiles
            TaggedInterval::new(0, 300, 1), // same sid
        ];
        let tiles = index_sweep(bounds, 100, &intervals);

        // All three tiles should have running count of 2 for sid 1
        assert_eq!(tiles[0].running_counts, vec![(1, 2)]);
        assert_eq!(tiles[1].running_counts, vec![(1, 2)]);
        assert_eq!(tiles[2].running_counts, vec![(1, 2)]);
    }

    #[test]
    fn test_index_sweep_different_sids() {
        let bounds = Interval::new(0, 200);
        let intervals = vec![
            TaggedInterval::new(0, 200, 1),
            TaggedInterval::new(0, 200, 2),
        ];
        let tiles = index_sweep(bounds, 100, &intervals);

        // Both tiles should have running counts for both sids
        assert!(tiles[0].running_counts.contains(&(1, 1)));
        assert!(tiles[0].running_counts.contains(&(2, 1)));
    }

    #[test]
    fn test_index_sweep_boundary_interval() {
        let bounds = Interval::new(0, 200);
        // Interval exactly at tile boundary
        let intervals = vec![TaggedInterval::new(100, 200, 1)];
        let tiles = index_sweep(bounds, 100, &intervals);

        // Tile 0 (0-100): no overlap (end is 100, exclusive)
        assert!(tiles[0].is_empty());

        // Tile 1 (100-200): spans entire tile
        assert_eq!(tiles[1].running_counts.len(), 1);
    }

    #[test]
    fn test_index_sweep_point_interval() {
        // While we don't allow empty intervals, a 1-bp interval is fine
        let bounds = Interval::new(0, 100);
        let intervals = vec![TaggedInterval::new(50, 51, 1)];
        let tiles = index_sweep(bounds, 100, &intervals);

        // Single tile, interval starts and ends in it
        assert_eq!(tiles[0].start_ivs.len(), 1);
        assert_eq!(tiles[0].start_ivs[0], (50, 1));
        assert_eq!(tiles[0].end_ivs.len(), 1);
        assert_eq!(tiles[0].end_ivs[0], (51, 1));
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
