mod index;
mod query;

pub use index::index_sweep;
pub use query::query_sweep;

use crate::core::Interval;

/// A generic line-sweep that yields a slice of all items intersecting each tile.
///
/// Both index_sweep and query_sweep share the same head/scan pointer advancement
/// pattern. This function extracts that logic into a reusable abstraction.
///
/// # Arguments
/// * `num_tiles` - Number of tiles in the chunk
/// * `chunk_start` - Start coordinate of the chunk
/// * `tile_size` - Size of each tile
/// * `sorted_batch` - Items sorted by start coordinate
/// * `extract_interval` - Function to extract an Interval from an item
/// * `process_tile` - Callback for each tile with (tile_idx, tile_start, tile_end, active_items)
#[inline]
pub fn sweep_tiles<'a, T>(
    num_tiles: usize,
    chunk_start: u32,
    tile_size: u32,
    sorted_batch: &'a [T],
    extract_interval: impl Fn(&T) -> Interval,
    mut process_tile: impl FnMut(usize, u32, u32, &'a [T]),
) {
    let mut head = 0;

    for tile_idx in 0..num_tiles {
        let tile_start = chunk_start + tile_idx as u32 * tile_size;
        let tile_end = tile_start + tile_size;

        // Advance head past items that completely precede this tile
        while head < sorted_batch.len() && extract_interval(&sorted_batch[head]).end <= tile_start {
            head += 1;
        }

        // Find the upper bound of items that start before this tile ends
        let mut scan = head;
        while scan < sorted_batch.len() && extract_interval(&sorted_batch[scan]).start < tile_end {
            scan += 1;
        }

        // Yield the contiguous slice of intersecting items
        if head < scan {
            process_tile(tile_idx, tile_start, tile_end, &sorted_batch[head..scan]);
        }
    }
}
