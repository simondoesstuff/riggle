use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

use memmap2::Mmap;

use crate::core::Interval;
use crate::io::layer::LayerError;

/// A read-only memory-mapped jump table for O(1) fast-forward indexing.
///
/// `table[t]` stores the count of intervals with `start < t * tile_size`,
/// which equals the index of the first interval with `start >= t * tile_size`.
///
/// Entries are `u64` to support layer sizes beyond the u32 range.
pub struct MappedJumpTable {
    mmap: Mmap,
    pub tile_size: u32,
}

impl MappedJumpTable {
    /// Open a `.idx` file as a zero-copy jump table.
    pub fn open(path: &Path, tile_size: u32) -> Result<Self, LayerError> {
        let file = File::open(path)?;
        let size = file.metadata()?.len();
        if size % 8 != 0 {
            return Err(LayerError::InvalidIndexSize(size));
        }
        let mmap = unsafe { Mmap::map(&file)? };
        Ok(Self { mmap, tile_size })
    }

    /// Look up the approximate first interval index for a target coordinate.
    ///
    /// Returns a conservative (possibly too-small) index: all intervals before
    /// the true position have `start < target`, so the caller must scan forward.
    #[inline]
    pub fn lookup(&self, target: u32) -> usize {
        let num_tiles = self.mmap.len() / 8;
        let tile_idx = (target / self.tile_size) as usize;
        let effective_idx = tile_idx.min(num_tiles - 1);
        let off = effective_idx * 8;
        u64::from_le_bytes(self.mmap[off..off + 8].try_into().unwrap()) as usize
    }
}

/// Build a dense jump table from a sorted batch of intervals.
///
/// `table[t]` = index of the first interval with `start >= t * tile_size`
///            = count of intervals with `start < t * tile_size`
///
/// Uses a single linear scan to detect tile transitions, collecting sparse
/// carry-update checkpoints, then densifies via a loop over just those
/// checkpoints with vectorized `fill`.  Efficient when the batch is sparse
/// over coordinate space.
///
/// Always returns at least one entry (`[0]`) so the result can be safely
/// memory-mapped.
pub fn build_jump_table(sorted_ivs: &[Interval], tile_size: u32) -> Vec<u64> {
    assert!(tile_size > 0, "tile_size must be > 0");

    if sorted_ivs.is_empty() {
        return vec![0]; // single conservative entry; layer is empty
    }

    let max_start = sorted_ivs.last().unwrap().start;
    let max_tile = (max_start / tile_size) as usize;

    // Sparse pass: detect tile transitions and record carry-update checkpoints.
    //
    // When the tile index changes from `prev_tile` to `curr_tile` at interval
    // index `i`, we record sparse[prev_tile + 1] = i.  This encodes:
    //   "from tile prev_tile+1 onwards, the running count is i"
    // which equals partition_point(sorted_ivs, curr_tile * tile_size).
    //
    // The first tile needs no entry: the initial carry (0) is already correct
    // for all tiles 0..=first_tile (no intervals have start < first_tile * tile_size).
    let mut sparse: HashMap<usize, u64> = HashMap::new();
    let mut prev_tile: Option<usize> = None;
    for (i, iv) in sorted_ivs.iter().enumerate() {
        let t = (iv.start / tile_size) as usize;
        if Some(t) != prev_tile {
            if let Some(pt) = prev_tile {
                sparse.insert(pt + 1, i as u64);
            }
            prev_tile = Some(t);
        }
    }

    // Densify: iterate over sorted checkpoints and fill each segment with the
    // current carry using a vectorized slice::fill.  Far fewer iterations than
    // a tile-by-tile loop when the batch is sparse over coordinate space.
    let mut checkpoints: Vec<(usize, u64)> = sparse.into_iter().collect();
    checkpoints.sort_unstable_by_key(|&(t, _)| t);

    let mut dense = vec![0u64; max_tile + 1];
    let mut carry = 0u64;
    let mut start = 0usize;
    for (t, v) in checkpoints {
        dense[start..t].fill(carry);
        carry = v;
        start = t;
    }
    dense[start..=max_tile].fill(carry);
    dense
}

/// Write a jump table to a `.idx` file (flat little-endian `u64[]`, no header).
pub fn write_jump_table(path: &Path, table: &[u64]) -> Result<(), LayerError> {
    let mut bytes = Vec::with_capacity(table.len() * 8);
    for v in table {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    std::fs::write(path, bytes)?;
    Ok(())
}

/// Merge new intervals' jump table contribution into the `.idx` file.
///
/// If the `.idx` file is absent, it is created from scratch.  Otherwise the
/// existing table is loaded and merged with the new batch via element-wise sum:
///
/// `merged[t] = old[t] + (new_local[t] ?? 0)`
///
/// This is exact for tiles within `new_local`'s range and conservative
/// (under-counts) for tiles beyond it; correctness is preserved by the
/// forward scan at query time.
pub fn extend_jump_table(
    idx_path: &Path,
    new_sorted: &[Interval],
    tile_size: u32,
) -> Result<(), LayerError> {
    if new_sorted.is_empty() {
        return Ok(());
    }

    // Load existing table, or start from empty if the file is absent.
    let old_table: Vec<u64> = if idx_path.exists() {
        let bytes = std::fs::read(idx_path)?;
        if bytes.len() % 8 != 0 {
            return Err(LayerError::InvalidIndexSize(bytes.len() as u64));
        }
        bytes
            .chunks_exact(8)
            .map(|b| u64::from_le_bytes(b.try_into().unwrap()))
            .collect()
    } else {
        Vec::new()
    };

    let new_local = build_jump_table(new_sorted, tile_size);

    // Element-wise sum, extending to the longer of the two tables.
    let merged_len = old_table.len().max(new_local.len());
    let mut merged = vec![0u64; merged_len];
    for t in 0..merged_len {
        merged[t] = old_table.get(t).copied().unwrap_or(0)
            + new_local.get(t).copied().unwrap_or(0);
    }

    write_jump_table(idx_path, &merged)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Interval;
    use tempfile::NamedTempFile;

    fn iv(start: u32, end: u32, sid: u32) -> Interval {
        Interval::new(start, end, sid)
    }

    #[test]
    fn test_build_empty() {
        // Empty input → single conservative entry.
        let table = build_jump_table(&[], 64);
        assert_eq!(table, vec![0]);
    }

    #[test]
    fn test_build_single_tile() {
        // tile_size=100; all intervals in tile 0 (start < 100)
        let ivs = vec![iv(10, 20, 0), iv(50, 60, 1), iv(80, 90, 2)];
        let table = build_jump_table(&ivs, 100);
        assert_eq!(table.len(), 1);
        assert_eq!(table[0], 0); // tile 0: count with start < 0 = 0
    }

    #[test]
    fn test_build_multiple_tiles() {
        // tile_size=100; ivs at starts 10, 150, 250
        let ivs = vec![iv(10, 20, 0), iv(150, 160, 1), iv(250, 260, 2)];
        let table = build_jump_table(&ivs, 100);
        assert_eq!(table.len(), 3);
        assert_eq!(table[0], 0); // start < 0: 0 intervals
        assert_eq!(table[1], 1); // start < 100: 1 interval (start=10)
        assert_eq!(table[2], 2); // start < 200: 2 intervals (start=10,150)
    }

    #[test]
    fn test_fill_forward_empty_tiles() {
        // tile_size=100; ivs jump from tile 0 to tile 3 (tiles 1 and 2 empty)
        let ivs = vec![iv(10, 20, 0), iv(350, 360, 1)];
        let table = build_jump_table(&ivs, 100);
        assert_eq!(table.len(), 4); // tiles 0..=3
        assert_eq!(table[0], 0); // start < 0: 0
        assert_eq!(table[1], 1); // start < 100: 1
        assert_eq!(table[2], 1); // start < 200: fill-forward
        assert_eq!(table[3], 1); // start < 300: fill-forward
    }

    #[test]
    fn test_build_matches_partition_point() {
        // Verify against the naive partition_point approach.
        let ivs = vec![
            iv(5, 10, 0),
            iv(15, 20, 0),
            iv(105, 110, 0),
            iv(205, 210, 0),
            iv(305, 310, 0),
        ];
        let tile_size = 100u32;
        let table = build_jump_table(&ivs, tile_size);
        for (t, &entry) in table.iter().enumerate() {
            let boundary = t as u32 * tile_size;
            let expected = ivs.partition_point(|x| x.start < boundary) as u64;
            assert_eq!(entry, expected, "mismatch at tile {t}");
        }
    }

    #[test]
    fn test_write_read_roundtrip() {
        let file = NamedTempFile::new().unwrap();
        let table = vec![0u64, 5, 10, 15];
        write_jump_table(file.path(), &table).unwrap();

        let jt = MappedJumpTable::open(file.path(), 100).unwrap();
        assert_eq!(jt.lookup(0), 0);
        assert_eq!(jt.lookup(100), 5);
        assert_eq!(jt.lookup(200), 10);
        assert_eq!(jt.lookup(300), 15);
    }

    #[test]
    fn test_lookup_beyond_table() {
        let file = NamedTempFile::new().unwrap();
        let table = vec![0u64, 5, 10];
        write_jump_table(file.path(), &table).unwrap();

        let jt = MappedJumpTable::open(file.path(), 100).unwrap();
        // tile_idx=5 is beyond table length 3 → clamps to last entry (10)
        assert_eq!(jt.lookup(500), 10);
    }

    #[test]
    fn test_extend_creates_when_absent() {
        // When no .idx exists, extend creates it from scratch.
        let dir = tempfile::tempdir().unwrap();
        let idx_path = dir.path().join("layer_0.idx");
        extend_jump_table(&idx_path, &[iv(10, 20, 0), iv(110, 120, 1)], 100).unwrap();
        assert!(idx_path.exists());
        let jt = MappedJumpTable::open(&idx_path, 100).unwrap();
        assert_eq!(jt.lookup(0), 0);
        assert_eq!(jt.lookup(100), 1);
    }

    #[test]
    fn test_extend_merge_sum() {
        let file = NamedTempFile::new().unwrap();
        // Old table: [0, 3] — 0 intervals before tile 0, 3 before tile 1.
        write_jump_table(file.path(), &[0u64, 3]).unwrap();

        // New batch: 2 intervals in tile 0, 1 in tile 1.
        let new_ivs = vec![iv(10, 20, 0), iv(20, 30, 1), iv(110, 120, 2)];
        extend_jump_table(file.path(), &new_ivs, 100).unwrap();

        let jt = MappedJumpTable::open(file.path(), 100).unwrap();
        assert_eq!(jt.lookup(0), 0);   // 0 + 0
        assert_eq!(jt.lookup(100), 5); // 3 + 2
    }

    #[test]
    fn test_extend_new_local_longer_than_old() {
        let file = NamedTempFile::new().unwrap();
        // Old table covers only 2 tiles.
        write_jump_table(file.path(), &[0u64, 4]).unwrap();

        // New batch has intervals reaching tile 3.
        let new_ivs = vec![iv(10, 20, 0), iv(110, 120, 1), iv(210, 220, 2), iv(310, 320, 3)];
        extend_jump_table(file.path(), &new_ivs, 100).unwrap();

        let jt = MappedJumpTable::open(file.path(), 100).unwrap();
        assert_eq!(jt.lookup(0), 0);   // 0 + 0
        assert_eq!(jt.lookup(100), 5); // 4 + 1
        assert_eq!(jt.lookup(200), 2); // 0 + 2
        assert_eq!(jt.lookup(300), 3); // 0 + 3
    }
}
