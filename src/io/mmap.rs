use std::fs::File;
use std::path::Path;

use memmap2::Mmap;
use rkyv::rancor::Error as RkyvError;
use thiserror::Error;

use super::header::ChunkHeader;
use crate::core::{ArchivedTile, Tile};

/// Errors that can occur during memory mapping
#[derive(Debug, Error)]
pub enum MmapError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Invalid chunk format: {0}")]
    InvalidFormat(String),

    #[error("Tile {tile_id} out of bounds (max: {max_tiles})")]
    TileOutOfBounds { tile_id: usize, max_tiles: usize },
}

/// A memory-mapped chunk file providing zero-copy access to tiles
pub struct MappedChunk {
    mmap: Mmap,
    header: ChunkHeader,
    /// Offset where tile data begins (after header)
    data_offset: usize,
}

impl MappedChunk {
    /// Open and memory-map a chunk file
    pub fn open(path: &Path) -> Result<Self, MmapError> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };

        // Read header length (first 4 bytes)
        if mmap.len() < 4 {
            return Err(MmapError::InvalidFormat("File too small".to_string()));
        }

        let header_len = u32::from_le_bytes([mmap[0], mmap[1], mmap[2], mmap[3]]) as usize;
        let header_start = 4;
        let header_end = header_start + header_len;

        if mmap.len() < header_end {
            return Err(MmapError::InvalidFormat(
                "Header length exceeds file size".to_string(),
            ));
        }

        // Deserialize header using rkyv
        let header_bytes = &mmap[header_start..header_end];
        let archived = rkyv::access::<rkyv::Archived<ChunkHeader>, RkyvError>(header_bytes)
            .map_err(|e| MmapError::InvalidFormat(format!("Failed to read header: {}", e)))?;

        // Deserialize the archived header to owned
        let header: ChunkHeader =
            rkyv::deserialize::<ChunkHeader, RkyvError>(archived).map_err(|e| {
                MmapError::InvalidFormat(format!("Failed to deserialize header: {}", e))
            })?;

        Ok(Self {
            mmap,
            header,
            data_offset: header_end,
        })
    }

    /// Create a MappedChunk from raw bytes (useful for testing)
    pub fn from_bytes(
        bytes: Vec<u8>,
        _header: ChunkHeader,
        _data_offset: usize,
    ) -> Result<Self, MmapError> {
        // For testing, we create a mmap from a file backed by the bytes
        // This is a simplified version - in real use, we'd read from disk
        let _ = bytes; // We can't actually create Mmap from bytes directly

        // Instead, we'll store the header and offset directly
        // This method signature is preserved for API compatibility
        // but the implementation differs
        Err(MmapError::InvalidFormat(
            "from_bytes not supported in production".to_string(),
        ))
    }

    /// Get the chunk header
    pub fn header(&self) -> &ChunkHeader {
        &self.header
    }

    /// Get the number of tiles in this chunk
    pub fn num_tiles(&self) -> usize {
        self.header.tile_offsets.len()
    }

    /// Access a tile by ID (zero-copy)
    pub fn get_tile(&self, tile_id: usize) -> Result<&ArchivedTile, MmapError> {
        if tile_id >= self.header.tile_offsets.len() {
            return Err(MmapError::TileOutOfBounds {
                tile_id,
                max_tiles: self.header.tile_offsets.len(),
            });
        }

        let offset = self.data_offset + self.header.tile_offsets[tile_id] as usize;

        // Determine tile end
        let end = if tile_id + 1 < self.header.tile_offsets.len() {
            self.data_offset + self.header.tile_offsets[tile_id + 1] as usize
        } else {
            self.mmap.len()
        };

        // TODO: this will never happen, make this a debug assert
        if end > self.mmap.len() {
            return Err(MmapError::InvalidFormat(format!(
                "Tile {} extends beyond file",
                tile_id
            )));
        }

        let tile_bytes = &self.mmap[offset..end];
        rkyv::access::<ArchivedTile, RkyvError>(tile_bytes).map_err(|e| {
            MmapError::InvalidFormat(format!("Failed to read tile {}: {}", tile_id, e))
        })
    }

    /// Deserialize a tile (owned copy)
    pub fn read_tile(&self, tile_id: usize) -> Result<Tile, MmapError> {
        let archived = self.get_tile(tile_id)?;
        rkyv::deserialize::<Tile, RkyvError>(archived).map_err(|e| {
            MmapError::InvalidFormat(format!("Failed to deserialize tile {}: {}", tile_id, e))
        })
    }

    /// Get the start coordinate of the chunk
    pub fn start_coord(&self) -> u32 {
        self.header.start_coord
    }

    /// Get the end coordinate of the chunk
    pub fn end_coord(&self) -> u32 {
        self.header.end_coord
    }

    /// Get the layer ID
    pub fn layer_id(&self) -> u8 {
        self.header.layer_id
    }

    /// Get the chunk ID
    pub fn chunk_id(&self) -> u32 {
        self.header.chunk_id
    }
}

/// Write a chunk to a file
pub fn write_chunk(path: &Path, header: &ChunkHeader, tiles: &[Tile]) -> Result<(), MmapError> {
    use std::io::Write;

    let mut file = File::create(path)?;

    // Serialize header
    let header_bytes = rkyv::to_bytes::<RkyvError>(header)
        .map_err(|e| MmapError::InvalidFormat(format!("Failed to serialize header: {}", e)))?;

    // Write header length and header
    file.write_all(&(header_bytes.len() as u32).to_le_bytes())?;
    file.write_all(&header_bytes)?;

    // Write each tile
    for tile in tiles {
        let tile_bytes = rkyv::to_bytes::<RkyvError>(tile)
            .map_err(|e| MmapError::InvalidFormat(format!("Failed to serialize tile: {}", e)))?;
        file.write_all(&tile_bytes)?;
    }

    Ok(())
}

/// Byte size of one archived `TaggedInterval`.
///
/// `TaggedInterval` is three `u32` fields (start, end, sid) with no padding —
/// rkyv archives primitive types at their natural size, so the archived form is
/// identical in layout and therefore identical in size.  Adding `N` intervals to a
/// tile grows its serialized byte count by exactly `N * ARCHIVED_INTERVAL_BYTES`.
const ARCHIVED_INTERVAL_BYTES: usize = std::mem::size_of::<crate::core::TaggedInterval>();

/// Merge new tile data into an existing chunk file in-place.
///
/// # Algorithm
///
/// Sizes are derived analytically: each additional interval grows the serialized tile
/// by exactly `ARCHIVED_INTERVAL_BYTES`, so new offsets are computed with no I/O.
///
/// After extending the file to the required size, working from the last tile backward
/// to the first, each tile's old bytes are shifted rightward to their new position via
/// `copy_within` (memmove semantics), then the old tile is deserialized from that
/// position, sorted-merged with the corresponding new intervals, and the merged tile
/// is written back.  Back-to-front ordering ensures each tile's source bytes are
/// intact when copied.
///
/// Finally the header is rewritten with the updated tile offsets (same byte length —
/// only values change, not the Vec length).
///
/// # Errors
/// Returns an error if the chunk file's tile count does not match `new_tiles`.
pub fn merge_chunk(path: &Path, new_tiles: &[Tile]) -> Result<(), MmapError> {
    use memmap2::MmapMut;

    let mapped = MappedChunk::open(path)?;
    let n_tiles = mapped.num_tiles();

    if n_tiles != new_tiles.len() {
        return Err(MmapError::InvalidFormat(format!(
            "tile count mismatch: existing={n_tiles}, new={}",
            new_tiles.len()
        )));
    }

    // Short-circuit when there are no new intervals anywhere.
    if new_tiles.iter().all(|t| t.is_empty()) {
        return Ok(());
    }

    let data_offset = mapped.data_offset;
    let old_file_len = mapped.mmap.len();
    let old_header = mapped.header().clone();

    // Compute old and merged tile sizes analytically.
    // Serialized Tile = fixed base + N * ARCHIVED_INTERVAL_BYTES, so adding M new
    // intervals grows the tile by exactly M * ARCHIVED_INTERVAL_BYTES.
    let mut old_sizes: Vec<usize> = Vec::with_capacity(n_tiles);
    let mut merged_sizes: Vec<usize> = Vec::with_capacity(n_tiles);

    for i in 0..n_tiles {
        let old_start = data_offset + old_header.tile_offsets[i] as usize;
        let old_end = if i + 1 < n_tiles {
            data_offset + old_header.tile_offsets[i + 1] as usize
        } else {
            old_file_len
        };
        let old_sz = old_end - old_start;
        old_sizes.push(old_sz);
        merged_sizes.push(old_sz + new_tiles[i].intervals.len() * ARCHIVED_INTERVAL_BYTES);
    }

    // Compute new tile offsets (relative to data_offset)
    let mut new_tile_offsets: Vec<u32> = Vec::with_capacity(n_tiles);
    let mut cursor = 0u32;
    for &sz in &merged_sizes {
        new_tile_offsets.push(cursor);
        cursor += sz as u32;
    }

    // Build new header (same structure, only tile_offsets values differ)
    let mut new_header = old_header.clone();
    new_header.tile_offsets = new_tile_offsets.clone();
    let new_header_bytes = rkyv::to_bytes::<RkyvError>(&new_header)
        .map_err(|e| MmapError::InvalidFormat(format!("header serialize: {e}")))?;

    let old_header_len = data_offset - 4; // 4-byte length prefix not included
    if new_header_bytes.len() != old_header_len {
        return Err(MmapError::InvalidFormat(format!(
            "header size changed: old={old_header_len} new={} (tile count invariant violated)",
            new_header_bytes.len()
        )));
    }

    let new_file_len = data_offset + merged_sizes.iter().sum::<usize>();

    // Drop read-only mmap before resizing
    drop(mapped);

    // --- Extend file ---
    {
        let f = std::fs::OpenOptions::new().write(true).open(path)?;
        f.set_len(new_file_len as u64)?;
    }

    // --- Pass 2: shift tiles back-to-front, then overwrite with merged data ---
    let file = std::fs::OpenOptions::new()
        .read(true)
        .write(true)
        .open(path)?;
    let mut mmap = unsafe { MmapMut::map_mut(&file)? };

    // Absolute positions
    let old_abs: Vec<usize> = old_header
        .tile_offsets
        .iter()
        .map(|&off| data_offset + off as usize)
        .collect();
    let new_abs: Vec<usize> = new_tile_offsets
        .iter()
        .map(|&off| data_offset + off as usize)
        .collect();

    for i in (0..n_tiles).rev() {
        let old_start = old_abs[i];
        let old_sz = old_sizes[i];
        let new_start = new_abs[i];

        // Shift old bytes rightward to their new position.
        // copy_within uses memmove semantics, so overlapping ranges are safe.
        mmap.copy_within(old_start..old_start + old_sz, new_start);

        // Deserialize the (now-shifted) old tile.
        let old_tile: Tile = {
            let tile_bytes = &mmap[new_start..new_start + old_sz];
            let archived = rkyv::access::<ArchivedTile, RkyvError>(tile_bytes)
                .map_err(|e| MmapError::InvalidFormat(format!("tile access: {e}")))?;
            rkyv::deserialize::<Tile, RkyvError>(archived)
                .map_err(|e| MmapError::InvalidFormat(format!("tile deserialize: {e}")))?
        };

        // Sorted-merge and write merged bytes at the new position.
        let merged = merge_tile(&old_tile, &new_tiles[i]);
        let merged_bytes = rkyv::to_bytes::<RkyvError>(&merged)
            .map_err(|e| MmapError::InvalidFormat(format!("merged tile serialize: {e}")))?;
        mmap[new_start..new_start + merged_bytes.len()].copy_from_slice(&merged_bytes);
    }

    // Rewrite header (same offset, same byte length, updated tile_offsets values)
    mmap[4..4 + new_header_bytes.len()].copy_from_slice(&new_header_bytes);

    mmap.flush()?;
    Ok(())
}

/// Sorted merge of two tiles' interval lists (merge-sort style).
///
/// Both `old.intervals` and `new_tile.intervals` are sorted by start coordinate.
/// The result is a new Tile with start_coord from `old` and all intervals merged.
fn merge_tile(old: &Tile, new_tile: &Tile) -> Tile {
    let mut out = Tile::new(old.start_coord);
    let (a, b) = (old.intervals.as_slice(), new_tile.intervals.as_slice());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        if a[i].iv.start <= b[j].iv.start {
            out.intervals.push(a[i]);
            i += 1;
        } else {
            out.intervals.push(b[j]);
            j += 1;
        }
    }
    out.intervals.extend_from_slice(&a[i..]);
    out.intervals.extend_from_slice(&b[j..]);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_write_and_read_chunk() {
        let mut header = ChunkHeader::new(2, 0, 0, 1000);

        // Create some tiles
        let mut tile0 = Tile::new(0);
        tile0
            .intervals
            .push(crate::core::TaggedInterval::new(10, 30, 2));
        tile0
            .intervals
            .push(crate::core::TaggedInterval::new(20, 40, 1));

        let mut tile1 = Tile::new(500);
        tile1
            .intervals
            .push(crate::core::TaggedInterval::new(550, 580, 3));

        let tiles = vec![tile0, tile1];

        // Calculate offsets
        let tile_bytes: Vec<_> = tiles
            .iter()
            .map(|t| rkyv::to_bytes::<RkyvError>(t).unwrap())
            .collect();

        let mut offset = 0u32;
        for bytes in &tile_bytes {
            header.tile_offsets.push(offset);
            offset += bytes.len() as u32;
        }

        // Write to temp file
        let file = NamedTempFile::new().unwrap();
        write_chunk(file.path(), &header, &tiles).unwrap();

        // Read back
        let mapped = MappedChunk::open(file.path()).unwrap();
        assert_eq!(mapped.num_tiles(), 2);
        assert_eq!(mapped.layer_id(), 2);
        assert_eq!(mapped.chunk_id(), 0);

        // Check tile data
        let read_tile0 = mapped.read_tile(0).unwrap();
        assert_eq!(read_tile0.start_coord, 0);
        assert_eq!(read_tile0.intervals.len(), 2);

        let read_tile1 = mapped.read_tile(1).unwrap();
        assert_eq!(read_tile1.start_coord, 500);
        assert_eq!(read_tile1.intervals.len(), 1);
    }

    #[test]
    fn test_merge_chunk_adds_intervals() {
        use crate::core::TaggedInterval;
        use crate::sweep::index_sweep;

        // Build an initial chunk: chunk 0 of layer 0, bounds [0, 1000), tile_size 500
        let chunk_bounds = crate::core::Interval::new(0, 1000);
        let tile_size = 500u32;

        let initial_intervals = vec![TaggedInterval::new(100, 200, 0)];
        let initial_tiles = index_sweep(chunk_bounds, tile_size, &initial_intervals);

        let mut header = ChunkHeader::new(0, 0, 0, 1000);
        let mut offset = 0u32;
        for tile in &initial_tiles {
            header.tile_offsets.push(offset);
            let b = rkyv::to_bytes::<RkyvError>(tile).unwrap();
            offset += b.len() as u32;
        }

        let file = NamedTempFile::new().unwrap();
        write_chunk(file.path(), &header, &initial_tiles).unwrap();

        // Merge new intervals (sid=1, different tile than sid=0)
        let new_intervals = vec![TaggedInterval::new(600, 700, 1)];
        let new_tiles = index_sweep(chunk_bounds, tile_size, &new_intervals);
        merge_chunk(file.path(), &new_tiles).unwrap();

        // Verify: tile 0 has sid=0 interval, tile 1 has sid=1 interval
        let mapped = MappedChunk::open(file.path()).unwrap();
        assert_eq!(mapped.num_tiles(), 2);

        let t0 = mapped.read_tile(0).unwrap();
        assert_eq!(t0.intervals.len(), 1);
        assert_eq!(t0.intervals[0].sid, 0);

        let t1 = mapped.read_tile(1).unwrap();
        assert_eq!(t1.intervals.len(), 1);
        assert_eq!(t1.intervals[0].sid, 1);
    }

    #[test]
    fn test_merge_chunk_sorted_merge() {
        use crate::core::TaggedInterval;
        use crate::sweep::index_sweep;

        let chunk_bounds = crate::core::Interval::new(0, 1000);
        let tile_size = 1000u32; // single tile

        // Initial: one interval at 200–300
        let initial = vec![TaggedInterval::new(200, 300, 0)];
        let initial_tiles = index_sweep(chunk_bounds, tile_size, &initial);

        let mut header = ChunkHeader::new(0, 0, 0, 1000);
        let mut offset = 0u32;
        for tile in &initial_tiles {
            header.tile_offsets.push(offset);
            let b = rkyv::to_bytes::<RkyvError>(tile).unwrap();
            offset += b.len() as u32;
        }

        let file = NamedTempFile::new().unwrap();
        write_chunk(file.path(), &header, &initial_tiles).unwrap();

        // Merge: two new intervals interleaved with the existing one
        let new_intervals = vec![
            TaggedInterval::new(100, 150, 1), // before existing
            TaggedInterval::new(400, 500, 2), // after existing
        ];
        let new_tiles = index_sweep(chunk_bounds, tile_size, &new_intervals);
        merge_chunk(file.path(), &new_tiles).unwrap();

        let mapped = MappedChunk::open(file.path()).unwrap();
        let t = mapped.read_tile(0).unwrap();
        assert_eq!(t.intervals.len(), 3);
        // Must be sorted by start coordinate
        assert!(
            t.intervals
                .windows(2)
                .all(|w| w[0].iv.start <= w[1].iv.start)
        );
        // Correct order: sid=1 (100), sid=0 (200), sid=2 (400)
        assert_eq!(t.intervals[0].sid, 1);
        assert_eq!(t.intervals[1].sid, 0);
        assert_eq!(t.intervals[2].sid, 2);
    }

    #[test]
    fn test_tile_out_of_bounds() {
        let header = ChunkHeader::new(0, 0, 0, 100);
        let tiles: Vec<Tile> = vec![];

        let file = NamedTempFile::new().unwrap();
        write_chunk(file.path(), &header, &tiles).unwrap();

        let mapped = MappedChunk::open(file.path()).unwrap();
        let result = mapped.get_tile(0);
        assert!(matches!(result, Err(MmapError::TileOutOfBounds { .. })));
    }
}
