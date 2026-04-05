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

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_write_and_read_chunk() {
        let mut header = ChunkHeader::new(2, 0, 0, 1000);

        // Create some tiles
        let mut tile0 = Tile::new(0);
        tile0.start_ivs.push(crate::core::OffsetSid::new(10, 2));
        tile0.start_ivs.push(crate::core::OffsetSid::new(20, 1));

        let mut tile1 = Tile::new(500);
        tile1.end_ivs.push(crate::core::OffsetSid::new(50, 3));

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
        assert_eq!(read_tile0.start_ivs.len(), 2);

        let read_tile1 = mapped.read_tile(1).unwrap();
        assert_eq!(read_tile1.start_coord, 500);
        assert_eq!(read_tile1.end_ivs.len(), 1);
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
