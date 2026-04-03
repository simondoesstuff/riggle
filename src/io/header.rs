use std::collections::HashMap;

use rkyv::{Archive, Deserialize, Serialize};
use serde::{Deserialize as SerdeDeserialize, Serialize as SerdeSerialize};

/// Metadata for a source ID (database entry)
#[derive(Debug, Clone, Archive, Serialize, Deserialize, SerdeSerialize, SerdeDeserialize)]
pub struct SidMetadata {
    /// Original filename or identifier
    pub name: String,
    /// Total number of intervals from this source
    pub interval_count: u32,
    /// Total base pairs covered
    pub total_bases: u64,
}

impl SidMetadata {
    pub fn new(name: String, interval_count: u32, total_bases: u64) -> Self {
        Self {
            name,
            interval_count,
            total_bases,
        }
    }
}

/// Configuration for a layer (exponential size bucket)
///
/// Chunks within a layer are spaced at predictable intervals, enabling O(1) lookup:
/// - chunk_id = coordinate / chunk_size
/// - chunk file path: layer_{layer_id}/chunk_{chunk_id}.bin
///
/// This predictable spacing means given any coordinate, we can immediately compute
/// which chunk file to access without scanning the filesystem.
#[derive(Debug, Clone, Archive, Serialize, Deserialize, SerdeSerialize, SerdeDeserialize)]
pub struct LayerConfig {
    /// Layer ID (log2 of size range)
    pub layer_id: u8,
    /// Minimum interval size for this layer
    pub min_size: u32,
    /// Maximum interval size for this layer
    pub max_size: u32,
    /// Size of each tile within chunks
    pub tile_size: u32,
    /// Size of each chunk (collection of tiles) - determines O(1) lookup spacing
    pub chunk_size: u32,
}

impl LayerConfig {
    /// Create a new layer configuration
    pub fn new(layer_id: u8, min_size: u32, max_size: u32, tile_size: u32, chunk_size: u32) -> Self {
        Self {
            layer_id,
            min_size,
            max_size,
            tile_size,
            chunk_size,
        }
    }

    /// Create default layer configurations for a given maximum coordinate
    /// Layers are based on powers of 2, with tile/chunk sizes scaled accordingly
    pub fn default_layers() -> Vec<LayerConfig> {
        // Each layer handles intervals of size [2^layer_id, 2^(layer_id+1))
        // Tile size is 2^(layer_id+2) to ensure each interval touches at most a few tiles
        // Chunk size is 2^(layer_id+10) to balance memory usage
        (0u8..=20)
            .map(|layer_id| {
                let min_size = 1u32 << layer_id;
                let max_size = 1u32 << (layer_id + 1);
                let tile_size = 1u32.checked_shl((layer_id + 4) as u32).unwrap_or(1 << 24);
                let chunk_size = 1u32.checked_shl((layer_id + 12) as u32).unwrap_or(1 << 28);
                LayerConfig::new(layer_id, min_size, max_size, tile_size, chunk_size)
            })
            .collect()
    }

    /// Get the number of tiles in a chunk
    pub fn tiles_per_chunk(&self) -> u32 {
        self.chunk_size / self.tile_size
    }
}

/// Master header containing global metadata
#[derive(Debug, Clone, Archive, Serialize, Deserialize, SerdeSerialize, SerdeDeserialize)]
pub struct MasterHeader {
    /// Mapping from source ID to metadata
    pub sid_map: HashMap<u32, SidMetadata>,
    /// Configuration for each layer
    pub layer_configs: Vec<LayerConfig>,
    /// Maximum coordinate seen across all data
    pub max_coord: u32,
}

impl MasterHeader {
    pub fn new() -> Self {
        Self {
            sid_map: HashMap::new(),
            layer_configs: LayerConfig::default_layers(),
            max_coord: 0,
        }
    }

    /// Add a source with its metadata
    pub fn add_source(&mut self, sid: u32, metadata: SidMetadata) {
        self.sid_map.insert(sid, metadata);
    }

    /// Update max coordinate
    pub fn update_max_coord(&mut self, coord: u32) {
        self.max_coord = self.max_coord.max(coord);
    }

    /// Get the number of sources
    pub fn num_sources(&self) -> usize {
        self.sid_map.len()
    }
}

impl Default for MasterHeader {
    fn default() -> Self {
        Self::new()
    }
}

/// Header for a chunk file (stored at the beginning of each chunk)
#[derive(Debug, Clone, Archive, Serialize, Deserialize)]
pub struct ChunkHeader {
    /// Layer ID this chunk belongs to
    pub layer_id: u8,
    /// Chunk ID within the layer
    pub chunk_id: u32,
    /// Start coordinate of this chunk
    pub start_coord: u32,
    /// End coordinate of this chunk
    pub end_coord: u32,
    /// Byte offsets for each tile within the chunk (tile_id -> offset)
    pub tile_offsets: Vec<u32>,
}

impl ChunkHeader {
    pub fn new(layer_id: u8, chunk_id: u32, start_coord: u32, end_coord: u32) -> Self {
        Self {
            layer_id,
            chunk_id,
            start_coord,
            end_coord,
            tile_offsets: Vec::new(),
        }
    }

    /// Get the number of tiles in this chunk
    pub fn num_tiles(&self) -> usize {
        self.tile_offsets.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sid_metadata() {
        let meta = SidMetadata::new("test.bed".to_string(), 100, 50000);
        assert_eq!(meta.name, "test.bed");
        assert_eq!(meta.interval_count, 100);
        assert_eq!(meta.total_bases, 50000);
    }

    #[test]
    fn test_layer_config() {
        let config = LayerConfig::new(5, 32, 64, 128, 4096);
        assert_eq!(config.layer_id, 5);
        assert_eq!(config.tiles_per_chunk(), 32);
    }

    #[test]
    fn test_default_layers() {
        let layers = LayerConfig::default_layers();
        assert!(!layers.is_empty());

        // Check first few layers
        assert_eq!(layers[0].min_size, 1);
        assert_eq!(layers[0].max_size, 2);
        assert_eq!(layers[1].min_size, 2);
        assert_eq!(layers[1].max_size, 4);
    }

    #[test]
    fn test_master_header() {
        let mut header = MasterHeader::new();
        assert_eq!(header.num_sources(), 0);

        header.add_source(1, SidMetadata::new("a.bed".to_string(), 10, 1000));
        header.add_source(2, SidMetadata::new("b.bed".to_string(), 20, 2000));
        assert_eq!(header.num_sources(), 2);

        header.update_max_coord(50000);
        header.update_max_coord(30000);
        assert_eq!(header.max_coord, 50000);
    }

    #[test]
    fn test_chunk_header() {
        let header = ChunkHeader::new(3, 5, 0, 4096);
        assert_eq!(header.layer_id, 3);
        assert_eq!(header.chunk_id, 5);
        assert_eq!(header.num_tiles(), 0);
    }

    #[test]
    fn test_header_serialization() {
        use rkyv::rancor::Error;

        let mut master = MasterHeader::new();
        master.add_source(1, SidMetadata::new("test.bed".to_string(), 100, 5000));
        master.update_max_coord(100000);

        // Serialize
        let bytes = rkyv::to_bytes::<Error>(&master).unwrap();

        // Deserialize
        let archived = rkyv::access::<ArchivedMasterHeader, Error>(&bytes).unwrap();
        assert_eq!(archived.max_coord, 100000);
        assert_eq!(archived.sid_map.len(), 1);
    }
}
