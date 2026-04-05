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
    pub fn new(
        layer_id: u8,
        min_size: u32,
        max_size: u32,
        tile_size: u32,
        chunk_size: u32,
    ) -> Self {
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
    ///
    /// Key invariant: tile_size > max_interval_size for each layer.
    /// This guarantees no interval can span an entire tile, eliminating the need
    /// for running_counts in tiles.
    pub fn default_layers() -> Vec<LayerConfig> {
        // tile_size(k) = 16e3 * 2^k
        //      starting at 16k, a coommon convention.
        // max_size(k) = tile_size(k) / 8
        //      this could use adjustment; we want ~2k interval events per tile
        //      to fit in cpu cache.
        // min_size(k) = max_size(k-1) = max_size / 2;
        // chunk_size(k) = 250e6 / threads ~= 2^21
        //      since rykv can read partial-chunks, we can make them as large as possible
        //      granted that each thread can get their own.
        //      we will account for 128 cores.
        (0u8..=15)
            .map(|k| {
                let tile_size = 2u32.pow(14 + k as u32);
                let max_size = tile_size / 8;
                let min_size = max_size / 2;
                let chunk_size = 2u32.pow(21).max(1); // ensure chunk_size is at least 1
                LayerConfig::new(k, min_size, max_size, tile_size, chunk_size)
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
    /// Maximum coordinate seen across all data (legacy, use shard_max_coords for sharded DBs)
    #[serde(default)]
    pub max_coord: u32,
    /// List of shard names (e.g., chromosome names)
    #[serde(default)]
    pub shards: Vec<String>,
    /// Maximum coordinate per shard
    #[serde(default)]
    pub shard_max_coords: HashMap<String, u32>,
}

impl MasterHeader {
    pub fn new() -> Self {
        Self {
            sid_map: HashMap::new(),
            layer_configs: LayerConfig::default_layers(),
            max_coord: 0,
            shards: Vec::new(),
            shard_max_coords: HashMap::new(),
        }
    }

    /// Add a source with its metadata
    pub fn add_source(&mut self, sid: u32, metadata: SidMetadata) {
        self.sid_map.insert(sid, metadata);
    }

    /// Update max coordinate (global)
    pub fn update_max_coord(&mut self, coord: u32) {
        self.max_coord = self.max_coord.max(coord);
    }

    /// Update max coordinate for a specific shard
    pub fn update_shard_max_coord(&mut self, shard: &str, coord: u32) {
        let entry = self.shard_max_coords.entry(shard.to_string()).or_insert(0);
        *entry = (*entry).max(coord);
    }

    /// Add a shard if it doesn't exist
    pub fn add_shard(&mut self, shard: String) {
        if !self.shards.contains(&shard) {
            self.shards.push(shard);
        }
    }

    /// Get the number of sources
    pub fn num_sources(&self) -> usize {
        self.sid_map.len()
    }

    /// Get the number of shards
    pub fn num_shards(&self) -> usize {
        self.shards.len()
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
        assert_eq!(layers[0].tile_size, 16_384);
        assert!(layers[0].max_size < layers[0].tile_size);
        assert_eq!(layers[1].tile_size, 32_768);
        assert_eq!(layers[1].max_size, layers[0].max_size * 2);
    }

    #[test]
    fn test_master_header() {
        let mut header = MasterHeader::new();
        assert_eq!(header.num_sources(), 0);
        assert_eq!(header.num_shards(), 0);

        header.add_source(1, SidMetadata::new("a.bed".to_string(), 10, 1000));
        header.add_source(2, SidMetadata::new("b.bed".to_string(), 20, 2000));
        assert_eq!(header.num_sources(), 2);

        header.update_max_coord(50000);
        header.update_max_coord(30000);
        assert_eq!(header.max_coord, 50000);

        // Test shard functionality
        header.add_shard("chr1".to_string());
        header.add_shard("chr2".to_string());
        header.add_shard("chr1".to_string()); // duplicate, should not add
        assert_eq!(header.num_shards(), 2);

        header.update_shard_max_coord("chr1", 100000);
        header.update_shard_max_coord("chr1", 50000); // should keep 100000
        assert_eq!(header.shard_max_coords["chr1"], 100000);
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
