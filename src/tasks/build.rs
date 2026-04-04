use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use rayon::prelude::*;
use thiserror::Error;
use voracious_radix_sort::RadixSort;

use crate::core::{ChunkID, Interval, LayerID, TaggedInterval};
use crate::io::{
    BedParseError, ChunkHeader, LayerConfig, MasterHeader, SidMetadata,
    mmap::MmapError, parse_bed_file, write_chunk,
};
use crate::sweep::index_sweep;

/// Errors that can occur during database build
#[derive(Debug, Error)]
pub enum BuildError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("BED parse error: {0}")]
    BedParse(#[from] BedParseError),

    #[error("No input files found")]
    NoInputFiles,

    #[error("Serialization error: {0}")]
    Serialization(String),

    #[error("Database not found at {0}")]
    DatabaseNotFound(PathBuf),

    #[error("Invalid database: {0}")]
    InvalidDatabase(String),

    #[error("Chunk write error: {0}")]
    ChunkWrite(#[from] MmapError),
}

/// Configuration for database build
#[derive(Debug, Clone)]
pub struct BuildConfig {
    /// Input directory containing BED files
    pub input_dir: PathBuf,
    /// Output directory for database
    pub output_dir: PathBuf,
    /// Layer configurations (use default if None)
    pub layer_configs: Option<Vec<LayerConfig>>,
}

impl BuildConfig {
    pub fn new(input_dir: PathBuf, output_dir: PathBuf) -> Self {
        Self {
            input_dir,
            output_dir,
            layer_configs: None,
        }
    }
}

/// Configuration for adding files to existing database
#[derive(Debug, Clone)]
pub struct AddConfig {
    /// Input directory containing BED files to add
    pub input_dir: PathBuf,
    /// Path to existing database directory
    pub db_dir: PathBuf,
}

impl AddConfig {
    pub fn new(input_dir: PathBuf, db_dir: PathBuf) -> Self {
        Self { input_dir, db_dir }
    }
}

// =============================================================================
// Shared Pipeline Helpers
// =============================================================================

/// Find all BED files in a directory
fn find_bed_files(dir: &Path) -> Result<Vec<PathBuf>, BuildError> {
    let files: Vec<_> = fs::read_dir(dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().map(|ext| ext == "bed").unwrap_or(false))
        .collect();

    if files.is_empty() {
        return Err(BuildError::NoInputFiles);
    }

    Ok(files)
}

/// Parse BED files, update header metadata, and return intervals grouped by shard
fn ingest_bed_files(
    bed_files: Vec<PathBuf>,
    start_sid: u32,
    header: &mut MasterHeader,
) -> Result<HashMap<String, Vec<TaggedInterval>>, BuildError> {
    // Parse all files first to collect intervals grouped by shard
    let parsed_files: Vec<(u32, PathBuf, HashMap<String, Vec<TaggedInterval>>)> = bed_files
        .into_iter()
        .enumerate()
        .map(|(idx, bed_path)| {
            let sid = start_sid + idx as u32;
            let shard_intervals = parse_bed_file(&bed_path, sid)?;
            Ok((sid, bed_path, shard_intervals))
        })
        .collect::<Result<Vec<_>, BedParseError>>()?;

    // Global shard map to aggregate across files
    let mut all_shards: HashMap<String, Vec<TaggedInterval>> = HashMap::new();

    // Process each file's intervals
    for (sid, bed_path, shard_intervals) in parsed_files {
        // Compute metadata across all shards for this file
        let interval_count: u32 = shard_intervals.values().map(|v| v.len() as u32).sum();
        let total_bases: u64 = shard_intervals
            .values()
            .flat_map(|v| v.iter())
            .map(|iv| iv.iv.len() as u64)
            .sum();
        let name = bed_path
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| format!("source_{}", sid));

        header.add_source(sid, SidMetadata::new(name, interval_count, total_bases));

        // Process each shard's intervals
        for (shard, intervals) in shard_intervals {
            // Update shard max coordinate
            for iv in &intervals {
                header.update_shard_max_coord(&shard, iv.iv.end);
                header.update_max_coord(iv.iv.end); // Keep global max for backwards compat
            }

            // Add shard to header
            header.add_shard(shard.clone());

            // Aggregate intervals by shard
            all_shards.entry(shard).or_default().extend(intervals);
        }
    }

    Ok(all_shards)
}

/// Partition intervals by layer, sort, and write chunks to disk for a single shard
fn flush_shard_to_disk(
    shard: &str,
    intervals: &[TaggedInterval],
    header: &MasterHeader,
    output_dir: &Path,
) -> Result<(), BuildError> {
    let shard_dir = output_dir.join(shard);
    let partitions = partition_by_layer(intervals, &header.layer_configs);

    for (layer_id, mut layer_intervals) in partitions.into_iter().enumerate() {
        if layer_intervals.is_empty() {
            continue;
        }

        let layer_config = &header.layer_configs[layer_id];

        // Sort by start coordinate using radix sort O(n * w)
        layer_intervals.voracious_sort();

        // Group by chunk - O(n) sequential, but chunk writing below is parallel
        let chunks = group_by_chunk(&layer_intervals, layer_config);

        // Create layer directory within shard
        let layer_dir = shard_dir.join(format!("layer_{}", layer_id));
        fs::create_dir_all(&layer_dir)?;

        // Process chunks in parallel
        chunks
            .into_par_iter()
            .try_for_each(|(chunk_id, chunk_intervals)| {
                write_chunk_file(&layer_dir, chunk_id, layer_config, &chunk_intervals)
            })?;
    }

    Ok(())
}

/// Partition intervals by shard, layer, and write chunks to disk
fn flush_intervals_to_disk(
    shard_intervals: &HashMap<String, Vec<TaggedInterval>>,
    header: &MasterHeader,
    output_dir: &Path,
) -> Result<(), BuildError> {
    // Process each shard (could parallelize at this level too if desired)
    for (shard, intervals) in shard_intervals {
        if intervals.is_empty() {
            continue;
        }
        flush_shard_to_disk(shard, intervals, header, output_dir)?;
    }

    Ok(())
}

// =============================================================================
// Public API
// =============================================================================

/// Build a database from BED files in a directory
pub fn build_database(config: &BuildConfig) -> Result<(), BuildError> {
    let bed_files = find_bed_files(&config.input_dir)?;

    // Create output directory
    fs::create_dir_all(&config.output_dir)?;

    // Initialize master header
    let mut master_header = MasterHeader::new();
    if let Some(ref configs) = config.layer_configs {
        master_header.layer_configs = configs.clone();
    }

    // Ingest BED files starting from sid 0 (returns grouped by shard)
    let shard_intervals = ingest_bed_files(bed_files, 0, &mut master_header)?;

    // Write intervals to disk (per shard)
    flush_intervals_to_disk(&shard_intervals, &master_header, &config.output_dir)?;

    // Write master header
    write_master_header(&config.output_dir, &master_header)?;

    Ok(())
}

/// Partition intervals by layer based on their size.
///
/// IMPORTANT: Each interval is assigned to exactly ONE layer based on its length.
/// The layer is determined by log2(length), so intervals of similar sizes are grouped.
/// An interval may span multiple tiles and chunks within its layer, but never
/// crosses layer boundaries. This is critical for the query algorithm's correctness.
fn partition_by_layer(
    intervals: &[TaggedInterval],
    layer_configs: &[LayerConfig],
) -> Vec<Vec<TaggedInterval>> {
    let mut partitions: Vec<Vec<TaggedInterval>> = vec![Vec::new(); layer_configs.len()];

    for iv in intervals {
        // Assign to exactly one layer based on interval length
        let layer_id = LayerID::from_length(iv.iv.len());
        let layer_idx = (layer_id.0 as usize).min(partitions.len() - 1);
        partitions[layer_idx].push(*iv);
    }

    partitions
}

/// Group intervals by chunk, duplicating boundary-crossers
fn group_by_chunk(
    intervals: &[TaggedInterval],
    config: &LayerConfig,
) -> HashMap<u32, Vec<TaggedInterval>> {
    let mut chunks: HashMap<u32, Vec<TaggedInterval>> = HashMap::new();

    for iv in intervals {
        let start_chunk = ChunkID::from_coord(iv.iv.start, config.chunk_size);
        let end_chunk = ChunkID::from_coord(iv.iv.end.saturating_sub(1), config.chunk_size);

        // Add to all chunks the interval touches
        for chunk_id in start_chunk.0..=end_chunk.0 {
            chunks.entry(chunk_id).or_default().push(*iv);
        }
    }

    chunks
}

/// Write a chunk file
fn write_chunk_file(
    layer_dir: &Path,
    chunk_id: u32,
    config: &LayerConfig,
    intervals: &[TaggedInterval],
) -> Result<(), BuildError> {
    let chunk_start = chunk_id * config.chunk_size;
    let chunk_end = chunk_start + config.chunk_size;
    let chunk_bounds = Interval::new(chunk_start, chunk_end);

    // Build tiles using sweep algorithm
    let tiles = index_sweep(chunk_bounds, config.tile_size, intervals);

    // Create header with tile offsets
    let mut header = ChunkHeader::new(config.layer_id, chunk_id, chunk_start, chunk_end);

    // Calculate tile offsets for header
    let mut offset = 0u32;
    for tile in &tiles {
        header.tile_offsets.push(offset);
        let tile_bytes = rkyv::to_bytes::<rkyv::rancor::Error>(tile)
            .map_err(|e| BuildError::Serialization(e.to_string()))?;
        offset += tile_bytes.len() as u32;
    }

    // Write to file using consolidated write_chunk from io module
    let chunk_path = layer_dir.join(format!("chunk_{}.bin", chunk_id));
    write_chunk(&chunk_path, &header, &tiles)?;

    Ok(())
}

/// Write master header to file
fn write_master_header(output_dir: &Path, header: &MasterHeader) -> Result<(), BuildError> {
    let header_path = output_dir.join("header.json");
    let json = serde_json::to_string_pretty(header)
        .map_err(|e| BuildError::Serialization(e.to_string()))?;
    fs::write(header_path, json)?;
    Ok(())
}

/// Add BED files to an existing database
///
/// Returns the number of sources added.
pub fn add_to_database(config: &AddConfig) -> Result<usize, BuildError> {
    // Load existing master header
    let header_path = config.db_dir.join("header.json");
    if !header_path.exists() {
        return Err(BuildError::DatabaseNotFound(config.db_dir.clone()));
    }

    let header_content = fs::read_to_string(&header_path)?;
    let mut master_header: MasterHeader = serde_json::from_str(&header_content)
        .map_err(|e| BuildError::InvalidDatabase(e.to_string()))?;

    // Find new BED files
    let bed_files = find_bed_files(&config.input_dir)?;
    let num_added = bed_files.len();

    // Determine starting sid (max existing + 1)
    let start_sid = master_header
        .sid_map
        .keys()
        .max()
        .map(|&max| max + 1)
        .unwrap_or(0);

    // Ingest BED files (returns grouped by shard)
    let shard_intervals = ingest_bed_files(bed_files, start_sid, &mut master_header)?;

    // Write intervals to disk (per shard)
    // Note: For simplicity, this rewrites affected chunks rather than merging
    flush_intervals_to_disk(&shard_intervals, &master_header, &config.db_dir)?;

    // Update master header
    write_master_header(&config.db_dir, &master_header)?;

    Ok(num_added)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write as IoWrite;
    use tempfile::TempDir;

    fn create_test_bed(dir: &Path, name: &str, content: &str) {
        let path = dir.join(name);
        let mut file = File::create(path).unwrap();
        file.write_all(content.as_bytes()).unwrap();
    }

    #[test]
    fn test_build_database_basic() {
        let input_dir = TempDir::new().unwrap();
        let output_dir = TempDir::new().unwrap();

        // Create test BED files
        create_test_bed(
            input_dir.path(),
            "a.bed",
            "chr1\t100\t200\nchr1\t300\t400\n",
        );
        create_test_bed(
            input_dir.path(),
            "b.bed",
            "chr1\t150\t250\nchr1\t350\t450\n",
        );

        let config = BuildConfig::new(
            input_dir.path().to_path_buf(),
            output_dir.path().to_path_buf(),
        );

        build_database(&config).unwrap();

        // Check that master header was created
        let header_path = output_dir.path().join("header.json");
        assert!(header_path.exists());

        // Read and verify header
        let header_content = fs::read_to_string(&header_path).unwrap();
        let header: MasterHeader = serde_json::from_str(&header_content).unwrap();
        assert_eq!(header.num_sources(), 2);
        assert_eq!(header.num_shards(), 1);
        assert!(header.shards.contains(&"chr1".to_string()));

        // Check shard directory was created
        assert!(output_dir.path().join("chr1").exists());
    }

    #[test]
    fn test_build_database_multi_shard() {
        let input_dir = TempDir::new().unwrap();
        let output_dir = TempDir::new().unwrap();

        // Create test BED file with multiple chromosomes
        create_test_bed(
            input_dir.path(),
            "a.bed",
            "chr1\t100\t200\nchr2\t300\t400\nchr1\t500\t600\n",
        );

        let config = BuildConfig::new(
            input_dir.path().to_path_buf(),
            output_dir.path().to_path_buf(),
        );

        build_database(&config).unwrap();

        let header_path = output_dir.path().join("header.json");
        let header: MasterHeader =
            serde_json::from_str(&fs::read_to_string(&header_path).unwrap()).unwrap();

        assert_eq!(header.num_shards(), 2);
        assert!(header.shards.contains(&"chr1".to_string()));
        assert!(header.shards.contains(&"chr2".to_string()));

        // Check shard directories
        assert!(output_dir.path().join("chr1").exists());
        assert!(output_dir.path().join("chr2").exists());
    }

    #[test]
    fn test_build_database_no_files() {
        let input_dir = TempDir::new().unwrap();
        let output_dir = TempDir::new().unwrap();

        let config = BuildConfig::new(
            input_dir.path().to_path_buf(),
            output_dir.path().to_path_buf(),
        );

        let result = build_database(&config);
        assert!(matches!(result, Err(BuildError::NoInputFiles)));
    }

    #[test]
    fn test_partition_by_layer() {
        let intervals = vec![
            TaggedInterval::new(0, 10, 0),   // length 10, layer 3
            TaggedInterval::new(0, 100, 1),  // length 100, layer 6
            TaggedInterval::new(0, 1000, 2), // length 1000, layer 9
        ];

        let layers = LayerConfig::default_layers();
        let partitions = partition_by_layer(&intervals, &layers);

        // Check intervals are in correct layers
        assert!(partitions[3].iter().any(|iv| iv.sid == 0));
        assert!(partitions[6].iter().any(|iv| iv.sid == 1));
        assert!(partitions[9].iter().any(|iv| iv.sid == 2));
    }

    #[test]
    fn test_group_by_chunk() {
        let config = LayerConfig::new(0, 1, 2, 100, 1000);
        let intervals = vec![
            TaggedInterval::new(100, 200, 0),   // chunk 0
            TaggedInterval::new(900, 1100, 1),  // crosses chunks 0 and 1
            TaggedInterval::new(2000, 2100, 2), // chunk 2
        ];

        let chunks = group_by_chunk(&intervals, &config);

        assert!(chunks.get(&0).unwrap().iter().any(|iv| iv.sid == 0));
        assert!(chunks.get(&0).unwrap().iter().any(|iv| iv.sid == 1));
        assert!(chunks.get(&1).unwrap().iter().any(|iv| iv.sid == 1));
        assert!(chunks.get(&2).unwrap().iter().any(|iv| iv.sid == 2));
    }

    #[test]
    fn test_add_to_database() {
        let input_dir = TempDir::new().unwrap();
        let add_dir = TempDir::new().unwrap();
        let db_dir = TempDir::new().unwrap();

        // Create initial database
        create_test_bed(input_dir.path(), "a.bed", "chr1\t100\t200\n");

        let build_config =
            BuildConfig::new(input_dir.path().to_path_buf(), db_dir.path().to_path_buf());
        build_database(&build_config).unwrap();

        // Verify initial state
        let header_path = db_dir.path().join("header.json");
        let initial_header: MasterHeader =
            serde_json::from_str(&fs::read_to_string(&header_path).unwrap()).unwrap();
        assert_eq!(initial_header.num_sources(), 1);
        assert_eq!(initial_header.num_shards(), 1);

        // Add new files (including a new shard)
        create_test_bed(add_dir.path(), "b.bed", "chr1\t300\t400\n");
        create_test_bed(add_dir.path(), "c.bed", "chr2\t500\t600\n");

        let add_config = AddConfig::new(add_dir.path().to_path_buf(), db_dir.path().to_path_buf());
        let added = add_to_database(&add_config).unwrap();

        assert_eq!(added, 2);

        // Verify updated state
        let updated_header: MasterHeader =
            serde_json::from_str(&fs::read_to_string(&header_path).unwrap()).unwrap();
        assert_eq!(updated_header.num_sources(), 3);
        assert_eq!(updated_header.num_shards(), 2);

        // Check sids are sequential
        assert!(updated_header.sid_map.contains_key(&0)); // original
        assert!(updated_header.sid_map.contains_key(&1)); // added
        assert!(updated_header.sid_map.contains_key(&2)); // added

        // Check shard directories
        assert!(db_dir.path().join("chr1").exists());
        assert!(db_dir.path().join("chr2").exists());
    }

    #[test]
    fn test_add_to_nonexistent_database() {
        let input_dir = TempDir::new().unwrap();
        create_test_bed(input_dir.path(), "a.bed", "chr1\t100\t200\n");

        let config = AddConfig::new(
            input_dir.path().to_path_buf(),
            PathBuf::from("/nonexistent/db"),
        );

        let result = add_to_database(&config);
        assert!(matches!(result, Err(BuildError::DatabaseNotFound(_))));
    }
}
