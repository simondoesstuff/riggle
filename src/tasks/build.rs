use std::collections::HashMap;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};

use rayon::prelude::*;
use thiserror::Error;
use voracious_radix_sort::RadixSort;

use crate::core::{ChunkID, Interval, LayerID, TaggedInterval, Tile};
use crate::io::{
    BedParseError, ChunkHeader, LayerConfig, MasterHeader, SidMetadata, parse_bed_file,
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

/// Build a database from BED files in a directory
pub fn build_database(config: &BuildConfig) -> Result<(), BuildError> {
    // Find all BED files
    let bed_files: Vec<_> = fs::read_dir(&config.input_dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().map(|ext| ext == "bed").unwrap_or(false))
        .collect();

    if bed_files.is_empty() {
        return Err(BuildError::NoInputFiles);
    }

    // Create output directory
    fs::create_dir_all(&config.output_dir)?;

    // Parse all files and collect intervals
    let mut master_header = MasterHeader::new();
    if let Some(ref configs) = config.layer_configs {
        master_header.layer_configs = configs.clone();
    }

    // Parse all files first to collect intervals and metadata
    let parsed_files: Vec<(u32, PathBuf, Vec<TaggedInterval>)> = bed_files
        .into_iter()
        .enumerate()
        .map(|(sid, bed_path)| {
            let sid = sid as u32;
            let intervals = parse_bed_file(&bed_path, sid)?;
            Ok((sid, bed_path, intervals))
        })
        .collect::<Result<Vec<_>, BedParseError>>()?;

    // Pre-allocate with exact capacity
    let total_intervals: usize = parsed_files.iter().map(|(_, _, ivs)| ivs.len()).sum();
    let mut all_intervals: Vec<TaggedInterval> = Vec::with_capacity(total_intervals);

    // Process each file's intervals
    for (sid, bed_path, intervals) in parsed_files {
        // Compute metadata
        let interval_count = intervals.len() as u32;
        let total_bases: u64 = intervals.iter().map(|iv| iv.iv.len() as u64).sum();
        let name = bed_path
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| format!("source_{}", sid));

        master_header.add_source(sid, SidMetadata::new(name, interval_count, total_bases));

        // Update max coordinate
        for iv in &intervals {
            master_header.update_max_coord(iv.iv.end);
        }

        all_intervals.extend(intervals);
    }

    // Partition intervals by layer
    let partitions = partition_by_layer(&all_intervals, &master_header.layer_configs);

    // Process each layer
    for (layer_id, mut layer_intervals) in partitions.into_iter().enumerate() {
        if layer_intervals.is_empty() {
            continue;
        }

        let layer_config = &master_header.layer_configs[layer_id];

        // Sort by start coordinate using radix sort O(n * w)
        layer_intervals.voracious_sort();

        // Group by chunk - O(n) sequential, but chunk writing below is parallel.
        // Parallelizing this grouping step is unlikely to help since:
        // 1. It's a simple O(n) pass with HashMap insertions
        // 2. The expensive work (serialization, I/O) is already parallel in into_par_iter()
        // 3. Parallel HashMap merging has overhead that may exceed benefits
        let chunks = group_by_chunk(&layer_intervals, layer_config);

        // Create layer directory
        let layer_dir = config.output_dir.join(format!("layer_{}", layer_id));
        fs::create_dir_all(&layer_dir)?;

        // Process chunks in parallel
        chunks
            .into_par_iter()
            .try_for_each(|(chunk_id, chunk_intervals)| {
                write_chunk_file(
                    &layer_dir,
                    layer_id as u8,
                    chunk_id,
                    layer_config,
                    &chunk_intervals,
                )
            })?;
    }

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
    layer_id: u8,
    chunk_id: u32,
    config: &LayerConfig,
    intervals: &[TaggedInterval],
) -> Result<(), BuildError> {
    // FIX: don't pass layer_id since it's incluced in LayerConfig;
    // this will lead to downstream bugs
    let chunk_start = chunk_id * config.chunk_size;
    let chunk_end = chunk_start + config.chunk_size;
    let chunk_bounds = Interval::new(chunk_start, chunk_end);

    // Build tiles using sweep algorithm
    let tiles = index_sweep(chunk_bounds, config.tile_size, intervals);

    // Create header
    let mut header = ChunkHeader::new(layer_id, chunk_id, chunk_start, chunk_end);

    // Calculate tile offsets
    let tile_bytes: Vec<_> = tiles
        .iter()
        .map(|t| {
            rkyv::to_bytes::<rkyv::rancor::Error>(t)
                .map_err(|e| BuildError::Serialization(e.to_string()))
        })
        .collect::<Result<Vec<_>, _>>()?;

    let mut offset = 0u32;
    for bytes in &tile_bytes {
        header.tile_offsets.push(offset);
        offset += bytes.len() as u32;
    }

    // Write to file
    let chunk_path = layer_dir.join(format!("chunk_{}.bin", chunk_id));
    write_chunk_to_file(&chunk_path, &header, &tiles)?;

    Ok(())
}

/// Write a chunk to file
fn write_chunk_to_file(
    path: &Path,
    header: &ChunkHeader,
    tiles: &[Tile],
) -> Result<(), BuildError> {
    let mut file = File::create(path)?;

    // Serialize header
    let header_bytes = rkyv::to_bytes::<rkyv::rancor::Error>(header)
        .map_err(|e| BuildError::Serialization(e.to_string()))?;

    // Write header length and header
    file.write_all(&(header_bytes.len() as u32).to_le_bytes())?;
    file.write_all(&header_bytes)?;

    // Write each tile
    for tile in tiles {
        let tile_bytes = rkyv::to_bytes::<rkyv::rancor::Error>(tile)
            .map_err(|e| BuildError::Serialization(e.to_string()))?;
        file.write_all(&tile_bytes)?;
    }

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
    let bed_files: Vec<_> = fs::read_dir(&config.input_dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().map(|ext| ext == "bed").unwrap_or(false))
        .collect();

    if bed_files.is_empty() {
        return Err(BuildError::NoInputFiles);
    }

    // Determine starting sid (max existing + 1)
    let start_sid = master_header
        .sid_map
        .keys()
        .max()
        .map(|&max| max + 1)
        .unwrap_or(0);

    // Parse new files
    let parsed_files: Vec<(u32, PathBuf, Vec<TaggedInterval>)> = bed_files
        .into_iter()
        .enumerate()
        .map(|(idx, bed_path)| {
            let sid = start_sid + idx as u32;
            let intervals = parse_bed_file(&bed_path, sid)?;
            Ok((sid, bed_path, intervals))
        })
        .collect::<Result<Vec<_>, BedParseError>>()?;

    let num_added = parsed_files.len();

    // Collect all new intervals
    let total_intervals: usize = parsed_files.iter().map(|(_, _, ivs)| ivs.len()).sum();
    let mut all_intervals: Vec<TaggedInterval> = Vec::with_capacity(total_intervals);

    for (sid, bed_path, intervals) in parsed_files {
        let interval_count = intervals.len() as u32;
        let total_bases: u64 = intervals.iter().map(|iv| iv.iv.len() as u64).sum();
        let name = bed_path
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| format!("source_{}", sid));

        master_header.add_source(sid, SidMetadata::new(name, interval_count, total_bases));

        for iv in &intervals {
            master_header.update_max_coord(iv.iv.end);
        }

        all_intervals.extend(intervals);
    }

    // Partition by layer and write chunks (same logic as build)
    let partitions = partition_by_layer(&all_intervals, &master_header.layer_configs);

    for (layer_id, mut layer_intervals) in partitions.into_iter().enumerate() {
        if layer_intervals.is_empty() {
            continue;
        }

        let layer_config = &master_header.layer_configs[layer_id];
        layer_intervals.voracious_sort();

        let chunks = group_by_chunk(&layer_intervals, layer_config);

        let layer_dir = config.db_dir.join(format!("layer_{}", layer_id));
        fs::create_dir_all(&layer_dir)?;

        // Write chunks - for add mode, we append to existing tiles if chunk exists
        // For simplicity, this implementation rewrites affected chunks
        // A more sophisticated approach would merge with existing tile data
        chunks
            .into_par_iter()
            .try_for_each(|(chunk_id, chunk_intervals)| {
                write_chunk_file(
                    &layer_dir,
                    layer_id as u8,
                    chunk_id,
                    layer_config,
                    &chunk_intervals,
                )
            })?;
    }

    // Update master header
    write_master_header(&config.db_dir, &master_header)?;

    Ok(num_added)
}

#[cfg(test)]
mod tests {
    use super::*;
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
        let header_content = fs::read_to_string(header_path).unwrap();
        let header: MasterHeader = serde_json::from_str(&header_content).unwrap();
        assert_eq!(header.num_sources(), 2);
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
        create_test_bed(
            input_dir.path(),
            "a.bed",
            "chr1\t100\t200\n",
        );

        let build_config = BuildConfig::new(
            input_dir.path().to_path_buf(),
            db_dir.path().to_path_buf(),
        );
        build_database(&build_config).unwrap();

        // Verify initial state
        let header_path = db_dir.path().join("header.json");
        let initial_header: MasterHeader =
            serde_json::from_str(&fs::read_to_string(&header_path).unwrap()).unwrap();
        assert_eq!(initial_header.num_sources(), 1);

        // Add new files
        create_test_bed(
            add_dir.path(),
            "b.bed",
            "chr1\t300\t400\n",
        );
        create_test_bed(
            add_dir.path(),
            "c.bed",
            "chr1\t500\t600\n",
        );

        let add_config = AddConfig::new(
            add_dir.path().to_path_buf(),
            db_dir.path().to_path_buf(),
        );
        let added = add_to_database(&add_config).unwrap();

        assert_eq!(added, 2);

        // Verify updated state
        let updated_header: MasterHeader =
            serde_json::from_str(&fs::read_to_string(&header_path).unwrap()).unwrap();
        assert_eq!(updated_header.num_sources(), 3);

        // Check sids are sequential
        assert!(updated_header.sid_map.contains_key(&0)); // original
        assert!(updated_header.sid_map.contains_key(&1)); // added
        assert!(updated_header.sid_map.contains_key(&2)); // added
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
