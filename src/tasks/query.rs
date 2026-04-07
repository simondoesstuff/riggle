use std::cell::RefCell;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use rayon::prelude::*;
use thiserror::Error;
use tqdm::{Style, tqdm};

use crate::core::TaggedInterval;
use crate::io::{BedParseError, MappedChunk, MasterHeader, parse_bed_file};
use crate::matrix::{
    BitwiseMask, DenseMatrix, SparseMatrix, condense_to_sparse, merge_sparse, zero_flagged_regions,
};
use crate::sweep::query_sweep;

// Thread-local buffer for query accumulation
// Reuses DenseMatrix and BitwiseMask across chunk tasks on the same thread
thread_local! {
    static THREAD_LOCAL_BUFFER: RefCell<Option<(DenseMatrix, BitwiseMask)>> = const { RefCell::new(None) };
}

/// Errors that can occur during query
#[derive(Debug, Error)]
pub enum QueryError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("BED parse error: {0}")]
    BedParse(#[from] BedParseError),

    #[error("Database not found at {0}")]
    DatabaseNotFound(PathBuf),

    #[error("Invalid database format: {0}")]
    InvalidDatabase(String),

    #[error("Mmap error: {0}")]
    Mmap(#[from] crate::io::mmap::MmapError),
}

/// Configuration for query execution
#[derive(Debug, Clone)]
pub struct QueryConfig {
    /// Path to database directory
    pub db_path: PathBuf,
    /// Path to query BED file or directory containing BED files
    pub query_path: PathBuf,
    /// Number of threads (default: all available)
    pub num_threads: Option<usize>,
}

impl QueryConfig {
    pub fn new(db_path: PathBuf, query_path: PathBuf) -> Self {
        Self {
            db_path,
            query_path,
            num_threads: None,
        }
    }
}

/// Metadata for a query source file
#[derive(Debug, Clone)]
pub struct QuerySource {
    /// File name
    pub name: String,
    /// Starting row index in results matrix
    pub start_idx: usize,
    /// Number of intervals
    pub count: usize,
}

/// Result of a query operation
#[derive(Debug)]
pub struct QueryResult {
    /// Sparse matrix of intersection counts (rows = queries, cols = db sources)
    pub counts: SparseMatrix,
    /// Query metadata (index -> name mapping)
    pub query_names: Vec<String>,
    /// Query source files (when querying with a directory)
    pub query_sources: Vec<QuerySource>,
    /// Database source metadata
    pub db_sources: HashMap<u32, String>,
}

/// Execute a query against the database
///
/// Accepts either a single BED file or a directory containing multiple BED files.
/// When given a directory, each file is treated as a separate query source.
///
/// Queries are processed per-shard: only shards present in both the query AND
/// the database are queried. Results are aggregated across all shards.
pub fn query_database(config: &QueryConfig) -> Result<QueryResult, QueryError> {
    // Load master header
    let header_path = config.db_path.join("header.json");
    if !header_path.exists() {
        return Err(QueryError::DatabaseNotFound(config.db_path.clone()));
    }

    let header_content = fs::read_to_string(&header_path)?;
    let master_header: MasterHeader = serde_json::from_str(&header_content)
        .map_err(|e| QueryError::InvalidDatabase(e.to_string()))?;

    // Parse query file(s) - handles both single file and directory
    let parsed = parse_query_path(&config.query_path)?;

    if parsed.total_count == 0 {
        // Return empty result
        let counts = sprs::CsMat::empty(sprs::CompressedStorage::CSR, master_header.num_sources());
        return Ok(QueryResult {
            counts,
            query_names: Vec::new(),
            query_sources: Vec::new(),
            db_sources: master_header
                .sid_map
                .iter()
                .map(|(k, v)| (*k, v.name.clone()))
                .collect(),
        });
    }

    let num_queries = parsed.total_count;
    let num_sources = master_header.num_sources();

    // Build a set of database shards for fast lookup
    let db_shards: std::collections::HashSet<&str> =
        master_header.shards.iter().map(|s| s.as_str()).collect();

    // Collect all sparse matrices from all shards
    let mut all_sparse_matrices: Vec<SparseMatrix> = Vec::new();

    // Process each shard that exists in both query and database
    for (shard, shard_queries) in tqdm(parsed.shard_intervals.into_iter())
        .style(Style::Balloon)
        .desc(Some("Processing shards"))
    {
        // Skip if shard doesn't exist in database
        if !db_shards.contains(shard.as_str()) {
            continue;
        }

        if shard_queries.is_empty() {
            continue;
        }

        // Sort queries by start coordinate while preserving global indices
        let mut indexed_queries: Vec<(usize, TaggedInterval)> = shard_queries.clone();
        indexed_queries.sort_by_key(|(_, iv)| iv.iv.start);

        // Compute which chunks to query for this shard
        let chunk_tasks = compute_chunk_tasks_for_shard(
            &config.db_path,
            &shard,
            &indexed_queries,
            &master_header,
        );

        // Process chunks in parallel
        let shard_matrices: Vec<SparseMatrix> = chunk_tasks
            .par_iter()
            .filter_map(|(layer_id, _chunk_id, path)| {
                process_chunk(
                    path,
                    *layer_id,
                    &master_header,
                    &indexed_queries,
                    num_queries,
                    num_sources,
                )
                .ok()
            })
            .collect();

        all_sparse_matrices.extend(shard_matrices);
    }

    // Merge all sparse matrices from all shards
    let final_counts = tree_reduce_sparse(all_sparse_matrices, num_queries, num_sources);

    Ok(QueryResult {
        counts: final_counts,
        query_names: parsed.query_names,
        query_sources: parsed.query_sources,
        db_sources: master_header
            .sid_map
            .iter()
            .map(|(k, v)| (*k, v.name.clone()))
            .collect(),
    })
}

/// Parsed query result with shard grouping
struct ParsedQueries {
    /// Intervals grouped by shard, with global indices preserved
    /// Each entry is (global_index, interval)
    shard_intervals: HashMap<String, Vec<(usize, TaggedInterval)>>,
    /// Query source metadata
    query_sources: Vec<QuerySource>,
    /// Flat list of query names (for result indexing)
    query_names: Vec<String>,
    /// Total number of query intervals
    total_count: usize,
}

/// Parse query BED file(s) from a path (file or directory)
fn parse_query_path(path: &Path) -> Result<ParsedQueries, QueryError> {
    // Collect BED files to process
    let bed_files: Vec<PathBuf> = if path.is_dir() {
        fs::read_dir(path)?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| p.extension().map(|ext| ext == "bed").unwrap_or(false))
            .collect()
    } else {
        vec![path.to_path_buf()]
    };

    if bed_files.is_empty() {
        return Ok(ParsedQueries {
            shard_intervals: HashMap::new(),
            query_sources: Vec::new(),
            query_names: Vec::new(),
            total_count: 0,
        });
    }

    // First pass: collect all intervals with global indices
    let mut shard_intervals: HashMap<String, Vec<(usize, TaggedInterval)>> = HashMap::new();
    let mut query_sources = Vec::new();
    let mut query_names = Vec::new();
    let mut global_idx = 0;

    for bed_path in bed_files {
        let file_shards = parse_bed_file(&bed_path, 0)?;

        // Count total intervals in this file
        let file_count: usize = file_shards.values().map(|v| v.len()).sum();

        if file_count == 0 {
            continue;
        }

        let name = bed_path
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| "query".to_string());

        query_names.extend((0..file_count).map(|i| format!("{}:{}", name, i)));
        query_sources.push(QuerySource {
            name,
            start_idx: global_idx,
            count: file_count,
        });

        // Merge into global shard map with global indices
        for (shard, intervals) in file_shards {
            let indexed: Vec<(usize, TaggedInterval)> = intervals
                .into_iter()
                .map(|iv| {
                    let idx = global_idx;
                    global_idx += 1;
                    (idx, iv)
                })
                .collect();
            shard_intervals.entry(shard).or_default().extend(indexed);
        }
    }

    Ok(ParsedQueries {
        shard_intervals,
        query_sources,
        query_names,
        total_count: global_idx,
    })
}

/// Compute chunk file paths to query for a specific shard
/// This provides O(1) lookup per query - we compute exactly which chunks need to be accessed
/// based on coordinates rather than scanning the filesystem.
fn compute_chunk_tasks_for_shard(
    db_path: &Path,
    shard: &str,
    indexed_queries: &[(usize, TaggedInterval)],
    master_header: &MasterHeader,
) -> Vec<(u8, u32, PathBuf)> {
    use std::collections::HashSet;

    let mut tasks: HashSet<(u8, u32)> = HashSet::new();

    // For each query, determine which layer and chunks it needs to access
    // Queries must check ALL layers since database intervals are partitioned by size
    for layer_config in &master_header.layer_configs {
        let layer_id = layer_config.layer_id;
        let chunk_size = layer_config.chunk_size;

        for (_, query) in indexed_queries {
            // Compute chunk range that this query intersects
            let start_chunk = query.iv.start / chunk_size;
            let end_chunk = query.iv.end.saturating_sub(1) / chunk_size;

            // Add all chunks in range
            for chunk_id in start_chunk..=end_chunk {
                tasks.insert((layer_id, chunk_id));
            }
        }
    }

    // Convert to paths, filtering to only existing files
    tasks
        .into_iter()
        .filter_map(|(layer_id, chunk_id)| {
            let path = chunk_path(db_path, shard, layer_id, chunk_id);
            if path.exists() {
                Some((layer_id, chunk_id, path))
            } else {
                None
            }
        })
        .collect()
}

/// Construct the path to a chunk file within a shard (O(1) - predictable location)
#[inline]
fn chunk_path(db_path: &Path, shard: &str, layer_id: u8, chunk_id: u32) -> PathBuf {
    db_path
        .join(shard)
        .join(format!("layer_{}", layer_id))
        .join(format!("chunk_{}.bin", chunk_id))
}

/// Process a single chunk using thread-local buffer reuse
fn process_chunk(
    path: &Path,
    layer_id: u8,
    master_header: &MasterHeader,
    indexed_queries: &[(usize, TaggedInterval)],
    num_queries: usize,
    num_sources: usize,
) -> Result<SparseMatrix, QueryError> {
    let mapped = MappedChunk::open(path)?;

    // Filter queries that intersect this chunk
    let chunk_start = mapped.start_coord();
    let chunk_end = mapped.end_coord();

    let relevant_queries: Vec<(usize, TaggedInterval)> = indexed_queries
        .iter()
        .filter(|(_, iv)| iv.iv.start < chunk_end && iv.iv.end > chunk_start)
        .cloned()
        .collect();

    if relevant_queries.is_empty() {
        return Ok(sprs::CsMat::empty(
            sprs::CompressedStorage::CSR,
            num_sources,
        ));
    }

    // Get tile size from layer config
    let tile_size = master_header
        .layer_configs
        .get(layer_id as usize)
        .map(|c| c.tile_size)
        .unwrap_or(1024);

    // Use thread-local buffer (reused across chunks on same thread)
    THREAD_LOCAL_BUFFER.with(|cell| {
        let mut borrow = cell.borrow_mut();

        // Initialize or resize buffers as needed
        let (dense, mask) = borrow.get_or_insert_with(|| {
            (
                DenseMatrix::new(num_queries, num_sources),
                BitwiseMask::new(num_queries, num_sources),
            )
        });

        // Resize if dimensions changed (only reallocates if larger)
        if dense.num_rows() != num_queries || dense.num_cols() != num_sources {
            dense.resize_and_zero(num_queries, num_sources);
            mask.resize_and_clear(num_queries, num_sources);
        }

        // Run query sweep
        query_sweep(&mapped, tile_size, &relevant_queries, dense, mask);

        // Condense to sparse result
        let result = condense_to_sparse(dense, mask);

        // Clear only flagged regions for next use (more efficient than zeroing all)
        zero_flagged_regions(dense, mask);
        mask.clear_all();

        Ok(result)
    })
}

/// Tree-reduce sparse matrices
fn tree_reduce_sparse(
    matrices: Vec<SparseMatrix>,
    num_queries: usize,
    num_sources: usize,
) -> SparseMatrix {
    if matrices.is_empty() {
        // Return properly sized empty matrix
        return sprs::CsMat::new(
            (num_queries, num_sources),
            vec![0; num_queries + 1],
            Vec::new(),
            Vec::new(),
        );
    }

    if matrices.len() == 1 {
        return matrices.into_iter().next().unwrap();
    }

    // Parallel tree reduction
    let mut current = matrices;

    while current.len() > 1 {
        let pairs: Vec<_> = current
            .chunks(2)
            .map(|chunk| {
                if chunk.len() == 2 {
                    merge_sparse(&chunk[0], &chunk[1])
                } else {
                    chunk[0].clone()
                }
            })
            .collect();
        current = pairs;
    }

    current.into_iter().next().unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tasks::build::{BuildConfig, build_database};
    use std::io::Write;
    use tempfile::TempDir;

    fn create_test_bed(dir: &Path, name: &str, content: &str) {
        let path = dir.join(name);
        let mut file = fs::File::create(path).unwrap();
        file.write_all(content.as_bytes()).unwrap();
    }

    #[test]
    fn test_query_database_basic() {
        let input_dir = TempDir::new().unwrap();
        let db_dir = TempDir::new().unwrap();
        let query_dir = TempDir::new().unwrap();

        // Create database
        create_test_bed(
            input_dir.path(),
            "a.bed",
            "chr1\t100\t200\nchr1\t300\t400\n",
        );
        create_test_bed(input_dir.path(), "b.bed", "chr1\t150\t250\n");

        let build_config =
            BuildConfig::new(input_dir.path().to_path_buf(), db_dir.path().to_path_buf());
        build_database(&build_config).unwrap();

        // Create query file (must include shard/chromosome)
        create_test_bed(
            query_dir.path(),
            "query.bed",
            "chr1\t100\t200\n", // Should overlap with a.bed interval 1 and b.bed
        );

        let query_config = QueryConfig::new(
            db_dir.path().to_path_buf(),
            query_dir.path().join("query.bed"),
        );

        let result = query_database(&query_config).unwrap();

        assert_eq!(result.counts.rows(), 1);
        assert_eq!(result.db_sources.len(), 2);
    }

    #[test]
    fn test_query_database_multi_shard() {
        let input_dir = TempDir::new().unwrap();
        let db_dir = TempDir::new().unwrap();
        let query_dir = TempDir::new().unwrap();

        // Create database with multiple shards
        create_test_bed(
            input_dir.path(),
            "a.bed",
            "chr1\t100\t200\nchr2\t300\t400\n",
        );

        let build_config =
            BuildConfig::new(input_dir.path().to_path_buf(), db_dir.path().to_path_buf());
        build_database(&build_config).unwrap();

        // Create query file targeting both shards
        create_test_bed(
            query_dir.path(),
            "query.bed",
            "chr1\t100\t200\nchr2\t300\t400\n",
        );

        let query_config = QueryConfig::new(
            db_dir.path().to_path_buf(),
            query_dir.path().join("query.bed"),
        );

        let result = query_database(&query_config).unwrap();

        // Should have 2 query intervals
        assert_eq!(result.counts.rows(), 2);
    }

    #[test]
    fn test_query_database_shard_not_in_db() {
        let input_dir = TempDir::new().unwrap();
        let db_dir = TempDir::new().unwrap();
        let query_dir = TempDir::new().unwrap();

        // Create database with chr1 only
        create_test_bed(input_dir.path(), "a.bed", "chr1\t100\t200\n");

        let build_config =
            BuildConfig::new(input_dir.path().to_path_buf(), db_dir.path().to_path_buf());
        build_database(&build_config).unwrap();

        // Query with chr2 (not in database) - should be silently skipped
        create_test_bed(query_dir.path(), "query.bed", "chr2\t100\t200\n");

        let query_config = QueryConfig::new(
            db_dir.path().to_path_buf(),
            query_dir.path().join("query.bed"),
        );

        let result = query_database(&query_config).unwrap();

        // Should have 1 query interval but no overlaps
        assert_eq!(result.counts.rows(), 1);
    }

    #[test]
    fn test_query_database_not_found() {
        let query_dir = TempDir::new().unwrap();
        create_test_bed(query_dir.path(), "query.bed", "chr1\t100\t200\n");

        let config = QueryConfig::new(
            PathBuf::from("/nonexistent/path"),
            query_dir.path().join("query.bed"),
        );

        let result = query_database(&config);
        assert!(matches!(result, Err(QueryError::DatabaseNotFound(_))));
    }

    #[test]
    fn test_tree_reduce_sparse() {
        use crate::matrix::{BitwiseMask, DenseMatrix};

        let mut matrices = Vec::new();

        for i in 0..5 {
            let mut dense = DenseMatrix::new(2, 3);
            let mut mask = BitwiseMask::new(2, 3);
            dense.set(0, 0, i + 1);
            mask.flag(0, 0);
            matrices.push(condense_to_sparse(&dense, &mask));
        }

        let result = tree_reduce_sparse(matrices, 2, 3);

        // Should sum to 1+2+3+4+5 = 15
        assert_eq!(result.get(0, 0), Some(&15));
    }

    #[test]
    fn test_query_database_batch() {
        let input_dir = TempDir::new().unwrap();
        let db_dir = TempDir::new().unwrap();
        let query_dir = TempDir::new().unwrap();

        // Create database
        create_test_bed(
            input_dir.path(),
            "a.bed",
            "chr1\t100\t200\nchr1\t300\t400\n",
        );

        let build_config =
            BuildConfig::new(input_dir.path().to_path_buf(), db_dir.path().to_path_buf());
        build_database(&build_config).unwrap();

        // Create multiple query files
        create_test_bed(query_dir.path(), "query1.bed", "chr1\t100\t200\n");
        create_test_bed(
            query_dir.path(),
            "query2.bed",
            "chr1\t300\t400\nchr1\t350\t450\n",
        );

        // Query with directory
        let query_config =
            QueryConfig::new(db_dir.path().to_path_buf(), query_dir.path().to_path_buf());

        let result = query_database(&query_config).unwrap();

        // Should have 3 total query intervals (1 from query1, 2 from query2)
        assert_eq!(result.counts.rows(), 3);
        assert_eq!(result.query_sources.len(), 2);

        // Check source metadata
        let total_count: usize = result.query_sources.iter().map(|s| s.count).sum();
        assert_eq!(total_count, 3);
    }
}
