use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use rayon::prelude::*;
use thiserror::Error;

use crate::core::TaggedInterval;
use crate::io::{parse_bed_file, BedParseError, MappedChunk, MasterHeader};
use crate::matrix::{allocate_dense_accumulator, condense_to_sparse, merge_sparse, SparseMatrix};
use crate::sweep::query_sweep;

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
    /// Path to query BED file
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

/// Result of a query operation
#[derive(Debug)]
pub struct QueryResult {
    /// Sparse matrix of intersection counts (rows = queries, cols = db sources)
    pub counts: SparseMatrix,
    /// Query metadata (index -> name mapping)
    pub query_names: Vec<String>,
    /// Database source metadata
    pub db_sources: HashMap<u32, String>,
}

/// Execute a query against the database
pub fn query_database(config: &QueryConfig) -> Result<QueryResult, QueryError> {
    // Load master header
    let header_path = config.db_path.join("header.json");
    if !header_path.exists() {
        return Err(QueryError::DatabaseNotFound(config.db_path.clone()));
    }

    let header_content = fs::read_to_string(&header_path)?;
    let master_header: MasterHeader = serde_json::from_str(&header_content)
        .map_err(|e| QueryError::InvalidDatabase(e.to_string()))?;

    // Parse query file
    let queries = parse_bed_file(&config.query_path, 0)?;
    if queries.is_empty() {
        // Return empty result
        let counts = sprs::CsMat::empty(sprs::CompressedStorage::CSR, master_header.num_sources());
        return Ok(QueryResult {
            counts,
            query_names: Vec::new(),
            db_sources: master_header
                .sid_map
                .iter()
                .map(|(k, v)| (*k, v.name.clone()))
                .collect(),
        });
    }

    // Sort queries and track original indices
    let mut indexed_queries: Vec<(usize, TaggedInterval)> = queries
        .into_iter()
        .enumerate()
        .collect();
    indexed_queries.sort_by_key(|(_, iv)| iv.iv.start);

    let num_queries = indexed_queries.len();
    let num_sources = master_header.num_sources();

    // Compute which chunks to query based on query coordinates (O(1) per query)
    // Each interval goes to exactly one layer based on its size
    let chunk_tasks = compute_chunk_tasks(&config.db_path, &indexed_queries, &master_header);

    // Process chunks in parallel, each with thread-local accumulator
    let sparse_matrices: Vec<SparseMatrix> = chunk_tasks
        .par_iter()
        .filter_map(|(layer_id, chunk_id, path)| {
            process_chunk(
                path,
                *layer_id,
                *chunk_id,
                &master_header,
                &indexed_queries,
                num_queries,
                num_sources,
            )
            .ok()
        })
        .collect();

    // Merge all sparse matrices
    let final_counts = tree_reduce_sparse(sparse_matrices, num_queries, num_sources);

    // Build query names
    let query_names: Vec<String> = (0..num_queries)
        .map(|i| format!("query_{}", i))
        .collect();

    Ok(QueryResult {
        counts: final_counts,
        query_names,
        db_sources: master_header
            .sid_map
            .iter()
            .map(|(k, v)| (*k, v.name.clone()))
            .collect(),
    })
}

/// Compute chunk file paths to query based on query coordinates
/// This provides O(1) lookup per query - we compute exactly which chunks need to be accessed
/// based on coordinates rather than scanning the filesystem.
fn compute_chunk_tasks(
    db_path: &Path,
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
            let path = chunk_path(db_path, layer_id, chunk_id);
            if path.exists() {
                Some((layer_id, chunk_id, path))
            } else {
                None
            }
        })
        .collect()
}

/// Construct the path to a chunk file (O(1) - predictable location)
#[inline]
fn chunk_path(db_path: &Path, layer_id: u8, chunk_id: u32) -> PathBuf {
    db_path
        .join(format!("layer_{}", layer_id))
        .join(format!("chunk_{}.bin", chunk_id))
}

/// Process a single chunk
fn process_chunk(
    path: &Path,
    layer_id: u8,
    _chunk_id: u32,
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
        return Ok(sprs::CsMat::empty(sprs::CompressedStorage::CSR, num_sources));
    }

    // Allocate thread-local accumulator
    let (mut dense, mut mask) = allocate_dense_accumulator(num_queries, num_sources);

    // Get tile size from layer config
    let tile_size = master_header
        .layer_configs
        .get(layer_id as usize)
        .map(|c| c.tile_size)
        .unwrap_or(1024);

    // Run query sweep
    query_sweep(&mapped, tile_size, &relevant_queries, &mut dense, &mut mask);

    // Condense to sparse
    Ok(condense_to_sparse(&dense, &mask))
}

/// Tree-reduce sparse matrices
fn tree_reduce_sparse(
    matrices: Vec<SparseMatrix>,
    _num_queries: usize,
    num_sources: usize,
) -> SparseMatrix {
    if matrices.is_empty() {
        return sprs::CsMat::empty(sprs::CompressedStorage::CSR, num_sources);
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
    use crate::tasks::build::{build_database, BuildConfig};
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
        create_test_bed(
            input_dir.path(),
            "b.bed",
            "chr1\t150\t250\n",
        );

        let build_config = BuildConfig::new(
            input_dir.path().to_path_buf(),
            db_dir.path().to_path_buf(),
        );
        build_database(&build_config).unwrap();

        // Create query file
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
}
