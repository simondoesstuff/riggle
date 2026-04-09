use std::cell::RefCell;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use rayon::prelude::*;
use thiserror::Error;
use voracious_radix_sort::RadixSort;

use crate::core::Interval;
use crate::io::{BedParseError, LayerError, Meta, MetaError, MappedLayer, is_bed_file, parse_bed_file};
use crate::matrix::{DenseMatrix, SparseMatrix, condense_to_sparse_no_mask};
use crate::sweep::query_sweep;

// ---------------------------------------------------------------------------
// Thread-local scratch buffer — reused across sweep calls on the same thread.
// ---------------------------------------------------------------------------
thread_local! {
    static THREAD_LOCAL_BUFFER: RefCell<Option<DenseMatrix>> = const { RefCell::new(None) };
}

/// Errors from query execution
#[derive(Debug, Error)]
pub enum QueryError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("BED parse error: {0}")]
    BedParse(#[from] BedParseError),
    #[error("Meta error: {0}")]
    Meta(#[from] MetaError),
    #[error("Layer error: {0}")]
    Layer(#[from] LayerError),
    #[error("Database not found at {0}")]
    DatabaseNotFound(PathBuf),
}

/// Configuration for a query operation
#[derive(Debug, Clone)]
pub struct QueryConfig {
    pub db_path: PathBuf,
    pub query_path: PathBuf,
    /// Number of parallel threads; defaults to all available.
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

/// Per-file metadata returned with results
#[derive(Debug, Clone)]
pub struct QuerySource {
    pub name: String,
    pub count: usize,
}

/// Output of a query operation
#[derive(Debug)]
pub struct QueryResult {
    /// Sparse intersection count matrix: rows = Q_SIDs, cols = D_SIDs
    pub counts: SparseMatrix,
    /// Query file names (indexed by Q_SID)
    pub query_names: Vec<String>,
    /// Query source metadata
    pub query_sources: Vec<QuerySource>,
    /// Database source names (D_SID → name)
    pub db_sources: HashMap<u32, String>,
}

/// Execute a query against the database.
///
/// Accepts a single BED file or a directory of BED files.  Each file becomes
/// one row (Q_SID) in the result matrix.
///
/// ### Pipeline
///
/// 1. Load `meta.json`.
/// 2. Parse query files; tag each file's intervals with its Q_SID.
/// 3. For each shard present in both the database and the query:
///    a. Sort the shard's query intervals by `start`.
///    b. For each layer 0..`meta.num_layers`:
///       - Open the layer file; skip if absent.
///       - Divide query intervals into blocks (one block per available thread).
///       - `par_iter` over blocks → each runs `query_sweep` → returns a local
///         `DenseMatrix`.
///       - Tree-reduce (`reduce_with`) the local matrices.
///    c. Add shard + layer result into a cross-shard accumulator.
/// 4. Convert the final dense accumulator to sparse.
pub fn query_database(config: &QueryConfig) -> Result<QueryResult, QueryError> {
    let meta_path = config.db_path.join("meta.json");
    if !meta_path.exists() {
        return Err(QueryError::DatabaseNotFound(config.db_path.clone()));
    }

    let meta = Meta::load(&config.db_path)?;
    let num_sources = meta.num_sources();
    let db_sources: HashMap<u32, String> = meta
        .sid_map
        .iter()
        .map(|(k, v)| (*k, v.name.clone()))
        .collect();

    // Parse query files
    let parsed = parse_query_path(&config.query_path)?;

    if parsed.total_count == 0 {
        let counts = sprs::CsMat::empty(sprs::CompressedStorage::CSR, num_sources);
        return Ok(QueryResult {
            counts,
            query_names: Vec::new(),
            query_sources: Vec::new(),
            db_sources,
        });
    }

    let num_queries = parsed.total_count;

    // Determine thread count for block sizing.
    let num_threads = config
        .num_threads
        .unwrap_or_else(rayon::current_num_threads);
    let num_threads = num_threads.max(1);

    // Cross-shard dense accumulator — never more than one live at a time.
    let mut accumulator = DenseMatrix::new(num_queries, num_sources);

    let db_shard_set: std::collections::HashSet<&str> =
        meta.shards.iter().map(|s| s.as_str()).collect();

    for (shard, mut shard_queries) in parsed.shard_intervals {
        if !db_shard_set.contains(shard.as_str()) {
            continue;
        }
        if shard_queries.is_empty() {
            continue;
        }

        // Sort query intervals by start for this shard.
        shard_queries.voracious_sort();

        for layer_idx in 0..meta.num_layers {
            let layer_path = config
                .db_path
                .join(&shard)
                .join(format!("layer_{}.bin", layer_idx));

            if !layer_path.exists() {
                continue;
            }

            let layer = MappedLayer::open(&layer_path)?;
            if layer.is_empty() {
                continue;
            }

            let layer_max_size = meta.layer_config.layer_max_size(layer_idx);
            let db_intervals = layer.intervals();

            // Divide query intervals into evenly-sized blocks.
            let block_size = (shard_queries.len() + num_threads - 1) / num_threads;
            let blocks: Vec<&[Interval]> = shard_queries.chunks(block_size).collect();

            // Parallel sweep over blocks; tree-reduce results.
            let layer_result = blocks
                .par_iter()
                .map(|block| {
                    run_sweep_block(db_intervals, layer_max_size, block, num_queries, num_sources)
                })
                .reduce_with(|mut a, b| {
                    a.add_dense(&b);
                    a
                });

            if let Some(layer_dense) = layer_result {
                accumulator.add_dense(&layer_dense);
            }
        }
    }

    let final_counts = condense_to_sparse_no_mask(&accumulator, num_queries, num_sources);

    Ok(QueryResult {
        counts: final_counts,
        query_names: parsed.query_names,
        query_sources: parsed.query_sources,
        db_sources,
    })
}

/// Run a single sweep block using the thread-local scratch buffer.
fn run_sweep_block(
    db_intervals: &[Interval],
    layer_max_size: u32,
    query_block: &[Interval],
    num_queries: usize,
    num_sources: usize,
) -> DenseMatrix {
    THREAD_LOCAL_BUFFER.with(|cell| {
        let mut borrow = cell.borrow_mut();
        let scratch = borrow.get_or_insert_with(|| DenseMatrix::new(num_queries, num_sources));

        if scratch.num_rows() != num_queries || scratch.num_cols() != num_sources {
            scratch.resize_and_zero(num_queries, num_sources);
        }

        query_sweep(db_intervals, layer_max_size, query_block, scratch);

        let result = scratch.clone();
        scratch.resize_and_zero(num_queries, num_sources);
        result
    })
}

// ---------------------------------------------------------------------------
// Internal: parse query path into shards
// ---------------------------------------------------------------------------

struct ParsedQueries {
    shard_intervals: HashMap<String, Vec<Interval>>,
    query_sources: Vec<QuerySource>,
    query_names: Vec<String>,
    total_count: usize,
}

fn parse_query_path(path: &Path) -> Result<ParsedQueries, QueryError> {
    let bed_files: Vec<PathBuf> = if path.is_dir() {
        fs::read_dir(path)?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| is_bed_file(p))
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

    let mut shard_intervals: HashMap<String, Vec<Interval>> = HashMap::new();
    let mut query_sources = Vec::new();
    let mut query_names = Vec::new();
    let mut file_sid = 0u32;

    for bed_path in bed_files {
        let file_shards = parse_bed_file(&bed_path, file_sid)?;
        let file_count: usize = file_shards.values().map(|v| v.len()).sum();

        if file_count == 0 {
            continue;
        }

        let name = bed_path
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| format!("query_{}", file_sid));

        query_names.push(name.clone());
        query_sources.push(QuerySource { name, count: file_count });

        for (shard, intervals) in file_shards {
            shard_intervals.entry(shard).or_default().extend(intervals);
        }

        file_sid += 1;
    }

    Ok(ParsedQueries {
        shard_intervals,
        query_sources,
        query_names,
        total_count: file_sid as usize,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tasks::build::{AddConfig, add_to_database};
    use std::io::Write;
    use tempfile::TempDir;

    fn write_bed(dir: &Path, name: &str, content: &str) {
        let path = dir.join(name);
        let mut f = fs::File::create(path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
    }

    #[test]
    fn test_query_basic() {
        let input = TempDir::new().unwrap();
        let db = TempDir::new().unwrap();
        let query_dir = TempDir::new().unwrap();

        write_bed(input.path(), "a.bed", "chr1\t100\t200\nchr1\t300\t400\n");
        write_bed(input.path(), "b.bed", "chr1\t150\t250\n");
        add_to_database(&AddConfig::new(input.path().to_path_buf(), db.path().to_path_buf())).unwrap();

        write_bed(query_dir.path(), "q.bed", "chr1\t100\t200\n");
        let config = QueryConfig::new(db.path().to_path_buf(), query_dir.path().join("q.bed"));
        let result = query_database(&config).unwrap();

        assert_eq!(result.counts.rows(), 1);
        assert_eq!(result.db_sources.len(), 2);
    }

    #[test]
    fn test_query_multi_shard() {
        let input = TempDir::new().unwrap();
        let db = TempDir::new().unwrap();

        write_bed(input.path(), "a.bed", "chr1\t100\t200\nchr2\t300\t400\n");
        add_to_database(&AddConfig::new(input.path().to_path_buf(), db.path().to_path_buf())).unwrap();

        let query_dir = TempDir::new().unwrap();
        write_bed(query_dir.path(), "q.bed", "chr1\t100\t200\nchr2\t300\t400\n");
        let config = QueryConfig::new(db.path().to_path_buf(), query_dir.path().join("q.bed"));
        let result = query_database(&config).unwrap();

        assert_eq!(result.counts.rows(), 1);
    }

    #[test]
    fn test_query_shard_not_in_db() {
        let input = TempDir::new().unwrap();
        let db = TempDir::new().unwrap();

        write_bed(input.path(), "a.bed", "chr1\t100\t200\n");
        add_to_database(&AddConfig::new(input.path().to_path_buf(), db.path().to_path_buf())).unwrap();

        let query_dir = TempDir::new().unwrap();
        write_bed(query_dir.path(), "q.bed", "chr2\t100\t200\n"); // chr2 not in db
        let config = QueryConfig::new(db.path().to_path_buf(), query_dir.path().join("q.bed"));
        let result = query_database(&config).unwrap();

        assert_eq!(result.counts.rows(), 1);
        // No overlaps since shard is absent
        assert_eq!(result.counts.nnz(), 0);
    }

    #[test]
    fn test_query_db_not_found() {
        let config = QueryConfig::new(
            PathBuf::from("/nonexistent/path"),
            PathBuf::from("/tmp/q.bed"),
        );
        assert!(matches!(query_database(&config), Err(QueryError::DatabaseNotFound(_))));
    }

    #[test]
    fn test_query_batch_directory() {
        let input = TempDir::new().unwrap();
        let db = TempDir::new().unwrap();

        write_bed(input.path(), "a.bed", "chr1\t100\t200\nchr1\t300\t400\n");
        add_to_database(&AddConfig::new(input.path().to_path_buf(), db.path().to_path_buf())).unwrap();

        let query_dir = TempDir::new().unwrap();
        write_bed(query_dir.path(), "q1.bed", "chr1\t100\t200\n");
        write_bed(query_dir.path(), "q2.bed", "chr1\t300\t400\nchr1\t350\t450\n");

        let config = QueryConfig::new(db.path().to_path_buf(), query_dir.path().to_path_buf());
        let result = query_database(&config).unwrap();

        assert_eq!(result.counts.rows(), 2);
        assert_eq!(result.query_sources.len(), 2);
        let total_ivs: usize = result.query_sources.iter().map(|s| s.count).sum();
        assert_eq!(total_ivs, 3);
    }

    #[test]
    fn test_dense_fold_and_sparse_conversion() {
        let matrices: Vec<DenseMatrix> = (1u32..=5)
            .map(|i| {
                let mut d = DenseMatrix::new(2, 3);
                d.set(0, 0, i);
                d
            })
            .collect();

        let folded = matrices
            .into_iter()
            .fold(DenseMatrix::new(2, 3), |mut acc, m| {
                acc.add_dense(&m);
                acc
            });

        let sparse = condense_to_sparse_no_mask(&folded, 2, 3);
        assert_eq!(sparse.get(0, 0), Some(&15)); // 1+2+3+4+5
        assert_eq!(sparse.get(1, 0), None);
    }
}
