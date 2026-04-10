use std::cell::RefCell;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use rayon::prelude::*;
use thiserror::Error;
use voracious_radix_sort::RadixSort;

use crate::core::Interval;
use crate::io::{BedParseError, LayerError, Meta, MetaError, MappedLayer, is_bed_file, parse_bed_file};
use crate::matrix::{DenseMatrix, SparseMatrix};
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
    /// Maximum number of query files to parse and hold in memory at once.
    /// `None` (default) processes all query files in a single batch.
    pub batch_size: Option<usize>,
}

impl QueryConfig {
    pub fn new(db_path: PathBuf, query_path: PathBuf) -> Self {
        Self {
            db_path,
            query_path,
            num_threads: None,
            batch_size: None,
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
/// When `config.batch_size` is set, query files are processed in sequential
/// batches, bounding peak memory to roughly `batch_size` files at once.  The
/// per-batch dense accumulator is `batch_size × num_sources` instead of
/// `num_queries × num_sources`.
///
/// ### Pipeline
///
/// 1. Load `meta.json`; enumerate query files.
/// 2. For each batch of `batch_size` query files:
///    a. Parse the batch; assign local Q_SIDs 0..batch_len.
///    b. For each shard × layer: sweep, reduce, add to a batch-local accumulator.
///    c. Drain non-zero entries from the accumulator, shifting rows by the
///       global batch offset, into `all_entries`.
/// 3. Build the final CSR matrix from `all_entries` (already (row,col)-sorted).
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

    let all_files = enumerate_query_files(&config.query_path)?;

    if all_files.is_empty() {
        let counts = sprs::CsMat::empty(sprs::CompressedStorage::CSR, num_sources);
        return Ok(QueryResult {
            counts,
            query_names: Vec::new(),
            query_sources: Vec::new(),
            db_sources,
        });
    }

    // Determine thread count for block sizing.
    let num_threads = config
        .num_threads
        .unwrap_or_else(rayon::current_num_threads)
        .max(1);

    let batch_size = config.batch_size.unwrap_or(all_files.len()).max(1);

    let db_shard_set: std::collections::HashSet<&str> =
        meta.shards.iter().map(|s| s.as_str()).collect();

    // Accumulated (global_row, col, val) triples — sorted by (row, col).
    let mut all_entries: Vec<(usize, usize, u32)> = Vec::new();
    let mut query_names: Vec<String> = Vec::new();
    let mut query_sources: Vec<QuerySource> = Vec::new();
    let mut global_q_offset = 0usize;

    for batch_files in all_files.chunks(batch_size) {
        let parsed = parse_file_batch(batch_files)?;
        let batch_len = parsed.total_count;
        if batch_len == 0 {
            continue;
        }

        // Per-batch dense accumulator: batch_len × num_sources.
        let mut accumulator = DenseMatrix::new(batch_len, num_sources);

        for (shard, mut shard_queries) in parsed.shard_intervals {
            if !db_shard_set.contains(shard.as_str()) {
                continue;
            }
            if shard_queries.is_empty() {
                continue;
            }

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

                let block_size = (shard_queries.len() + num_threads - 1) / num_threads;
                let blocks: Vec<&[Interval]> = shard_queries.chunks(block_size).collect();

                let layer_result = blocks
                    .par_iter()
                    .map(|block| {
                        run_sweep_block(db_intervals, layer_max_size, block, batch_len, num_sources)
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

        // Drain non-zero entries into all_entries with global row indices.
        for local_row in 0..batch_len {
            let row_slice = accumulator.row(local_row);
            for (col, &val) in row_slice.iter().enumerate() {
                if val > 0 {
                    all_entries.push((global_q_offset + local_row, col, val));
                }
            }
        }

        global_q_offset += batch_len;
        query_names.extend(parsed.query_names);
        query_sources.extend(parsed.query_sources);
    }

    let num_queries = global_q_offset;
    let final_counts = build_csr_from_sorted_entries(&all_entries, num_queries, num_sources);

    Ok(QueryResult {
        counts: final_counts,
        query_names,
        query_sources,
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
// Internal helpers
// ---------------------------------------------------------------------------

struct ParsedQueries {
    shard_intervals: HashMap<String, Vec<Interval>>,
    query_sources: Vec<QuerySource>,
    query_names: Vec<String>,
    /// Number of non-empty files in this batch (= number of Q_SID rows used).
    total_count: usize,
}

/// Return all BED file paths under `path` (or `path` itself if it is a file).
fn enumerate_query_files(path: &Path) -> Result<Vec<PathBuf>, QueryError> {
    if path.is_dir() {
        let mut files: Vec<PathBuf> = fs::read_dir(path)?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| is_bed_file(p))
            .collect();
        files.sort(); // deterministic ordering
        Ok(files)
    } else {
        Ok(vec![path.to_path_buf()])
    }
}

/// Parse a slice of BED files, assigning local Q_SIDs 0..N (one per non-empty file).
fn parse_file_batch(files: &[PathBuf]) -> Result<ParsedQueries, QueryError> {
    let mut shard_intervals: HashMap<String, Vec<Interval>> = HashMap::new();
    let mut query_sources = Vec::new();
    let mut query_names = Vec::new();
    let mut file_sid = 0u32;

    for bed_path in files {
        let file_shards = parse_bed_file(bed_path, file_sid)?;
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

/// Build a CSR sparse matrix from entries that are already sorted by (row, col).
fn build_csr_from_sorted_entries(
    entries: &[(usize, usize, u32)],
    num_rows: usize,
    num_cols: usize,
) -> SparseMatrix {
    let mut indptr = vec![0usize; num_rows + 1];
    let mut indices = Vec::with_capacity(entries.len());
    let mut data = Vec::with_capacity(entries.len());

    for &(row, col, val) in entries {
        indices.push(col);
        data.push(val);
        indptr[row + 1] += 1;
    }
    // Convert per-row counts to prefix sums.
    for i in 0..num_rows {
        indptr[i + 1] += indptr[i];
    }

    sprs::CsMat::new((num_rows, num_cols), indptr, indices, data)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::condense_to_sparse_no_mask;
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
    fn test_batch_query_matches_unbatched() {
        let input = TempDir::new().unwrap();
        let db = TempDir::new().unwrap();

        // 3 DB source files
        write_bed(input.path(), "a.bed", "chr1\t100\t200\nchr1\t300\t400\n");
        write_bed(input.path(), "b.bed", "chr1\t150\t250\n");
        write_bed(input.path(), "c.bed", "chr1\t350\t500\n");
        add_to_database(&AddConfig::new(input.path().to_path_buf(), db.path().to_path_buf())).unwrap();

        // 4 query files
        let query_dir = TempDir::new().unwrap();
        write_bed(query_dir.path(), "q1.bed", "chr1\t100\t200\n");
        write_bed(query_dir.path(), "q2.bed", "chr1\t200\t350\n");
        write_bed(query_dir.path(), "q3.bed", "chr1\t300\t450\n");
        write_bed(query_dir.path(), "q4.bed", "chr1\t50\t600\n");

        let unbatched = {
            let config = QueryConfig::new(db.path().to_path_buf(), query_dir.path().to_path_buf());
            query_database(&config).unwrap()
        };

        // batch_size=1 forces one file per batch
        let batched = {
            let mut config = QueryConfig::new(db.path().to_path_buf(), query_dir.path().to_path_buf());
            config.batch_size = Some(1);
            query_database(&config).unwrap()
        };

        assert_eq!(unbatched.counts.rows(), batched.counts.rows());
        assert_eq!(unbatched.counts.cols(), batched.counts.cols());
        for r in 0..unbatched.counts.rows() {
            for c in 0..unbatched.counts.cols() {
                assert_eq!(
                    unbatched.counts.get(r, c),
                    batched.counts.get(r, c),
                    "mismatch at ({r},{c})"
                );
            }
        }
        assert_eq!(unbatched.query_names, batched.query_names);
    }

    #[test]
    fn test_batch_add_matches_unbatched() {
        let input = TempDir::new().unwrap();
        let db_unbatched = TempDir::new().unwrap();
        let db_batched = TempDir::new().unwrap();

        write_bed(input.path(), "a.bed", "chr1\t100\t200\nchr1\t300\t400\n");
        write_bed(input.path(), "b.bed", "chr1\t150\t250\n");
        write_bed(input.path(), "c.bed", "chr1\t350\t500\n");

        add_to_database(&AddConfig::new(
            input.path().to_path_buf(),
            db_unbatched.path().to_path_buf(),
        )).unwrap();

        let mut batched_config = AddConfig::new(
            input.path().to_path_buf(),
            db_batched.path().to_path_buf(),
        );
        batched_config.batch_size = Some(1);
        add_to_database(&batched_config).unwrap();

        let query_dir = TempDir::new().unwrap();
        write_bed(query_dir.path(), "q.bed", "chr1\t100\t500\n");

        let result_u = query_database(&QueryConfig::new(
            db_unbatched.path().to_path_buf(),
            query_dir.path().join("q.bed"),
        )).unwrap();

        let result_b = query_database(&QueryConfig::new(
            db_batched.path().to_path_buf(),
            query_dir.path().join("q.bed"),
        )).unwrap();

        assert_eq!(result_u.counts.rows(), result_b.counts.rows());
        for r in 0..result_u.counts.rows() {
            for c in 0..result_u.counts.cols() {
                assert_eq!(
                    result_u.counts.get(r, c),
                    result_b.counts.get(r, c),
                    "mismatch at ({r},{c})"
                );
            }
        }
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
