mod files;
mod intervals;
mod pvalues;
mod sweep;

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use rayon::prelude::*;
use thiserror::Error;

use crate::fourier::QueryChromData;
use crate::io::{BedParseError, LayerError, MappedJumpTable, MappedLayer, Meta, MetaError, Position};
use crate::matrix::{DenseMatrix, OverlapMatrix, SparseMatrix};
use crate::stats::IntervalHit;
use crate::sweep::query_sweep_windowed;

use files::{enumerate_query_files, parse_file_batch};
use intervals::collect_interval_pairs;
use pvalues::compute_analytic_pvalues;
use sweep::build_csr_from_sorted_entries;

const SWEEP_WINDOW_SIZE: usize = 64;

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

/// Output mode for a query operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum QueryMode {
    /// Return file-level overlap counts only (default).
    #[default]
    Counts,
    /// Compute analytic NB/Poisson p-values and LLR for every (query, DB) pair.
    /// Requires the database to have been built with `--stats`.
    Stats,
    /// Return the exact set of overlapping interval pairs with intersection lengths.
    /// Mutually exclusive with `Stats`.
    Intervals,
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
    pub mode: QueryMode,
}

impl QueryConfig {
    pub fn new(db_path: PathBuf, query_path: PathBuf) -> Self {
        Self {
            db_path,
            query_path,
            num_threads: None,
            batch_size: None,
            mode: QueryMode::default(),
        }
    }
}

/// Per-file metadata returned with results
#[derive(Debug, Clone)]
pub struct QuerySource {
    pub name: String,
    pub count: usize,
}

/// Analytic NB statistics for one (query, DB-source) pair.
#[derive(Debug, Clone)]
pub struct PValueResult {
    /// Row index in the overlap matrix (= Q_SID).
    pub query_id: usize,
    /// Database source identifier.
    pub db_sid: u32,
    /// Base-pair overlap measured in 100 bp bins (dot product of coverage arrays).
    pub observed_bins: f64,
    /// Right-tailed p-value under the NB null (Wilks approximation).
    pub p_value: f64,
    /// Saddlepoint LLR under the NB null.  None when the NB fit is not feasible.
    pub llr: Option<f64>,
}

/// Output of a query operation
#[derive(Debug)]
pub struct QueryResult {
    /// Sparse intersection count matrix: rows = Q_SIDs, cols = D_SIDs.
    /// Empty when `mode == QueryMode::Intervals`.
    pub counts: SparseMatrix,
    /// Query file names (indexed by Q_SID)
    pub query_names: Vec<String>,
    /// Query source metadata
    pub query_sources: Vec<QuerySource>,
    /// Database source names (D_SID → name)
    pub db_sources: HashMap<u32, String>,
    /// Analytic p-values; populated only when `mode == QueryMode::Stats`.
    pub pvalues: Vec<PValueResult>,
    /// Individual interval pairs; populated only when `mode == QueryMode::Intervals`.
    pub interval_hits: Vec<IntervalHit>,
}

/// Execute a query against the database.
///
/// Accepts a single BED file or a directory of BED files.  Each file becomes
/// one row (Q_SID) in the result matrix.
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

    // Intervals mode is a wholly separate sweep path: no count matrix is built.
    if config.mode == QueryMode::Intervals {
        let hits = collect_interval_pairs(config, &meta, &db_sources)?;
        let empty_counts = sprs::CsMat::empty(sprs::CompressedStorage::CSR, num_sources);
        return Ok(QueryResult {
            counts: empty_counts,
            query_names: Vec::new(),
            query_sources: Vec::new(),
            db_sources,
            pvalues: Vec::new(),
            interval_hits: hits,
        });
    }

    let all_files = enumerate_query_files(&config.query_path)?;

    if all_files.is_empty() {
        let counts = sprs::CsMat::empty(sprs::CompressedStorage::CSR, num_sources);
        return Ok(QueryResult {
            counts,
            query_names: Vec::new(),
            query_sources: Vec::new(),
            db_sources,
            pvalues: Vec::new(),
            interval_hits: Vec::new(),
        });
    }

    let _num_threads = config
        .num_threads
        .unwrap_or_else(rayon::current_num_threads)
        .max(1);

    let batch_size = config.batch_size.unwrap_or(all_files.len()).max(1);

    let db_shard_set: HashSet<&str> = meta.shards.iter().map(|s| s.as_str()).collect();

    let mut all_entries: Vec<(usize, usize, u32)> = Vec::new();
    let mut query_names: Vec<String> = Vec::new();
    let mut query_sources: Vec<QuerySource> = Vec::new();
    let mut all_query_chrom_data: Vec<Vec<QueryChromData>> = Vec::new();
    let mut global_q_offset = 0usize;
    // Exact overlap in bins for each (q_sid, d_sid) pair; populated for every hit.
    let mut overlap_map: HashMap<(usize, usize), f32> = HashMap::new();

    for batch_files in all_files.chunks(batch_size) {
        let parsed = parse_file_batch(batch_files)?;
        let batch_len = parsed.total_count;
        if batch_len == 0 {
            continue;
        }

        let mut count_accumulator = DenseMatrix::new(batch_len, num_sources);
        let mut overlap_accumulator = OverlapMatrix::new(batch_len, num_sources);


        for (shard, mut shard_queries) in parsed.shard_intervals {
            if !db_shard_set.contains(shard.as_str()) {
                continue;
            }
            if shard_queries.is_empty() {
                continue;
            }

            // Sort by (q_sid, start) so windows are contiguous slices.
            shard_queries.sort_unstable_by_key(|iv| (iv.sid, iv.start));

            // Pre-load all layers for this shard before entering the parallel section.
            let mut layers: Vec<(MappedLayer, u32, MappedJumpTable)> = Vec::new();
            for layer_idx in 0..meta.num_layers {
                let pos_path = config
                    .db_path
                    .join("shards")
                    .join(&shard)
                    .join(format!("layer_{}.pos", layer_idx));
                if !pos_path.exists() {
                    continue;
                }
                let sid_path = config
                    .db_path
                    .join("shards")
                    .join(&shard)
                    .join(format!("layer_{}.sid", layer_idx));
                let layer = MappedLayer::open(&pos_path, &sid_path)?;
                if layer.is_empty() {
                    continue;
                }
                let layer_max_size = meta.layer_config.layer_max_size(layer_idx);
                let tile_size = meta.layer_config.tile_size(layer_idx);
                let idx_path = config
                    .db_path
                    .join("shards")
                    .join(&shard)
                    .join(format!("layer_{}.idx", layer_idx));
                let jump_table = MappedJumpTable::open(&idx_path, tile_size)?;
                layers.push((layer, layer_max_size, jump_table));
            }
            if layers.is_empty() {
                continue;
            }

            // Extract position/sid slices so the par_iter closure captures only
            // &[Position], &[u32], and &MappedJumpTable (all Send+Sync).
            let layer_refs: Vec<(&[Position], &[u32], u32, &MappedJumpTable)> =
                layers.iter().map(|(l, ms, jt)| (l.positions(), l.sids(), *ms, jt)).collect();

            let n_windows = (batch_len + SWEEP_WINDOW_SIZE - 1) / SWEEP_WINDOW_SIZE;

            // Parallel sweep: one task per Q_SID window.
            // Each window operates on a sub-matrix of `window_rows × num_sources`
            // (typically 64 × 1905 = ~0.5 MB) which fits in per-core cache,
            // eliminating the DRAM misses of the full 1905×1905 matrix.
            let window_results: Vec<(DenseMatrix, OverlapMatrix, usize)> = (0..n_windows)
                .into_par_iter()
                .map(|w| {
                    let q_sid_lo = (w * SWEEP_WINDOW_SIZE) as u32;
                    let q_sid_hi = ((w + 1) * SWEEP_WINDOW_SIZE).min(batch_len) as u32;
                    let window_rows = (q_sid_hi - q_sid_lo) as usize;

                    let lo = shard_queries.partition_point(|iv| iv.sid < q_sid_lo);
                    let hi = shard_queries.partition_point(|iv| iv.sid < q_sid_hi);
                    let window_queries = &shard_queries[lo..hi];

                    let mut sub_counts = DenseMatrix::new(window_rows, num_sources);
                    let mut sub_overlap = OverlapMatrix::new(window_rows, num_sources);

                    for &(db_positions, db_sids, layer_max_size, jump_table) in &layer_refs {
                        query_sweep_windowed(
                            db_positions,
                            db_sids,
                            layer_max_size,
                            jump_table,
                            window_queries,
                            q_sid_lo,
                            &mut sub_counts,
                            &mut sub_overlap,
                        );
                    }

                    (sub_counts, sub_overlap, w)
                })
                .collect();

            // Serial merge: each window's sub-matrix lands in the right rows.
            for (sub_counts, sub_overlap, w) in window_results {
                let row_offset = w * SWEEP_WINDOW_SIZE;
                count_accumulator.add_submatrix_rows(&sub_counts, row_offset);
                overlap_accumulator.add_submatrix_rows(&sub_overlap, row_offset);
            }
        }

        for local_row in 0..batch_len {
            let row_slice = count_accumulator.row(local_row);
            let global_row = global_q_offset + local_row;
            for (col, &val) in row_slice.iter().enumerate() {
                if val > 0 {
                    all_entries.push((global_row, col, val));
                    let ov = overlap_accumulator.get(local_row, col);
                    *overlap_map.entry((global_row, col)).or_default() += ov;
                }
            }
        }

        global_q_offset += batch_len;
        query_names.extend(parsed.query_names);
        query_sources.extend(parsed.query_sources);
        all_query_chrom_data.extend(parsed.query_chrom_data);
    }

    let num_queries = global_q_offset;
    let final_counts = build_csr_from_sorted_entries(&all_entries, num_queries, num_sources);

    let pvalues = if config.mode == QueryMode::Stats {
        compute_analytic_pvalues(
            &final_counts,
            &all_query_chrom_data,
            &config.db_path,
            &overlap_map,
        )
    } else {
        Vec::new()
    };

    Ok(QueryResult {
        counts: final_counts,
        query_names,
        query_sources,
        db_sources,
        pvalues,
        interval_hits: Vec::new(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::condense_to_sparse_no_mask;
    use crate::tasks::build::{AddConfig, add_to_database};
    use std::fs;
    use std::io::Write;
    use std::path::Path;
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
        add_to_database(&AddConfig::new(
            input.path().to_path_buf(),
            db.path().to_path_buf(),
        ))
        .unwrap();

        write_bed(query_dir.path(), "q.bed", "chr1\t100\t200\n");
        let config = QueryConfig::new(db.path().to_path_buf(), query_dir.path().join("q.bed"));
        let result = query_database(&config).unwrap();

        assert_eq!(result.counts.rows(), 1);
        assert_eq!(result.db_sources.len(), 2);
        assert!(result.pvalues.is_empty());
    }

    #[test]
    fn test_query_multi_shard() {
        let input = TempDir::new().unwrap();
        let db = TempDir::new().unwrap();

        write_bed(input.path(), "a.bed", "chr1\t100\t200\nchr2\t300\t400\n");
        add_to_database(&AddConfig::new(
            input.path().to_path_buf(),
            db.path().to_path_buf(),
        ))
        .unwrap();

        let query_dir = TempDir::new().unwrap();
        write_bed(
            query_dir.path(),
            "q.bed",
            "chr1\t100\t200\nchr2\t300\t400\n",
        );
        let config = QueryConfig::new(db.path().to_path_buf(), query_dir.path().join("q.bed"));
        let result = query_database(&config).unwrap();

        assert_eq!(result.counts.rows(), 1);
    }

    #[test]
    fn test_query_shard_not_in_db() {
        let input = TempDir::new().unwrap();
        let db = TempDir::new().unwrap();

        write_bed(input.path(), "a.bed", "chr1\t100\t200\n");
        add_to_database(&AddConfig::new(
            input.path().to_path_buf(),
            db.path().to_path_buf(),
        ))
        .unwrap();

        let query_dir = TempDir::new().unwrap();
        write_bed(query_dir.path(), "q.bed", "chr2\t100\t200\n");
        let config = QueryConfig::new(db.path().to_path_buf(), query_dir.path().join("q.bed"));
        let result = query_database(&config).unwrap();

        assert_eq!(result.counts.rows(), 1);
        assert_eq!(result.counts.nnz(), 0);
    }

    #[test]
    fn test_query_db_not_found() {
        let config = QueryConfig::new(
            PathBuf::from("/nonexistent/path"),
            PathBuf::from("/tmp/q.bed"),
        );
        assert!(matches!(
            query_database(&config),
            Err(QueryError::DatabaseNotFound(_))
        ));
    }

    #[test]
    fn test_query_batch_directory() {
        let input = TempDir::new().unwrap();
        let db = TempDir::new().unwrap();

        write_bed(input.path(), "a.bed", "chr1\t100\t200\nchr1\t300\t400\n");
        add_to_database(&AddConfig::new(
            input.path().to_path_buf(),
            db.path().to_path_buf(),
        ))
        .unwrap();

        let query_dir = TempDir::new().unwrap();
        write_bed(query_dir.path(), "q1.bed", "chr1\t100\t200\n");
        write_bed(
            query_dir.path(),
            "q2.bed",
            "chr1\t300\t400\nchr1\t350\t450\n",
        );

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

        write_bed(input.path(), "a.bed", "chr1\t100\t200\nchr1\t300\t400\n");
        write_bed(input.path(), "b.bed", "chr1\t150\t250\n");
        write_bed(input.path(), "c.bed", "chr1\t350\t500\n");
        add_to_database(&AddConfig::new(
            input.path().to_path_buf(),
            db.path().to_path_buf(),
        ))
        .unwrap();

        let query_dir = TempDir::new().unwrap();
        write_bed(query_dir.path(), "q1.bed", "chr1\t100\t200\n");
        write_bed(query_dir.path(), "q2.bed", "chr1\t200\t350\n");
        write_bed(query_dir.path(), "q3.bed", "chr1\t300\t450\n");
        write_bed(query_dir.path(), "q4.bed", "chr1\t50\t600\n");

        let unbatched = {
            let config = QueryConfig::new(db.path().to_path_buf(), query_dir.path().to_path_buf());
            query_database(&config).unwrap()
        };

        let batched = {
            let mut config =
                QueryConfig::new(db.path().to_path_buf(), query_dir.path().to_path_buf());
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
        ))
        .unwrap();

        let mut batched_config =
            AddConfig::new(input.path().to_path_buf(), db_batched.path().to_path_buf());
        batched_config.batch_size = Some(1);
        add_to_database(&batched_config).unwrap();

        let query_dir = TempDir::new().unwrap();
        write_bed(query_dir.path(), "q.bed", "chr1\t100\t500\n");

        let result_u = query_database(&QueryConfig::new(
            db_unbatched.path().to_path_buf(),
            query_dir.path().join("q.bed"),
        ))
        .unwrap();

        let result_b = query_database(&QueryConfig::new(
            db_batched.path().to_path_buf(),
            query_dir.path().join("q.bed"),
        ))
        .unwrap();

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
        assert_eq!(sparse.get(0, 0), Some(&15));
        assert_eq!(sparse.get(1, 0), None);
    }

    #[test]
    fn test_stats_only_overlapping_pairs_emitted() {
        let input = TempDir::new().unwrap();
        let db = TempDir::new().unwrap();

        // Use chr22 (510K bins) to keep the FFT moment computation fast in tests.
        write_bed(input.path(), "a.bed", "chr22\t100\t200\n");
        write_bed(input.path(), "b.bed", "chr22\t5000000\t5000100\n");
        let mut add_cfg = AddConfig::new(input.path().to_path_buf(), db.path().to_path_buf());
        add_cfg.compute_stats = true;
        add_to_database(&add_cfg).unwrap();

        let query_dir = TempDir::new().unwrap();
        write_bed(query_dir.path(), "q.bed", "chr22\t100\t200\n");

        let mut config =
            QueryConfig::new(db.path().to_path_buf(), query_dir.path().join("q.bed"));
        config.mode = QueryMode::Stats;
        let result = query_database(&config).unwrap();

        // Only the overlapping pair (q vs a) should appear; zero-overlap pairs are absent.
        assert_eq!(result.pvalues.len(), 1);
        assert!(result.pvalues[0].observed_bins > 0.0);
        // The non-overlapping pair (q vs b) is not emitted.
        assert!(!result.pvalues.iter().any(|pv| pv.observed_bins == 0.0));
    }
}
