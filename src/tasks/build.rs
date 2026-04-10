use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

use rayon::prelude::*;
use thiserror::Error;
use voracious_radix_sort::RadixSort;

use crate::core::Interval;
use crate::io::{
    BedParseError, LayerConfig, LayerError, Meta, MetaError, SidEntry, extend_layer, is_bed_file,
    parse_bed_file, write_layer,
};

/// Errors from the add / build pipeline
#[derive(Debug, Error)]
pub enum AddError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("BED parse error: {0}")]
    BedParse(#[from] BedParseError),
    #[error("Meta error: {0}")]
    Meta(#[from] MetaError),
    #[error("Layer error: {0}")]
    Layer(#[from] LayerError),
    #[error("Input path not found: {0}")]
    InputNotFound(PathBuf),
    #[error("No BED files found in {0}")]
    NoBedFiles(PathBuf),
}

/// Configuration for `add_to_database`.
#[derive(Debug, Clone)]
pub struct AddConfig {
    /// Directory containing BED files to ingest
    pub input_path: PathBuf,
    /// Database directory (created if it does not exist)
    pub db_path: PathBuf,
    /// Layer configuration — used only when creating a new database.
    /// If `None` and the database is new, defaults to `LayerConfig::default_genomic()`.
    pub layer_config: Option<LayerConfig>,
    /// Maximum number of BED files to parse and hold in memory at once.
    /// `None` (default) processes all files in a single batch.
    pub batch_size: Option<usize>,
}

impl AddConfig {
    pub fn new(input_path: PathBuf, db_path: PathBuf) -> Self {
        Self {
            input_path,
            db_path,
            layer_config: None,
            batch_size: None,
        }
    }
}

/// Ingest BED files into the database.
///
/// If the database does not yet exist it is created.  If it already exists the
/// new intervals are merged into the existing layer files.
///
/// When `config.batch_size` is set, files are processed in sequential batches
/// of that size, bounding peak memory to roughly `batch_size` files at once.
///
/// ### Pipeline (per batch)
///
/// 1. Load or create `meta.json`.
/// 2. Collect BED files from `config.input_path`; chunk by `batch_size`.
/// 3. Parse the batch in parallel, tagging each file with a new D_SID.
/// 4. Partition intervals into `(shard, layer_idx)` buckets.
/// 5. Sort each bucket by `start` using voracious radix sort.
/// 6. Write new buckets / merge into existing layer files (parallel per bucket).
/// 7. Update and save `meta.json`.
pub fn add_to_database(config: &AddConfig) -> Result<(), AddError> {
    if !config.input_path.exists() {
        return Err(AddError::InputNotFound(config.input_path.clone()));
    }

    // Collect BED files
    let bed_files: Vec<PathBuf> = if config.input_path.is_dir() {
        fs::read_dir(&config.input_path)?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| is_bed_file(p))
            .collect()
    } else if is_bed_file(&config.input_path) {
        vec![config.input_path.clone()]
    } else {
        return Err(AddError::NoBedFiles(config.input_path.clone()));
    };

    if bed_files.is_empty() {
        return Err(AddError::NoBedFiles(config.input_path.clone()));
    }

    // Load or create database
    fs::create_dir_all(&config.db_path)?;
    let meta_path = config.db_path.join("meta.json");
    let mut meta = if meta_path.exists() {
        Meta::load(&config.db_path)?
    } else {
        let layer_config = config
            .layer_config
            .clone()
            .unwrap_or_else(LayerConfig::default_genomic);
        Meta::new(layer_config)
    };

    let batch_size = config.batch_size.unwrap_or(bed_files.len()).max(1);
    for batch in bed_files.chunks(batch_size) {
        process_file_batch(batch, &mut meta, &config.db_path)?;
        meta.save(&config.db_path)?;
    }

    Ok(())
}

/// Process one batch of BED files: parse → bucket → sort → write/merge → update meta.
fn process_file_batch(
    batch: &[PathBuf],
    meta: &mut Meta,
    db_path: &std::path::Path,
) -> Result<(), AddError> {
    let layer_config = meta.layer_config.clone();
    let next_sid = meta.next_sid();

    // Parse all files in parallel.
    // Each file is assigned a unique D_SID starting from `next_sid`.
    let parse_results: Vec<Result<(u32, String, HashMap<String, Vec<Interval>>), AddError>> =
        batch
            .par_iter()
            .enumerate()
            .map(|(i, path)| {
                let sid = next_sid + i as u32;
                let name = path
                    .file_name()
                    .map(|n| n.to_string_lossy().to_string())
                    .unwrap_or_else(|| format!("file_{}", sid));
                let shards = parse_bed_file(path, sid)?;
                Ok((sid, name, shards))
            })
            .collect();

    // Partition all intervals into (shard, layer_idx) buckets.
    let mut buckets: HashMap<(String, usize), Vec<Interval>> = HashMap::new();
    let mut new_sids: Vec<(u32, String)> = Vec::new();

    for result in parse_results {
        let (sid, name, shard_map) = result?;
        new_sids.push((sid, name));
        for (shard, intervals) in shard_map {
            for iv in intervals {
                let layer = layer_config.layer_for_size(iv.size());
                buckets
                    .entry((shard.clone(), layer))
                    .or_default()
                    .push(iv);
            }
        }
    }

    // Sort each bucket by start coordinate using voracious radix sort.
    let mut bucket_vec: Vec<((String, usize), Vec<Interval>)> = buckets.into_iter().collect();
    bucket_vec.par_iter_mut().for_each(|(_, ivs)| {
        ivs.voracious_sort();
    });

    // Write / merge each bucket into the database (parallel across buckets).
    let write_errors: Vec<AddError> = bucket_vec
        .par_iter()
        .filter_map(|((shard, layer_idx), sorted_ivs)| {
            if sorted_ivs.is_empty() {
                return None;
            }
            let shard_dir = db_path.join(shard);
            if let Err(e) = fs::create_dir_all(&shard_dir) {
                return Some(AddError::Io(e));
            }
            let layer_path = shard_dir.join(format!("layer_{}.bin", layer_idx));
            let result = if layer_path.exists() {
                extend_layer(&layer_path, sorted_ivs)
            } else {
                write_layer(&layer_path, sorted_ivs)
            };
            result.err().map(AddError::Layer)
        })
        .collect();

    if let Some(err) = write_errors.into_iter().next() {
        return Err(err);
    }

    // Update meta: add new SID entries, shards, num_layers.
    let mut max_layer_used = meta.num_layers.saturating_sub(1);

    for (sid, name) in new_sids {
        meta.sid_map.insert(sid, SidEntry { name });
    }

    for (shard, layer_idx) in bucket_vec.iter().map(|((s, l), _)| (s, *l)) {
        if !meta.shards.contains(shard) {
            meta.shards.push(shard.clone());
        }
        if layer_idx > max_layer_used {
            max_layer_used = layer_idx;
        }
    }

    meta.num_layers = max_layer_used + 1;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::{MappedLayer, Meta};
    use std::io::Write;
    use tempfile::TempDir;

    fn write_bed(dir: &std::path::Path, name: &str, content: &str) {
        let path = dir.join(name);
        let mut f = fs::File::create(path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
    }

    #[test]
    fn test_add_creates_database() {
        let input = TempDir::new().unwrap();
        let db = TempDir::new().unwrap();

        write_bed(input.path(), "a.bed", "chr1\t100\t200\nchr1\t300\t400\n");

        add_to_database(&AddConfig::new(
            input.path().to_path_buf(),
            db.path().to_path_buf(),
        ))
        .unwrap();

        assert!(db.path().join("meta.json").exists());
        let meta = Meta::load(db.path()).unwrap();
        assert_eq!(meta.sid_map.len(), 1);
        assert!(meta.shards.contains(&"chr1".to_string()));
    }

    #[test]
    fn test_add_twice_accumulates_sids() {
        let input1 = TempDir::new().unwrap();
        let input2 = TempDir::new().unwrap();
        let db = TempDir::new().unwrap();

        write_bed(input1.path(), "a.bed", "chr1\t100\t200\n");
        write_bed(input2.path(), "b.bed", "chr1\t300\t400\n");

        add_to_database(&AddConfig::new(
            input1.path().to_path_buf(),
            db.path().to_path_buf(),
        ))
        .unwrap();
        add_to_database(&AddConfig::new(
            input2.path().to_path_buf(),
            db.path().to_path_buf(),
        ))
        .unwrap();

        let meta = Meta::load(db.path()).unwrap();
        assert_eq!(meta.sid_map.len(), 2);
        assert!(meta.sid_map.contains_key(&0));
        assert!(meta.sid_map.contains_key(&1));
    }

    #[test]
    fn test_add_layer_file_sorted() {
        let input = TempDir::new().unwrap();
        let db = TempDir::new().unwrap();

        // Two files with interleaved coords — merged layer must be sorted
        write_bed(input.path(), "a.bed", "chr1\t300\t400\nchr1\t100\t200\n");
        write_bed(input.path(), "b.bed", "chr1\t150\t250\nchr1\t350\t450\n");

        add_to_database(&AddConfig::new(
            input.path().to_path_buf(),
            db.path().to_path_buf(),
        ))
        .unwrap();

        let meta = Meta::load(db.path()).unwrap();
        // All these intervals have size 100, which falls in layer 0 (size < 128)
        let layer_path = db.path().join("chr1").join("layer_0.bin");
        assert!(layer_path.exists());
        let mapped = MappedLayer::open(&layer_path).unwrap();
        let ivs = mapped.intervals();
        for w in ivs.windows(2) {
            assert!(w[0].start <= w[1].start);
        }
    }

    #[test]
    fn test_add_multi_shard() {
        let input = TempDir::new().unwrap();
        let db = TempDir::new().unwrap();

        write_bed(
            input.path(),
            "a.bed",
            "chr1\t100\t200\nchr2\t300\t400\n",
        );

        add_to_database(&AddConfig::new(
            input.path().to_path_buf(),
            db.path().to_path_buf(),
        ))
        .unwrap();

        let meta = Meta::load(db.path()).unwrap();
        assert_eq!(meta.shards.len(), 2);
        assert!(db.path().join("chr1").exists());
        assert!(db.path().join("chr2").exists());
    }
}
