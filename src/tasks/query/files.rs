use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use crate::core::Interval;
use crate::fourier::{BedMap, QueryChromData, build_query_chrom_data};
use crate::io::{is_bed_file, parse_bed_file};

use super::{QueryError, QuerySource};

pub(super) struct ParsedQueries {
    pub shard_intervals: HashMap<String, Vec<Interval>>,
    pub query_sources: Vec<QuerySource>,
    pub query_names: Vec<String>,
    /// Pre-built per-chrom interval data for each non-empty file in this batch,
    /// indexed by local Q_SID.  Avoids re-parsing at stats time.
    pub query_chrom_data: Vec<Vec<QueryChromData>>,
    /// Number of non-empty files in this batch (= number of Q_SID rows used).
    pub total_count: usize,
}

/// Return all BED file paths under `path` (or `path` itself if it is a file).
pub(super) fn enumerate_query_files(path: &Path) -> Result<Vec<PathBuf>, QueryError> {
    if path.is_dir() {
        let mut files: Vec<PathBuf> = fs::read_dir(path)?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| is_bed_file(p))
            .collect();
        files.sort();
        Ok(files)
    } else {
        Ok(vec![path.to_path_buf()])
    }
}

/// Parse a slice of BED files, assigning local Q_SIDs 0..N (one per non-empty file).
pub(super) fn parse_file_batch(files: &[PathBuf]) -> Result<ParsedQueries, QueryError> {
    let mut shard_intervals: HashMap<String, Vec<Interval>> = HashMap::new();
    let mut query_sources = Vec::new();
    let mut query_names = Vec::new();
    let mut query_chrom_data = Vec::new();
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

        let bed_map: BedMap = file_shards
            .iter()
            .map(|(chrom, ivs)| {
                (chrom.clone(), ivs.iter().map(|iv| (iv.start, iv.end)).collect())
            })
            .collect();
        query_chrom_data.push(build_query_chrom_data(&bed_map));

        for (shard, intervals) in file_shards {
            shard_intervals.entry(shard).or_default().extend(intervals);
        }

        file_sid += 1;
    }

    Ok(ParsedQueries {
        shard_intervals,
        query_sources,
        query_names,
        query_chrom_data,
        total_count: file_sid as usize,
    })
}

