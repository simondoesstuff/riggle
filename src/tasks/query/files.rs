use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use rayon::prelude::*;

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
    // Parse all files in parallel; use position index as a temporary SID so
    // that intervals from each file are tagged distinctly before we know which
    // files are non-empty.
    type RawEntry = (String, HashMap<String, Vec<Interval>>, Vec<QueryChromData>, usize);
    let raw: Vec<Result<RawEntry, QueryError>> = files
        .par_iter()
        .enumerate()
        .map(|(i, bed_path)| {
            let file_shards = parse_bed_file(bed_path, i as u32)?;
            let file_count: usize = file_shards.values().map(|v| v.len()).sum();
            let name = bed_path
                .file_name()
                .map(|n| n.to_string_lossy().to_string())
                .unwrap_or_else(|| format!("query_{}", i));
            let chrom_data = if file_count > 0 {
                let bed_map: BedMap = file_shards
                    .iter()
                    .map(|(chrom, ivs)| {
                        (chrom.clone(), ivs.iter().map(|iv| (iv.start, iv.end)).collect())
                    })
                    .collect();
                build_query_chrom_data(&bed_map)
            } else {
                Vec::new()
            };
            Ok((name, file_shards, chrom_data, file_count))
        })
        .collect();

    let mut shard_intervals: HashMap<String, Vec<Interval>> = HashMap::new();
    let mut query_sources = Vec::new();
    let mut query_names = Vec::new();
    let mut query_chrom_data = Vec::new();
    let mut file_sid = 0u32;

    for (position, result) in raw.into_iter().enumerate() {
        let (name, file_shards, chrom_data, file_count) = result?;

        if file_count == 0 {
            continue;
        }

        query_names.push(name.clone());
        query_sources.push(QuerySource { name, count: file_count });
        query_chrom_data.push(chrom_data);

        for (shard, intervals) in file_shards {
            // Remap SID from position index to compact Q_SID.  For datasets
            // with no empty files, position == file_sid so this is a no-op.
            let remapped: Vec<Interval> = if position as u32 == file_sid {
                intervals
            } else {
                intervals
                    .into_iter()
                    .map(|iv| Interval { sid: file_sid, ..iv })
                    .collect()
            };
            shard_intervals.entry(shard).or_default().extend(remapped);
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

