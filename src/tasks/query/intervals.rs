use std::collections::HashSet;

use rayon::prelude::*;
use voracious_radix_sort::RadixSort;

use crate::io::{MappedJumpTable, MappedLayer, Meta};
use crate::stats::IntervalHit;
use crate::sweep::query_sweep_pairs;

use super::files::{enumerate_query_files, parse_file_batch};
use super::{QueryConfig, QueryError};

/// Run the interval-pair sweep: return every (query, DB) interval overlap across
/// all query files and all shards.  No count matrix is built; the sweep result
/// is a flat list of [`IntervalHit`]s ready for serialization.
pub(super) fn collect_interval_pairs(
    config: &QueryConfig,
    meta: &Meta,
    db_sources: &std::collections::HashMap<u32, String>,
) -> Result<Vec<IntervalHit>, QueryError> {
    let all_files = enumerate_query_files(&config.query_path)?;
    if all_files.is_empty() {
        return Ok(Vec::new());
    }

    let num_threads = config
        .num_threads
        .unwrap_or_else(rayon::current_num_threads)
        .max(1);
    let batch_size = config.batch_size.unwrap_or(all_files.len()).max(1);
    let db_shard_set: HashSet<&str> = meta.shards.iter().map(|s| s.as_str()).collect();

    let mut all_hits: Vec<IntervalHit> = Vec::new();

    for batch_files in all_files.chunks(batch_size) {
        let parsed = parse_file_batch(batch_files)?;
        if parsed.total_count == 0 {
            continue;
        }

        for (shard, mut shard_queries) in parsed.shard_intervals {
            if !db_shard_set.contains(shard.as_str()) || shard_queries.is_empty() {
                continue;
            }
            shard_queries.voracious_sort();

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
                let db_positions = layer.positions();
                let db_sids = layer.sids();

                let block_size = (shard_queries.len() + num_threads - 1) / num_threads;
                let blocks: Vec<&[_]> = shard_queries.chunks(block_size).collect();

                let pairs: Vec<(_, _, u32)> = blocks
                    .par_iter()
                    .flat_map(|block| {
                        query_sweep_pairs(db_positions, db_sids, layer_max_size, &jump_table, block)
                    })
                    .collect();

                for (q, d, intersection_bp) in pairs {
                    let query_name = parsed
                        .query_names
                        .get(q.sid as usize)
                        .cloned()
                        .unwrap_or_default();
                    let db_name = db_sources
                        .get(&d.sid)
                        .cloned()
                        .unwrap_or_default();
                    all_hits.push(IntervalHit {
                        query_name,
                        db_name,
                        chrom: shard.clone(),
                        query_start: q.start,
                        query_end: q.end,
                        db_start: d.start,
                        db_end: d.end,
                        intersection_bp,
                    });
                }
            }
        }
    }

    Ok(all_hits)
}
