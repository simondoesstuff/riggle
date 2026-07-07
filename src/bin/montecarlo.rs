//! montecarlo — Monte Carlo shift permutation test for interval overlap.
//!
//! Null model: the entire query set (all chromosomes simultaneously) is shifted
//! by a random integer offset drawn uniformly from [-max_chrom_size, max_chrom_size].
//! Intervals that land below position 0 or past the chromosome end are trimmed or
//! removed.  The test statistic is the overlap count per database source.
//!
//! Usage:
//!   cargo run --release --bin montecarlo -- \
//!     --db data/rme_chuckle \
//!     --query data/gwas/trait.bed \
//!     --trials 100,1000,10000
//!
//! At each trial milestone, empirical right-tailed p-values are printed for
//! every database source.

use std::collections::HashMap;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;
use rand::Rng;
use rand::SeedableRng;
use rand::rngs::StdRng;
use tempfile::TempDir;

use riggle::fourier::{BedMap, hg38_chrom_sizes, parse_bed_as_map};
use riggle::tasks::{QueryConfig, query_database};

#[derive(Parser)]
#[command(name = "montecarlo")]
#[command(about = "Monte Carlo shift permutation p-value for interval overlap")]
struct Cli {
    /// Path to the riggle database directory
    #[arg(short, long)]
    db: PathBuf,

    /// Query BED file
    #[arg(short, long)]
    query: PathBuf,

    /// Comma-separated trial milestones at which to print p-values
    #[arg(short, long, value_delimiter = ',')]
    trials: Vec<usize>,

    /// Number of trials per batch (limits peak temp-file disk usage)
    #[arg(short, long, default_value = "100")]
    batch_size: usize,

    /// Random seed (default: from OS entropy)
    #[arg(short, long)]
    seed: Option<u64>,

    /// Whitelist BED: restrict shifted intervals to these regions (matches chuckle behaviour)
    #[arg(short = 'w', long)]
    whitelist: Option<PathBuf>,
}

/// Returns the sub-intervals of [q_start, q_end) that overlap sorted merged whitelist intervals.
fn intersect_whitelist(q_start: u32, q_end: u32, wl: &[(u32, u32)]) -> impl Iterator<Item = (u32, u32)> + '_ {
    let idx = wl.partition_point(|&(_, wl_end)| wl_end <= q_start);
    wl[idx..]
        .iter()
        .take_while(move |&&(wl_start, _)| wl_start < q_end)
        .map(move |&(wl_start, wl_end)| (wl_start.max(q_start), wl_end.min(q_end)))
}

fn write_shifted_bed<W: Write>(
    bed: &BedMap,
    chrom_sizes: &HashMap<&str, u32>,
    whitelist: Option<&BedMap>,
    shift: i64,
    writer: &mut W,
) -> bool {
    let mut wrote = false;
    let mut chroms: Vec<&String> = bed.keys().collect();
    chroms.sort();

    for chrom in chroms {
        let intervals = &bed[chrom];
        let chrom_size = match chrom_sizes.get(chrom.as_str()) {
            Some(&s) => s as i64,
            None => continue,
        };
        let wl_chrom = whitelist.and_then(|wl| wl.get(chrom.as_str()));

        for &(start, end) in intervals {
            let new_start = start as i64 + shift;
            let new_end = end as i64 + shift;
            if new_end <= 0 || new_start >= chrom_size {
                continue;
            }
            let trimmed_start = new_start.max(0) as u32;
            let trimmed_end = new_end.min(chrom_size) as u32;
            if trimmed_start >= trimmed_end {
                continue;
            }
            if let Some(wl) = wl_chrom {
                for (is, ie) in intersect_whitelist(trimmed_start, trimmed_end, wl) {
                    writeln!(writer, "{}\t{}\t{}", chrom, is, ie).unwrap();
                    wrote = true;
                }
            } else {
                writeln!(writer, "{}\t{}\t{}", chrom, trimmed_start, trimmed_end).unwrap();
                wrote = true;
            }
        }
    }
    wrote
}

fn print_pvalues(
    milestone: usize,
    total_trials: usize,
    exceed_counts: &HashMap<u32, u64>,
    observed_counts: &HashMap<u32, u32>,
    db_sources: &HashMap<u32, String>,
) {
    println!("# milestone={} trials={}", milestone, total_trials);
    println!("db_name\tobserved\tp_value");

    let mut d_sids: Vec<u32> = db_sources.keys().copied().collect();
    d_sids.sort();

    for d_sid in d_sids {
        let name = &db_sources[&d_sid];
        let obs = observed_counts.get(&d_sid).copied().unwrap_or(0);
        let p_value = if obs == 0 {
            1.0_f64
        } else {
            let exceed = exceed_counts.get(&d_sid).copied().unwrap_or(0);
            exceed as f64 / total_trials as f64
        };
        println!("{}\t{}\t{:.4e}", name, obs, p_value);
    }
    println!();
}

fn main() {
    let cli = Cli::parse();

    let mut milestones = cli.trials.clone();
    if milestones.is_empty() {
        eprintln!("Error: specify at least one trial milestone with --trials");
        std::process::exit(1);
    }
    milestones.sort_unstable();
    milestones.dedup();

    let total_target = *milestones.last().unwrap();
    let batch_size = cli.batch_size.max(1);

    // Parse query BED
    let query_bed: BedMap = parse_bed_as_map(&cli.query)
        .unwrap_or_else(|e| panic!("Failed to parse query BED: {e}"));

    let n_ivs: usize = query_bed.values().map(|v| v.len()).sum();
    eprintln!(
        "Query: {} intervals across {} chromosomes",
        n_ivs,
        query_bed.len()
    );

    let whitelist: Option<BedMap> = cli.whitelist.as_ref().map(|p| {
        parse_bed_as_map(p).unwrap_or_else(|e| panic!("Failed to parse whitelist: {e}"))
    });
    if whitelist.is_some() {
        eprintln!("Whitelist: active — shifts clipped to annotated regions");
    }

    // Chromosome sizes for shift range and trimming
    let chrom_sizes: HashMap<&str, u32> = hg38_chrom_sizes().iter().copied().collect();
    let max_shift = chrom_sizes.values().copied().max().unwrap_or(250_000_000) as i64;

    // Observed counts: one query file → one row in the result matrix
    let observed_result = {
        let config = QueryConfig::new(cli.db.clone(), cli.query.clone());
        query_database(&config).expect("Failed to query database for observed counts")
    };

    let db_sources = observed_result.db_sources.clone();

    let mut observed_counts: HashMap<u32, u32> = HashMap::new();
    if let Some(row_view) = observed_result.counts.outer_view(0) {
        for (col, &val) in row_view.iter() {
            if val > 0 {
                observed_counts.insert(col as u32, val);
            }
        }
    }

    eprintln!(
        "Database: {} sources; {} have non-zero overlap with query",
        db_sources.len(),
        observed_counts.len()
    );

    // d_sids where observed == 0: every trial trivially exceeds (0 >= 0)
    let zero_obs_sids: Vec<u32> = db_sources
        .keys()
        .filter(|&&d| observed_counts.get(&d).copied().unwrap_or(0) == 0)
        .copied()
        .collect();

    let mut rng: StdRng = match cli.seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_entropy(),
    };

    let mut exceed_counts: HashMap<u32, u64> = HashMap::new();
    let mut total_trials: usize = 0;
    let mut milestone_idx = 0;

    while total_trials < total_target {
        let next_milestone = milestones[milestone_idx];
        let to_run = (next_milestone - total_trials).min(batch_size);

        let tmp = TempDir::new().expect("Failed to create temp directory");

        // Write shifted BED files for non-empty trials
        let mut n_written = 0usize;
        for _ in 0..to_run {
            let shift: i64 = rng.gen_range(-max_shift..=max_shift);
            let path = tmp.path().join(format!("trial_{:08}.bed", n_written));
            let f = fs::File::create(&path).expect("Failed to create trial BED file");
            let mut writer = BufWriter::new(f);
            if write_shifted_bed(&query_bed, &chrom_sizes, whitelist.as_ref(), shift, &mut writer) {
                n_written += 1;
            } else {
                // All intervals trimmed off-chromosome: remove the empty file
                // so it doesn't appear as a query row (query_database skips it,
                // but an empty .bed file is still a .bed file that may be listed).
                drop(writer);
                let _ = fs::remove_file(path);
            }
        }

        // Zero-obs sources: every trial in this batch counts (0 >= 0)
        let n_empty = to_run - n_written;
        for &d_sid in &zero_obs_sids {
            *exceed_counts.entry(d_sid).or_insert(0) += to_run as u64;
        }

        // Query non-empty trials
        if n_written > 0 {
            let config = QueryConfig::new(cli.db.clone(), tmp.path().to_path_buf());
            let batch_result = query_database(&config).expect("Batch query failed");

            // Sanity check: result rows should match non-empty files.
            // (Empty files are skipped by parse_file_batch.)
            let result_rows = batch_result.counts.rows();
            // Empty-trial rows contribute 0 overlap → don't count for positive-obs sids.
            // We only iterate over actually-queried rows.
            let _ = n_empty; // used above for zero-obs sids

            for row in 0..result_rows {
                if let Some(row_view) = batch_result.counts.outer_view(row) {
                    for (col, &count) in row_view.iter() {
                        let d_sid = col as u32;
                        let obs = observed_counts.get(&d_sid).copied().unwrap_or(0);
                        if obs > 0 && count >= obs {
                            *exceed_counts.entry(d_sid).or_insert(0) += 1;
                        }
                    }
                }
            }
        }

        total_trials += to_run;
        // tmp is dropped here; temp files are cleaned up

        // Print p-values for any milestones we've now reached
        while milestone_idx < milestones.len() && total_trials >= milestones[milestone_idx] {
            let milestone = milestones[milestone_idx];
            eprintln!("Reached milestone {} after {} trials", milestone, total_trials);
            print_pvalues(
                milestone,
                total_trials,
                &exceed_counts,
                &observed_counts,
                &db_sources,
            );
            milestone_idx += 1;
        }

        if milestone_idx >= milestones.len() {
            break;
        }
    }
}
