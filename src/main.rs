use std::path::PathBuf;

use clap::{Parser, Subcommand};

use chuckle::io::Meta;
use chuckle::stats::{IntervalsOutput, StatResult, StatsOutput};
use chuckle::tasks::{AddConfig, QueryConfig, QueryMode, add_to_database, query_database};

#[derive(Parser)]
#[command(name = "chuckle")]
#[command(about = "Statistical interval intersection engine", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Add BED files to the database (creates the database if it does not exist)
    Add {
        /// Input directory containing BED files to add
        #[arg(short, long)]
        input: PathBuf,

        /// Path to database directory
        #[arg(short, long)]
        db: PathBuf,

        /// Maximum number of BED files to hold in memory at once (default: all)
        #[arg(long)]
        batch_size: Option<usize>,

        /// Precompute FFT moment tables for analytic p-values at query time.
        /// Adds O(N log N) per chromosome at index time; required for --stats at query time.
        #[arg(long)]
        stats: bool,

        /// Only compute FFT moments for these chromosomes (default: all).
        /// Example: --chroms chr1 chr2 chr3
        #[arg(long, num_args(1..))]
        chroms: Vec<String>,
    },

    /// Query a database with BED file(s)
    Query {
        /// Path to the database directory
        #[arg(short, long)]
        db: PathBuf,

        /// Query BED file or directory containing BED files
        #[arg(short, long)]
        query: PathBuf,

        /// Output JSON file for results
        #[arg(short, long)]
        output: PathBuf,

        /// Compute analytic NB p-values and LLR for all overlapping pairs.
        /// Requires the database to have been built with --stats.
        #[arg(long, conflicts_with = "intervals")]
        stats: bool,

        /// Return the exact set of overlapping interval pairs with intersection lengths.
        #[arg(long, conflicts_with = "stats")]
        intervals: bool,

        /// Maximum number of query files to hold in memory at once (default: all)
        #[arg(long)]
        batch_size: Option<usize>,
    },

    /// Print database summary (shards, sources, layers)
    Info {
        /// Path to the database directory
        #[arg(short, long)]
        db: PathBuf,
    },
}

fn main() {
    let cli = Cli::parse();

    let result = match cli.command {
        Commands::Add { input, db, batch_size, stats, chroms } => run_add(input, db, batch_size, stats, chroms),
        Commands::Query { db, query, output, stats, intervals, batch_size } => {
            run_query(db, query, output, stats, intervals, batch_size)
        }
        Commands::Info { db } => run_info(db),
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

fn run_add(
    input: PathBuf,
    db: PathBuf,
    batch_size: Option<usize>,
    stats: bool,
    chroms: Vec<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut config = AddConfig::new(input, db);
    config.batch_size = batch_size;
    config.compute_stats = stats;
    if !chroms.is_empty() {
        config.chrom_whitelist = Some(chroms.into_iter().collect());
    }
    add_to_database(&config)?;
    Ok(())
}

fn run_query(
    db: PathBuf,
    query: PathBuf,
    output: PathBuf,
    stats: bool,
    intervals: bool,
    batch_size: Option<usize>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut config = QueryConfig::new(db, query);
    config.batch_size = batch_size;
    config.mode = if stats {
        QueryMode::Stats
    } else if intervals {
        QueryMode::Intervals
    } else {
        QueryMode::Counts
    };
    let result = query_database(&config)?;

    if intervals {
        let n_hits = result.interval_hits.len();
        let json = IntervalsOutput { hits: result.interval_hits }.to_json()?;
        std::fs::write(&output, json)?;
        println!("Query complete: {} interval hits → {}", n_hits, output.display());
        return Ok(());
    }

    let db_sources = &result.db_sources;
    let query_names = &result.query_names;

    // When stats are enabled, pvalues covers every (q_sid, d_sid) pair.
    // Without stats, fall back to the sparse matrix which only records non-zero overlaps.
    let mut stat_results: Vec<StatResult> = if stats {
        result
            .pvalues
            .iter()
            .map(|pv| {
                let q_name = query_names.get(pv.query_id).cloned().unwrap_or_default();
                let db_name = db_sources.get(&pv.db_sid).cloned().unwrap_or_default();
                let overlap_count = result
                    .counts
                    .get(pv.query_id, pv.db_sid as usize)
                    .copied()
                    .unwrap_or(0);
                StatResult {
                    query_name: q_name,
                    db_name,
                    overlap_count,
                    observed_bins: Some(pv.observed_bins),
                    p_value: Some(pv.p_value),
                    llr: pv.llr,
                }
            })
            .collect()
    } else {
        result
            .counts
            .outer_iterator()
            .enumerate()
            .flat_map(|(q_sid, row)| {
                let q_name = query_names.get(q_sid).cloned().unwrap_or_default();
                row.iter()
                    .map(|(d_sid, &overlap_count)| {
                        let db_name = db_sources
                            .get(&(d_sid as u32))
                            .cloned()
                            .unwrap_or_default();
                        StatResult {
                            query_name: q_name.clone(),
                            db_name,
                            overlap_count,
                            observed_bins: None,
                            p_value: None,
                            llr: None,
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect()
    };

    // Sort by p-value when available, otherwise by descending overlap count.
    stat_results.sort_by(|a, b| match (a.p_value, b.p_value) {
        (Some(pa), Some(pb)) => pa.partial_cmp(&pb).unwrap_or(std::cmp::Ordering::Equal),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => b.overlap_count.cmp(&a.overlap_count),
    });

    let json = StatsOutput { results: stat_results }.to_json()?;
    std::fs::write(&output, json)?;

    println!(
        "Query complete: {} query files × {} database sources → {}",
        result.query_sources.len(),
        result.db_sources.len(),
        output.display()
    );

    Ok(())
}

fn run_info(db: PathBuf) -> Result<(), Box<dyn std::error::Error>> {
    let meta = Meta::load(&db)?;

    println!("Database: {}", db.display());
    println!(
        "Layer config: min_size={}, growth_factor={}",
        meta.layer_config.min_size, meta.layer_config.growth_factor
    );
    println!("Layers used: {}", meta.num_layers);
    println!("Shards ({}): {}", meta.shards.len(), meta.shards.join(", "));
    println!("\nSources ({}):", meta.sid_map.len());
    let mut sids: Vec<u32> = meta.sid_map.keys().copied().collect();
    sids.sort();
    for sid in sids {
        println!("  SID {:4}: {}", sid, meta.sid_map[&sid].name);
    }

    Ok(())
}
