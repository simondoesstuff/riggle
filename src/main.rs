use std::collections::HashMap;
use std::path::PathBuf;

use clap::{Parser, Subcommand};

use riggle::io::Meta;
use riggle::stats::{StatResult, StatsOutput};
use riggle::tasks::{AddConfig, QueryConfig, add_to_database, query_database};

#[derive(Parser)]
#[command(name = "riggle")]
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

        /// Compute FFT-based p-values for all overlapping pairs.
        /// Requires that the database was built with this version of riggle
        /// (Fourier spectra must be cached under {db}/fourier/).
        #[arg(long)]
        stats: bool,

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
        Commands::Add { input, db, batch_size } => run_add(input, db, batch_size),
        Commands::Query { db, query, output, stats, batch_size } => {
            run_query(db, query, output, stats, batch_size)
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
) -> Result<(), Box<dyn std::error::Error>> {
    let mut config = AddConfig::new(input, db);
    config.batch_size = batch_size;
    add_to_database(&config)?;
    Ok(())
}

fn run_query(
    db: PathBuf,
    query: PathBuf,
    output: PathBuf,
    stats: bool,
    batch_size: Option<usize>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut config = QueryConfig::new(db, query);
    config.batch_size = batch_size;
    config.stats = stats;
    let result = query_database(&config)?;

    // Build a fast lookup from (q_sid, d_sid) → (observed_bins, p_value).
    let pvalue_map: HashMap<(usize, u32), (f64, f64)> = result
        .pvalues
        .iter()
        .map(|pv| ((pv.query_id, pv.db_sid), (pv.observed_bins, pv.p_value)))
        .collect();

    // Flatten the sparse overlap matrix into StatResult records.
    let db_sources = &result.db_sources;
    let query_names = &result.query_names;
    let mut stat_results: Vec<StatResult> = result
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
                    let (observed_bins, p_value) = pvalue_map
                        .get(&(q_sid, d_sid as u32))
                        .map(|&(o, p)| (Some(o), Some(p)))
                        .unwrap_or((None, None));
                    StatResult {
                        query_name: q_name.clone(),
                        db_name,
                        overlap_count,
                        observed_bins,
                        p_value,
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect();

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
