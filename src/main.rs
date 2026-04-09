use std::path::PathBuf;

use clap::{Parser, Subcommand};

use riggle::io::Meta;
use riggle::stats::compute_statistics;
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

        /// Genome size for statistical calculations (default: 3 billion for human)
        #[arg(long, default_value = "3000000000")]
        genome_size: u64,
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
        Commands::Add { input, db } => run_add(input, db),
        Commands::Query {
            db,
            query,
            output,
            genome_size,
        } => run_query(db, query, output, genome_size),
        Commands::Info { db } => run_info(db),
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

fn run_add(input: PathBuf, db: PathBuf) -> Result<(), Box<dyn std::error::Error>> {
    let config = AddConfig::new(input, db);
    add_to_database(&config)?;
    Ok(())
}

fn run_query(
    db: PathBuf,
    query: PathBuf,
    output: PathBuf,
    genome_size: u64,
) -> Result<(), Box<dyn std::error::Error>> {
    let config = QueryConfig::new(db, query);
    let result = query_database(&config)?;

    // Compute sizes for Fisher's exact test
    let query_sizes: Vec<u64> = result
        .query_sources
        .iter()
        .map(|s| s.count as u64)
        .collect();

    let db_sizes: std::collections::HashMap<u32, u64> = result
        .db_sources
        .iter()
        .map(|(k, _)| {
            // Count intervals per D_SID from the sparse matrix column sums
            let col_sum: u32 = result
                .counts
                .outer_iterator()
                .map(|row| row.get(*k as usize).copied().unwrap_or(0))
                .sum();
            (*k, col_sum as u64)
        })
        .collect();

    let stats = compute_statistics(&result.counts, &query_sizes, &db_sizes, genome_size);
    let json = serde_json::to_string_pretty(&stats)?;
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
