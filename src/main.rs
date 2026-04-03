use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

use clap::{Parser, Subcommand};
use indicatif::{ProgressBar, ProgressStyle};

use riggle::io::MasterHeader;
use riggle::stats::compute_statistics;
use riggle::tasks::{build_database, query_database, BuildConfig, QueryConfig};

#[derive(Parser)]
#[command(name = "riggle")]
#[command(about = "Statistical interval intersection engine", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Build a database from BED files
    Build {
        /// Input directory containing BED files
        #[arg(short, long)]
        input: PathBuf,

        /// Output directory for the database
        #[arg(short, long)]
        output: PathBuf,
    },

    /// Query a database with a BED file
    Query {
        /// Path to the database directory
        #[arg(short, long)]
        db: PathBuf,

        /// Query BED file
        #[arg(short, long)]
        query: PathBuf,

        /// Output JSON file for results
        #[arg(short, long)]
        output: PathBuf,

        /// Genome size for statistical calculations (default: 3 billion for human)
        #[arg(long, default_value = "3000000000")]
        genome_size: u64,
    },

    /// Show database information
    Info {
        /// Path to the database directory
        #[arg(short, long)]
        db: PathBuf,
    },
}

fn main() {
    let cli = Cli::parse();

    let result = match cli.command {
        Commands::Build { input, output } => run_build(input, output),
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

fn run_build(input: PathBuf, output: PathBuf) -> Result<(), Box<dyn std::error::Error>> {
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} {msg}")
            .unwrap(),
    );
    pb.set_message("Building database...");
    pb.enable_steady_tick(std::time::Duration::from_millis(100));

    let config = BuildConfig::new(input, output.clone());
    build_database(&config)?;

    pb.finish_with_message(format!("Database built at {}", output.display()));
    Ok(())
}

fn run_query(
    db: PathBuf,
    query: PathBuf,
    output: PathBuf,
    genome_size: u64,
) -> Result<(), Box<dyn std::error::Error>> {
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} {msg}")
            .unwrap(),
    );
    pb.set_message("Querying database...");
    pb.enable_steady_tick(std::time::Duration::from_millis(100));

    let config = QueryConfig::new(db.clone(), query);
    let result = query_database(&config)?;

    pb.set_message("Computing statistics...");

    // Load master header for database sizes
    let header_path = db.join("header.json");
    let header_content = fs::read_to_string(&header_path)?;
    let master_header: MasterHeader = serde_json::from_str(&header_content)?;

    // Compute query sizes (assuming uniform for now - could be improved)
    let query_sizes: Vec<u64> = vec![1000; result.counts.rows()]; // Placeholder

    // Get database sizes
    let db_sizes: HashMap<u32, u64> = master_header
        .sid_map
        .iter()
        .map(|(k, v)| (*k, v.total_bases))
        .collect();

    let stats = compute_statistics(&result.counts, &query_sizes, &db_sizes, genome_size);

    // Write results
    let json = stats.to_json()?;
    fs::write(&output, json)?;

    pb.finish_with_message(format!(
        "Results written to {} ({} significant pairs)",
        output.display(),
        stats.results.len()
    ));

    Ok(())
}

fn run_info(db: PathBuf) -> Result<(), Box<dyn std::error::Error>> {
    let header_path = db.join("header.json");
    let header_content = fs::read_to_string(&header_path)?;
    let header: MasterHeader = serde_json::from_str(&header_content)?;

    println!("Database: {}", db.display());
    println!("Sources: {}", header.num_sources());
    println!("Max coordinate: {}", header.max_coord);
    println!("Layers: {}", header.layer_configs.len());
    println!();
    println!("Source files:");

    let mut sources: Vec<_> = header.sid_map.iter().collect();
    sources.sort_by_key(|(k, _)| *k);

    for (sid, meta) in sources {
        println!(
            "  [{}] {} - {} intervals, {} bp",
            sid, meta.name, meta.interval_count, meta.total_bases
        );
    }

    Ok(())
}
