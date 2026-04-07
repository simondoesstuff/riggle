//! Stress test binary for large-scale benchmarking
//!
//! Usage:
//!   cargo run --release --bin stress -- build --intervals 100000000
//!   cargo run --release --bin stress -- query --db /path/to/db --queries 100000
//!   cargo run --release --bin stress -- size --intervals 10000000

use std::path::PathBuf;
use std::time::Instant;

use clap::{Parser, Subcommand};
use tempfile::TempDir;

use riggle::bench::{
    BedGenConfig, TimingResult, analyze_index_size, generate_bed_file, generate_bed_files_parallel,
};
use riggle::tasks::{BuildConfig, QueryConfig, build_database, query_database};

#[derive(Parser)]
#[command(name = "stress")]
#[command(about = "Stress test Riggle at scale")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Benchmark database building at scale
    Build {
        /// Total number of intervals to index
        #[arg(short, long, default_value = "10_000_000")]
        intervals: usize,

        /// Number of source files to generate
        #[arg(short, long, default_value = "100")]
        files: usize,

        /// Minimum interval length
        #[arg(long, default_value = "100")]
        min_len: u32,

        /// Maximum interval length
        #[arg(long, default_value = "10000")]
        max_len: u32,

        /// Output database path (uses temp dir if not specified)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Keep the database after benchmark (requires --output)
        #[arg(short, long)]
        keep: bool,
    },

    /// Benchmark query performance at scale
    Query {
        /// Path to existing database
        #[arg(short, long)]
        db: PathBuf,

        /// Number of query intervals
        #[arg(short, long, default_value = "10000")]
        queries: usize,

        /// Minimum query length
        #[arg(long, default_value = "500")]
        min_len: u32,

        /// Maximum query length
        #[arg(long, default_value = "50000")]
        max_len: u32,

        /// Number of iterations
        #[arg(short, long, default_value = "3")]
        iterations: usize,
    },

    /// Analyze index size at various scales
    Size {
        /// Total number of intervals
        #[arg(short, long, default_value = "10000000")]
        intervals: usize,

        /// Number of source files
        #[arg(short, long, default_value = "100")]
        files: usize,

        /// Minimum interval length
        #[arg(long, default_value = "100")]
        min_len: u32,

        /// Maximum interval length
        #[arg(long, default_value = "10000")]
        max_len: u32,

        /// Show layer breakdown
        #[arg(short, long)]
        breakdown: bool,
    },

    /// Run full benchmark suite at specified scale
    Suite {
        /// Scale factor (intervals = scale * 1_000_000)
        #[arg(short, long, default_value = "10")]
        scale: usize,
    },
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Build {
            intervals,
            files,
            min_len,
            max_len,
            output,
            keep,
        } => {
            run_build_benchmark(intervals, files, min_len, max_len, output, keep);
        }
        Commands::Query {
            db,
            queries,
            min_len,
            max_len,
            iterations,
        } => {
            run_query_benchmark(&db, queries, min_len, max_len, iterations);
        }
        Commands::Size {
            intervals,
            files,
            min_len,
            max_len,
            breakdown,
        } => {
            run_size_analysis(intervals, files, min_len, max_len, breakdown);
        }
        Commands::Suite { scale } => {
            run_full_suite(scale);
        }
    }
}

fn run_build_benchmark(
    total_intervals: usize,
    num_files: usize,
    min_len: u32,
    max_len: u32,
    output: Option<PathBuf>,
    keep: bool,
) {
    let intervals_per_file = total_intervals / num_files;

    println!("=== Build Benchmark ===");
    println!("Total intervals: {}", total_intervals);
    println!("Files: {} x {} intervals", num_files, intervals_per_file);
    println!("Interval size: {}-{} bp", min_len, max_len);
    println!();

    // Generate input data
    let input_dir = TempDir::new().unwrap();
    println!("Generating {} BED files...", num_files);
    let gen_start = Instant::now();
    generate_bed_files_parallel(
        input_dir.path(),
        num_files,
        intervals_per_file,
        min_len,
        max_len,
        42,
        false,
        false,
    );
    let gen_time = gen_start.elapsed();
    TimingResult::new("Data generation", total_intervals, gen_time.as_secs_f64()).print();

    // Build database
    let db_dir = if let Some(ref path) = output {
        std::fs::create_dir_all(path).unwrap();
        None
    } else {
        Some(TempDir::new().unwrap())
    };
    let db_path = output
        .clone()
        .unwrap_or_else(|| db_dir.as_ref().unwrap().path().to_path_buf());

    println!("\nBuilding database...");
    let build_start = Instant::now();
    let config = BuildConfig::new(input_dir.path().to_path_buf(), db_path.clone());
    build_database(&config).unwrap();
    let build_time = build_start.elapsed();

    let result = TimingResult::new("Database build", total_intervals, build_time.as_secs_f64());
    result.print();

    // Report index size
    let report = analyze_index_size(&db_path, input_dir.path(), total_intervals);
    println!(
        "\nIndex size: {:.2} MB ({:.1} bytes/interval, {:.2}x input)",
        report.index_size_bytes as f64 / 1_000_000.0,
        report.bytes_per_interval,
        report.expansion_ratio
    );

    if keep && output.is_some() {
        println!("\nDatabase saved to: {}", db_path.display());
    }
}

fn run_query_benchmark(
    db_path: &PathBuf,
    num_queries: usize,
    min_len: u32,
    max_len: u32,
    iterations: usize,
) {
    println!("=== Query Benchmark ===");
    println!("Database: {}", db_path.display());
    println!("Queries: {}", num_queries);
    println!("Query size: {}-{} bp", min_len, max_len);
    println!("Iterations: {}", iterations);
    println!();

    // Generate query file
    let query_dir = TempDir::new().unwrap();
    let query_path = query_dir.path().join("query.bed");
    let config = BedGenConfig::default()
        .with_intervals(num_queries)
        .with_size_range(min_len, max_len)
        .with_seed(12345);
    generate_bed_file(&query_path, &config);

    // Run queries
    let mut times = Vec::with_capacity(iterations);
    for i in 0..iterations {
        let start = Instant::now();
        let query_config = QueryConfig::new(db_path.clone(), query_path.clone());
        let result = query_database(&query_config).unwrap();
        let elapsed = start.elapsed();
        times.push(elapsed.as_secs_f64());

        println!(
            "  Run {}: {:.3}s ({} queries, {} sources, {} non-zero)",
            i + 1,
            elapsed.as_secs_f64(),
            result.counts.rows(),
            result.db_sources.len(),
            result.counts.nnz()
        );
    }

    let avg_time = times.iter().sum::<f64>() / times.len() as f64;
    let min_time = times.iter().cloned().fold(f64::INFINITY, f64::min);

    println!();
    TimingResult::new("Query (avg)", num_queries, avg_time).print();
    TimingResult::new("Query (best)", num_queries, min_time).print();
}

fn run_size_analysis(
    total_intervals: usize,
    num_files: usize,
    min_len: u32,
    max_len: u32,
    show_breakdown: bool,
) {
    let intervals_per_file = total_intervals / num_files;

    println!("=== Index Size Analysis ===");
    println!("Total intervals: {}", total_intervals);
    println!("Interval size: {}-{} bp", min_len, max_len);
    println!();

    // Generate and build
    let input_dir = TempDir::new().unwrap();
    let db_dir = TempDir::new().unwrap();

    print!("Generating data... ");
    std::io::Write::flush(&mut std::io::stdout()).unwrap();
    generate_bed_files_parallel(
        input_dir.path(),
        num_files,
        intervals_per_file,
        min_len,
        max_len,
        42,
        false,
        false,
    );
    println!("done");

    print!("Building index... ");
    std::io::Write::flush(&mut std::io::stdout()).unwrap();
    let config = BuildConfig::new(input_dir.path().to_path_buf(), db_dir.path().to_path_buf());
    build_database(&config).unwrap();
    println!("done");
    println!();

    let report = analyze_index_size(db_dir.path(), input_dir.path(), total_intervals);
    report.print();

    if show_breakdown {
        report.print_breakdown();
    }
}

fn run_full_suite(scale: usize) {
    let base_intervals = scale * 1_000_000;

    println!("╔══════════════════════════════════════════════════════════════╗");
    println!(
        "║           RIGGLE STRESS TEST SUITE - {}M SCALE             ║",
        scale
    );
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    // Test different interval size distributions
    let size_configs = [
        ("small (50-500bp)", 50, 500),
        ("medium (500-5K)", 500, 5_000),
        ("large (5K-50K)", 5_000, 50_000),
        ("mixed (100-10K)", 100, 10_000),
    ];

    println!("=== Build Performance by Interval Size ===\n");

    for (name, min_len, max_len) in size_configs {
        let input_dir = TempDir::new().unwrap();
        let db_dir = TempDir::new().unwrap();

        let num_files = 100;
        let intervals_per_file = base_intervals / num_files;

        // Generate
        generate_bed_files_parallel(
            input_dir.path(),
            num_files,
            intervals_per_file,
            min_len,
            max_len,
            42,
            false,
            false,
        );

        // Build
        let start = Instant::now();
        let config = BuildConfig::new(input_dir.path().to_path_buf(), db_dir.path().to_path_buf());
        build_database(&config).unwrap();
        let elapsed = start.elapsed();

        let report = analyze_index_size(db_dir.path(), input_dir.path(), base_intervals);

        let throughput = base_intervals as f64 / elapsed.as_secs_f64() / 1_000_000.0;
        println!(
            "{:20} | {:>6.2}s | {:>5.2}M/s | {:>6.1} MB | {:>5.1} B/iv | {:.2}x",
            name,
            elapsed.as_secs_f64(),
            throughput,
            report.index_size_bytes as f64 / 1_000_000.0,
            report.bytes_per_interval,
            report.expansion_ratio
        );
    }

    // Query scaling test
    println!("\n=== Query Performance Scaling ===\n");

    // Build a reference database
    let input_dir = TempDir::new().unwrap();
    let db_dir = TempDir::new().unwrap();
    generate_bed_files_parallel(
        input_dir.path(),
        100,
        base_intervals / 100,
        100,
        10_000,
        42,
        false,
        false,
    );
    let config = BuildConfig::new(input_dir.path().to_path_buf(), db_dir.path().to_path_buf());
    build_database(&config).unwrap();

    let query_sizes = [1_000, 10_000, 100_000];

    for num_queries in query_sizes {
        let query_dir = TempDir::new().unwrap();
        let query_path = query_dir.path().join("query.bed");
        let qconfig = BedGenConfig::default()
            .with_intervals(num_queries)
            .with_size_range(500, 50_000)
            .with_seed(12345);
        generate_bed_file(&query_path, &qconfig);

        let start = Instant::now();
        let query_config = QueryConfig::new(db_dir.path().to_path_buf(), query_path);
        let result = query_database(&query_config).unwrap();
        let elapsed = start.elapsed();

        let throughput = num_queries as f64 / elapsed.as_secs_f64() / 1_000.0;
        println!(
            "{:>8} queries | {:>6.3}s | {:>6.1}K/s | {:>8} nnz",
            num_queries,
            elapsed.as_secs_f64(),
            throughput,
            result.counts.nnz()
        );
    }

    println!("\n=== Complete ===");
}
