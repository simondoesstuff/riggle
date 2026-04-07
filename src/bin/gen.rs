use clap::Parser;
use std::fs;
use std::path::PathBuf;
use tqdm::{Style, tqdm};

use riggle::bench::{BedGenConfig, TimingResult, generate_bed_file, generate_bed_files_parallel};

/// CLI tool for generating synthetic BED files for benchmarks and stress tests
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Output path. If --num-files > 1, this must be a directory.
    #[arg(short, long)]
    pub output: PathBuf,

    /// Number of BED files to generate (generates in parallel if > 1)
    #[arg(short = 'f', long, default_value_t = 1)]
    pub num_files: usize,

    /// Number of intervals per file
    #[arg(short = 'n', long, default_value_t = 10_000)]
    pub num_intervals: usize,

    /// Genome size (Note: fixed to 250Mb in parallel generation due to function signature)
    #[arg(short = 'g', long, default_value_t = 250_000_000)]
    pub genome_size: u32,

    /// Minimum interval length
    #[arg(long, default_value_t = 100)]
    pub min_len: u32,

    /// Maximum interval length
    #[arg(long, default_value_t = 10_000)]
    pub max_len: u32,

    /// Base random seed
    #[arg(short, long, default_value_t = 42)]
    pub seed: u64,

    /// Compress using bgzip
    #[arg(short, long, default_value_t = false)]
    pub compress: bool,
}

fn main() {
    // Parse command line arguments
    let cli = Cli::parse();

    println!("Starting BED file generation...");

    if cli.num_files == 1 {
        // Single file generation
        if let Some(parent) = cli.output.parent() {
            if !parent.as_os_str().is_empty() {
                fs::create_dir_all(parent).expect("Failed to create parent directories");
            }
        }

        let config = BedGenConfig {
            num_intervals: cli.num_intervals,
            genome_size: cli.genome_size,
            min_len: cli.min_len,
            max_len: cli.max_len,
            seed: cli.seed,
        };

        let start = std::time::Instant::now();

        for _ in tqdm(0..1)
            .style(Style::Balloon)
            .desc(Some(format!("Generating BED file: {:?}", cli.output)))
        {
            generate_bed_file(&cli.output, &config);
        }

        let duration = start.elapsed().as_secs_f64();

        println!("Successfully generated 1 file at: {:?}", cli.output);
        let timing = TimingResult::new("generate_single_bed", cli.num_intervals, duration);
        timing.print();
    } else {
        // Parallel generation
        fs::create_dir_all(&cli.output).expect("Failed to create output directory");

        let start = std::time::Instant::now();

        for _ in tqdm(0..1)
            .style(Style::Balloon)
            .desc(Some(format!(
                "Generating {} files in parallel to {:?}",
                cli.num_files, cli.output
            )))
        {
            generate_bed_files_parallel(
                &cli.output,
                cli.num_files,
                cli.num_intervals,
                cli.min_len,
                cli.max_len,
                cli.seed,
                cli.compress,
            );
        }

        let duration = start.elapsed().as_secs_f64();

        let total_intervals = cli.num_files * cli.num_intervals;
        println!(
            "Successfully generated {} files in {:?}",
            cli.num_files, cli.output
        );

        let timing = TimingResult::new("generate_parallel_beds", total_intervals, duration);
        timing.print();
    }
}
