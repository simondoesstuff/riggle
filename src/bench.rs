//! Benchmark utilities for data generation and analysis
//!
//! Used by both criterion benchmarks and stress test binaries.

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

/// Configuration for synthetic BED file generation
#[derive(Debug, Clone)]
pub struct BedGenConfig {
    pub num_intervals: usize,
    pub genome_size: u32,
    pub min_len: u32,
    pub max_len: u32,
    pub seed: u64,
}

impl Default for BedGenConfig {
    fn default() -> Self {
        Self {
            num_intervals: 10_000,
            genome_size: 100_000_000, // 100Mb
            min_len: 100,
            max_len: 10_000,
            seed: 42,
        }
    }
}

impl BedGenConfig {
    pub fn with_intervals(mut self, n: usize) -> Self {
        self.num_intervals = n;
        self
    }

    pub fn with_size_range(mut self, min: u32, max: u32) -> Self {
        self.min_len = min;
        self.max_len = max;
        self
    }

    pub fn with_seed(mut self, seed: u64) -> Self {
        self.seed = seed;
        self
    }
}

/// Generate a synthetic BED file with random intervals
pub fn generate_bed_file(path: &Path, config: &BedGenConfig) {
    let mut rng = StdRng::seed_from_u64(config.seed);
    let file = File::create(path).unwrap();
    let mut writer = BufWriter::with_capacity(256 * 1024, file);

    for _ in 0..config.num_intervals {
        let len = rng.gen_range(config.min_len..=config.max_len);
        let start = rng.gen_range(0..config.genome_size.saturating_sub(len));
        let end = start + len;
        writeln!(writer, "chr1\t{}\t{}", start, end).unwrap();
    }
    writer.flush().unwrap();
}

/// Generate multiple BED files in parallel
pub fn generate_bed_files_parallel(
    dir: &Path,
    num_files: usize,
    intervals_per_file: usize,
    min_len: u32,
    max_len: u32,
    base_seed: u64,
) {
    (0..num_files).into_par_iter().for_each(|i| {
        let path = dir.join(format!("source_{}.bed", i));
        let config = BedGenConfig::default()
            .with_intervals(intervals_per_file)
            .with_size_range(min_len, max_len)
            .with_seed(base_seed + i as u64);
        generate_bed_file(&path, &config);
    });
}

/// Calculate total size of a directory recursively
pub fn dir_size(path: &Path) -> u64 {
    let mut size = 0;
    if path.is_dir() {
        for entry in fs::read_dir(path).unwrap() {
            let entry = entry.unwrap();
            let path = entry.path();
            if path.is_dir() {
                size += dir_size(&path);
            } else {
                size += entry.metadata().map(|m| m.len()).unwrap_or(0);
            }
        }
    }
    size
}

/// Index size analysis result
#[derive(Debug, Clone)]
pub struct IndexSizeReport {
    pub total_intervals: usize,
    pub input_size_bytes: u64,
    pub index_size_bytes: u64,
    pub bytes_per_interval: f64,
    pub expansion_ratio: f64,
    pub layer_breakdown: Vec<LayerSizeInfo>,
}

/// Size information for a single layer
#[derive(Debug, Clone)]
pub struct LayerSizeInfo {
    pub layer_id: u8,
    pub size_bytes: u64,
    pub num_chunks: usize,
    pub avg_chunk_bytes: f64,
}

/// Analyze index size and structure
pub fn analyze_index_size(db_path: &Path, input_path: &Path, total_intervals: usize) -> IndexSizeReport {
    let input_size = dir_size(input_path);
    let index_size = dir_size(db_path);

    let mut layer_breakdown = Vec::new();

    for entry in fs::read_dir(db_path).unwrap() {
        let entry = entry.unwrap();
        let path = entry.path();
        let name = entry.file_name().to_string_lossy().to_string();

        if path.is_dir() && name.starts_with("layer_") {
            let layer_id: u8 = name.strip_prefix("layer_").unwrap().parse().unwrap_or(0);
            let size = dir_size(&path);
            let chunk_count = fs::read_dir(&path).map(|d| d.count()).unwrap_or(0);

            layer_breakdown.push(LayerSizeInfo {
                layer_id,
                size_bytes: size,
                num_chunks: chunk_count,
                avg_chunk_bytes: if chunk_count > 0 { size as f64 / chunk_count as f64 } else { 0.0 },
            });
        }
    }

    layer_breakdown.sort_by_key(|l| l.layer_id);

    IndexSizeReport {
        total_intervals,
        input_size_bytes: input_size,
        index_size_bytes: index_size,
        bytes_per_interval: index_size as f64 / total_intervals as f64,
        expansion_ratio: index_size as f64 / input_size as f64,
        layer_breakdown,
    }
}

impl IndexSizeReport {
    /// Print a formatted report
    pub fn print(&self) {
        println!(
            "{:>12} intervals | index: {:>8.2} MB | {:.1} bytes/interval | {:.2}x input",
            self.total_intervals,
            self.index_size_bytes as f64 / 1_000_000.0,
            self.bytes_per_interval,
            self.expansion_ratio
        );
    }

    /// Print detailed layer breakdown
    pub fn print_breakdown(&self) {
        println!("\nLayer Breakdown:");
        for layer in &self.layer_breakdown {
            println!(
                "  layer_{:<2} | {:>8.2} MB | {:>5} chunks | {:>8.2} KB/chunk avg",
                layer.layer_id,
                layer.size_bytes as f64 / 1_000_000.0,
                layer.num_chunks,
                layer.avg_chunk_bytes / 1000.0
            );
        }
        println!(
            "  {:9} | {:>8.2} MB",
            "TOTAL",
            self.index_size_bytes as f64 / 1_000_000.0
        );
    }
}

/// Timing result for a benchmark run
#[derive(Debug, Clone)]
pub struct TimingResult {
    pub operation: String,
    pub total_elements: usize,
    pub duration_secs: f64,
    pub throughput_per_sec: f64,
}

impl TimingResult {
    pub fn new(operation: &str, total_elements: usize, duration_secs: f64) -> Self {
        Self {
            operation: operation.to_string(),
            total_elements,
            duration_secs,
            throughput_per_sec: total_elements as f64 / duration_secs,
        }
    }

    pub fn print(&self) {
        let (throughput, unit) = if self.throughput_per_sec >= 1_000_000.0 {
            (self.throughput_per_sec / 1_000_000.0, "M")
        } else if self.throughput_per_sec >= 1_000.0 {
            (self.throughput_per_sec / 1_000.0, "K")
        } else {
            (self.throughput_per_sec, "")
        };

        println!(
            "{:20} | {:>12} elements | {:>8.3}s | {:>8.2} {}elem/s",
            self.operation,
            self.total_elements,
            self.duration_secs,
            throughput,
            unit
        );
    }
}
