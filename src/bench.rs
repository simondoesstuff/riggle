//! Benchmark utilities for data generation and analysis
//!
//! Used by both criterion benchmarks and stress test binaries.

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::process::Command;

use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

/// Standard human chromosomes for random selection
const CHROMOSOMES: &[&str] = &[
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
    "chr20", "chr21", "chr22", "chrX", "chrY",
];

/// Configuration for synthetic BED file generation
#[derive(Debug, Clone)]
pub struct BedGenConfig {
    pub num_intervals: usize,
    pub genome_size: u32,
    pub min_len: u32,
    pub max_len: u32,
    pub seed: u64,
    pub sort: bool,
}

impl Default for BedGenConfig {
    fn default() -> Self {
        Self {
            num_intervals: 10_000,
            genome_size: 250_000_000, // 100Mb
            min_len: 100,
            max_len: 10_000,
            seed: 42,
            sort: false,
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

    pub fn with_sort(mut self, sort: bool) -> Self {
        self.sort = sort;
        self
    }
}

/// Generate a synthetic BED file with random intervals
pub fn generate_bed_file(path: &Path, config: &BedGenConfig) {
    let mut rng = StdRng::seed_from_u64(config.seed);
    let file = File::create(path).unwrap();
    let mut writer = BufWriter::with_capacity(256 * 1024, file);

    if config.sort {
        // Collect intervals, sort by (chrom, start), then write
        let mut intervals: Vec<(&str, u32, u32)> = Vec::with_capacity(config.num_intervals);
        for _ in 0..config.num_intervals {
            let chrom = CHROMOSOMES.choose(&mut rng).unwrap();
            let len = rng.gen_range(config.min_len..=config.max_len);
            let start = rng.gen_range(0..config.genome_size.saturating_sub(len));
            let end = start + len;
            intervals.push((chrom, start, end));
        }
        intervals.sort_by(|a, b| a.0.cmp(b.0).then(a.1.cmp(&b.1)));
        for (chrom, start, end) in intervals {
            writeln!(writer, "{}\t{}\t{}", chrom, start, end).unwrap();
        }
    } else {
        for _ in 0..config.num_intervals {
            let chrom = CHROMOSOMES.choose(&mut rng).unwrap();
            let len = rng.gen_range(config.min_len..=config.max_len);
            let start = rng.gen_range(0..config.genome_size.saturating_sub(len));
            let end = start + len;
            writeln!(writer, "{}\t{}\t{}", chrom, start, end).unwrap();
        }
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
    bgzip: bool,
    sort: bool,
) {
    (0..num_files).into_par_iter().for_each(|i| {
        let path = dir.join(format!("source_{}.bed", i));
        let config = BedGenConfig::default()
            .with_intervals(intervals_per_file)
            .with_size_range(min_len, max_len)
            .with_seed(base_seed + i as u64)
            .with_sort(sort);
        generate_bed_file(&path, &config);

        if bgzip {
            let status = Command::new("bgzip").arg("-f").arg(&path).status().unwrap();
            if !status.success() {
                panic!("bgzip failed for file: {:?}", path);
            }
        }
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

/// Size information for a single layer across all shards
#[derive(Debug, Clone)]
pub struct LayerSizeInfo {
    pub layer_id: u8,
    pub size_bytes: u64,
    pub num_files: usize,
}

/// Analyze index size and structure
pub fn analyze_index_size(
    db_path: &Path,
    input_path: &Path,
    total_intervals: usize,
) -> IndexSizeReport {
    let input_size = dir_size(input_path);
    let index_size = dir_size(db_path);

    // Collect per-layer stats by scanning shard subdirectories for layer_*.bin files.
    let mut layer_map: std::collections::HashMap<u8, (u64, usize)> =
        std::collections::HashMap::new();

    if let Ok(shards) = fs::read_dir(db_path) {
        for entry in shards.filter_map(|e| e.ok()) {
            let shard_path = entry.path();
            if !shard_path.is_dir() {
                continue;
            }
            if let Ok(files) = fs::read_dir(&shard_path) {
                for fentry in files.filter_map(|e| e.ok()) {
                    let fname = fentry.file_name().to_string_lossy().to_string();
                    if let Some(rest) = fname.strip_prefix("layer_") {
                        if let Some(id_str) = rest.strip_suffix(".bin") {
                            if let Ok(layer_id) = id_str.parse::<u8>() {
                                let size = fentry.metadata().map(|m| m.len()).unwrap_or(0);
                                let entry = layer_map.entry(layer_id).or_insert((0, 0));
                                entry.0 += size;
                                entry.1 += 1;
                            }
                        }
                    }
                }
            }
        }
    }

    let mut layer_breakdown: Vec<LayerSizeInfo> = layer_map
        .into_iter()
        .map(|(layer_id, (size_bytes, num_files))| LayerSizeInfo {
            layer_id,
            size_bytes,
            num_files,
        })
        .collect();

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
                "  layer_{:<2} | {:>8.2} MB | {:>5} shards",
                layer.layer_id,
                layer.size_bytes as f64 / 1_000_000.0,
                layer.num_files,
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
            self.operation, self.total_elements, self.duration_secs, throughput, unit
        );
    }
}
