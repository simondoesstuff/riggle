//! Benchmarks for Riggle build and query pipelines
//!
//! Run with: cargo bench
//! Results saved to target/criterion/
//!
//! For large-scale stress tests, use the stress binary instead:
//!   cargo run --release --bin stress -- suite --scale 100

use std::path::Path;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main, Throughput};
use tempfile::TempDir;

use riggle::bench::{
    BedGenConfig, analyze_index_size, generate_bed_file, generate_bed_files_parallel,
};
use riggle::tasks::{BuildConfig, QueryConfig, build_database, query_database};

/// Benchmark database building at various scales
fn bench_build(c: &mut Criterion) {
    let mut group = c.benchmark_group("build");
    group.sample_size(10);

    // Test configurations: (num_files, intervals_per_file)
    let configs = [
        (5, 1_000),     // 5k intervals
        (10, 5_000),    // 50k intervals
        (20, 10_000),   // 200k intervals
        (50, 20_000),   // 1M intervals
    ];

    for (num_files, intervals_per_file) in configs {
        let total_intervals = num_files * intervals_per_file;
        group.throughput(Throughput::Elements(total_intervals as u64));

        group.bench_with_input(
            BenchmarkId::new("intervals", total_intervals),
            &(num_files, intervals_per_file),
            |b, &(nf, ipf)| {
                b.iter_with_setup(
                    || {
                        let input_dir = TempDir::new().unwrap();
                        let db_dir = TempDir::new().unwrap();
                        generate_bed_files_parallel(input_dir.path(), nf, ipf, 100, 10_000, 42);
                        (input_dir, db_dir)
                    },
                    |(input_dir, db_dir)| {
                        let config = BuildConfig::new(
                            input_dir.path().to_path_buf(),
                            db_dir.path().to_path_buf(),
                        );
                        build_database(&config).unwrap();
                    },
                );
            },
        );
    }

    group.finish();
}

/// Benchmark querying at various scales
fn bench_query(c: &mut Criterion) {
    let mut group = c.benchmark_group("query");

    // Build a fixed database for query benchmarks
    let input_dir = TempDir::new().unwrap();
    let db_dir = TempDir::new().unwrap();

    // 20 source files with 10k intervals each = 200k total database intervals
    generate_bed_files_parallel(input_dir.path(), 20, 10_000, 100, 10_000, 42);
    let build_config = BuildConfig::new(
        input_dir.path().to_path_buf(),
        db_dir.path().to_path_buf(),
    );
    build_database(&build_config).unwrap();

    // Test various query sizes
    let query_sizes = [100, 500, 1_000, 5_000, 10_000];

    for num_queries in query_sizes {
        group.throughput(Throughput::Elements(num_queries as u64));

        group.bench_with_input(
            BenchmarkId::new("queries", num_queries),
            &num_queries,
            |b, &nq| {
                b.iter_with_setup(
                    || {
                        let query_dir = TempDir::new().unwrap();
                        let query_path = query_dir.path().join("query.bed");
                        let config = BedGenConfig::default()
                            .with_intervals(nq)
                            .with_size_range(500, 50_000)
                            .with_seed(12345);
                        generate_bed_file(&query_path, &config);
                        (query_dir, query_path)
                    },
                    |(_query_dir, query_path)| {
                        let config = QueryConfig::new(
                            db_dir.path().to_path_buf(),
                            query_path,
                        );
                        query_database(&config).unwrap();
                    },
                );
            },
        );
    }

    group.finish();
}

/// Benchmark with varying interval size distributions
fn bench_interval_sizes(c: &mut Criterion) {
    let mut group = c.benchmark_group("interval_sizes");

    // Different size distributions to test layer partitioning
    let size_configs = [
        ("small", 50, 500),       // All in early layers
        ("medium", 500, 5_000),   // Middle layers
        ("large", 5_000, 50_000), // Later layers
        ("mixed", 50, 50_000),    // Spans all layers
    ];

    for (name, min_len, max_len) in size_configs {
        let input_dir = TempDir::new().unwrap();

        // Generate 10 files with 5k intervals each (parallel)
        generate_bed_files_parallel(input_dir.path(), 10, 5_000, min_len, max_len, 42);

        group.throughput(Throughput::Elements(50_000));

        group.bench_function(BenchmarkId::new("build", name), |b| {
            b.iter_with_setup(
                || TempDir::new().unwrap(),
                |out_dir| {
                    let config = BuildConfig::new(
                        input_dir.path().to_path_buf(),
                        out_dir.path().to_path_buf(),
                    );
                    build_database(&config).unwrap();
                },
            );
        });
    }

    group.finish();
}

/// Benchmark query scaling with database size
fn bench_query_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("query_scaling");
    group.sample_size(20);

    // Build databases of increasing size
    let db_sizes = [
        (10, 5_000),   // 50k intervals
        (20, 10_000),  // 200k intervals
        (50, 10_000),  // 500k intervals
    ];

    for (num_files, intervals_per_file) in db_sizes {
        let total = num_files * intervals_per_file;
        let input_dir = TempDir::new().unwrap();
        let db_dir = TempDir::new().unwrap();

        generate_bed_files_parallel(input_dir.path(), num_files, intervals_per_file, 100, 10_000, 42);
        let build_config = BuildConfig::new(
            input_dir.path().to_path_buf(),
            db_dir.path().to_path_buf(),
        );
        build_database(&build_config).unwrap();

        // Fixed query size
        let query_dir = TempDir::new().unwrap();
        let query_path = query_dir.path().join("query.bed");
        let config = BedGenConfig::default()
            .with_intervals(1_000)
            .with_size_range(500, 50_000)
            .with_seed(12345);
        generate_bed_file(&query_path, &config);

        group.throughput(Throughput::Elements(total as u64));

        group.bench_function(BenchmarkId::new("db_intervals", total), |b| {
            b.iter(|| {
                let config = QueryConfig::new(
                    db_dir.path().to_path_buf(),
                    query_path.clone(),
                );
                query_database(&config).unwrap();
            });
        });
    }

    group.finish();
}

/// Analyze index size (prints summary, tracked by criterion)
fn bench_index_size(c: &mut Criterion) {
    let mut group = c.benchmark_group("index_size");
    group.sample_size(10);

    println!("\n=== Index Size Analysis ===\n");

    // Test configurations: (num_files, intervals_per_file, min_len, max_len, description)
    let configs: [(usize, usize, u32, u32, &str); 5] = [
        (10, 5_000, 100, 10_000, "50k mixed"),
        (20, 10_000, 100, 10_000, "200k mixed"),
        (50, 20_000, 100, 10_000, "1M mixed"),
        (10, 5_000, 50, 500, "50k small"),
        (10, 5_000, 5_000, 50_000, "50k large"),
    ];

    for (num_files, intervals_per_file, min_len, max_len, desc) in configs {
        let total_intervals = num_files * intervals_per_file;
        let input_dir = TempDir::new().unwrap();
        let db_dir = TempDir::new().unwrap();

        generate_bed_files_parallel(input_dir.path(), num_files, intervals_per_file, min_len, max_len, 42);

        let config = BuildConfig::new(
            input_dir.path().to_path_buf(),
            db_dir.path().to_path_buf(),
        );
        build_database(&config).unwrap();

        let report = analyze_index_size(db_dir.path(), input_dir.path(), total_intervals);

        println!(
            "{:15} | {:>10} intervals | index: {:>8.2} MB | {:.1} bytes/interval | {:.2}x input",
            desc,
            total_intervals,
            report.index_size_bytes as f64 / 1_000_000.0,
            report.bytes_per_interval,
            report.expansion_ratio
        );

        group.throughput(Throughput::Elements(total_intervals as u64));
        group.bench_function(BenchmarkId::new("bytes", desc), |b| {
            b.iter(|| report.index_size_bytes)
        });
    }

    println!();
    group.finish();
}

criterion_group!(
    benches,
    bench_build,
    bench_query,
    bench_interval_sizes,
    bench_query_scaling,
    bench_index_size,
);
criterion_main!(benches);
