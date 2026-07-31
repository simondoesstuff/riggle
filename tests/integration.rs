/// Integration tests at real scale using seeded synthetic BED data.
///
/// Ground-truth totals were computed by running each test once and printing
/// the observed values, then hardcoded here.  The seed-deterministic generator
/// guarantees reproducibility across platforms.

use chuckle::bench::{BedGenConfig, generate_bed_file};
use chuckle::tasks::build::AddConfig;
use chuckle::tasks::{QueryConfig, QueryMode, add_to_database, query_database};
use tempfile::TempDir;

/// Build a seeded scenario: `n_db` source files + `n_query` query files, each
/// with the given interval counts, on a reduced genome size for test speed.
///
/// DB seeds: 100, 101, …  Query seeds: 200, 201, …
fn setup(
    n_db: usize,
    db_intervals: usize,
    n_query: usize,
    query_intervals: usize,
    genome_size: u32,
    stats: bool,
) -> (TempDir, TempDir, TempDir) {
    let db_src = TempDir::new().unwrap();
    for i in 0..n_db {
        generate_bed_file(
            &db_src.path().join(format!("source_{i}.bed")),
            &BedGenConfig::default()
                .with_intervals(db_intervals)
                .with_genome_size(genome_size)
                .with_seed(100 + i as u64),
        );
    }

    let db = TempDir::new().unwrap();
    let mut add_cfg = AddConfig::new(db_src.path().to_path_buf(), db.path().to_path_buf());
    add_cfg.compute_stats = stats;
    add_to_database(&add_cfg).unwrap();

    let query_dir = TempDir::new().unwrap();
    for i in 0..n_query {
        generate_bed_file(
            &query_dir.path().join(format!("query_{i}.bed")),
            &BedGenConfig::default()
                .with_intervals(query_intervals)
                .with_genome_size(genome_size)
                .with_seed(200 + i as u64),
        );
    }

    (db_src, db, query_dir)
}

// ── Intervals mode ────────────────────────────────────────────────────────────

#[test]
fn test_intervals_total_hits_and_bp() {
    // 10 DB files × 2000 intervals; 1 query × 1000 intervals; 50 Mbp genome.
    let (_src, db, query_dir) = setup(10, 2000, 1, 1000, 50_000_000, false);

    let mut cfg = QueryConfig::new(db.path().to_path_buf(), query_dir.path().to_path_buf());
    cfg.mode = QueryMode::Intervals;
    let result = query_database(&cfg).unwrap();
    let hits = &result.interval_hits;

    assert!(
        hits.iter().all(|h| h.intersection_bp > 0),
        "every hit must have positive intersection"
    );

    let total_bp: u64 = hits.iter().map(|h| h.intersection_bp as u64).sum();

    assert_eq!(hits.len(), 158);
    assert_eq!(total_bp, 413_652);
}

#[test]
fn test_intervals_all_hits_valid() {
    // Smoke-test the structural invariants across a small scenario.
    let (_src, db, query_dir) = setup(3, 500, 2, 300, 20_000_000, false);

    let mut cfg = QueryConfig::new(db.path().to_path_buf(), query_dir.path().to_path_buf());
    cfg.mode = QueryMode::Intervals;
    let result = query_database(&cfg).unwrap();

    for h in &result.interval_hits {
        assert!(!h.chrom.is_empty());
        assert!(h.query_start < h.query_end, "query interval malformed");
        assert!(h.db_start < h.db_end, "db interval malformed");
        assert!(h.intersection_bp > 0);
        // Intersection must be contained within both intervals.
        let lo = h.query_start.max(h.db_start);
        let hi = h.query_end.min(h.db_end);
        assert_eq!(hi - lo, h.intersection_bp, "intersection_bp mismatch");
        assert!(!h.query_name.is_empty());
        assert!(!h.db_name.is_empty());
    }
}

// ── Counts mode ───────────────────────────────────────────────────────────────

#[test]
fn test_counts_nonzero_pairs() {
    // 5 DB files × 1000 intervals; 2 query files × 500 intervals.
    let (_src, db, query_dir) = setup(5, 1000, 2, 500, 30_000_000, false);

    let cfg = QueryConfig::new(db.path().to_path_buf(), query_dir.path().to_path_buf());
    let result = query_database(&cfg).unwrap();

    let nnz = result.counts.nnz();
    let total_count: u32 = result.counts.data().iter().copied().sum();

    assert_eq!(nnz, 10);          // all 2×5 pairs have at least one overlap
    assert_eq!(total_count, 58);
}

// ── Stats mode ────────────────────────────────────────────────────────────────

#[test]
fn test_stats_self_query_is_most_enriched() {
    // Build with one source file, query with that same file: should be the
    // most significant (lowest p-value) among all results.
    let db_src = TempDir::new().unwrap();
    let self_bed = db_src.path().join("self.bed");
    generate_bed_file(
        &self_bed,
        &BedGenConfig::default()
            .with_intervals(2000)
            .with_genome_size(30_000_000)
            .with_seed(77),
    );
    // Add a second, unrelated source file so there's something to compare against.
    generate_bed_file(
        &db_src.path().join("other.bed"),
        &BedGenConfig::default()
            .with_intervals(2000)
            .with_genome_size(30_000_000)
            .with_seed(88),
    );

    let db = TempDir::new().unwrap();
    let mut add_cfg = AddConfig::new(db_src.path().to_path_buf(), db.path().to_path_buf());
    add_cfg.compute_stats = true;
    add_to_database(&add_cfg).unwrap();

    let mut cfg = QueryConfig::new(db.path().to_path_buf(), self_bed.clone());
    cfg.mode = QueryMode::Stats;
    let result = query_database(&cfg).unwrap();

    // The self-pair should have the lowest p-value of the two results.
    let self_pv = result
        .pvalues
        .iter()
        .filter(|pv| {
            result
                .db_sources
                .get(&pv.db_sid)
                .map_or(false, |n| n == "self.bed")
        })
        .map(|pv| pv.p_value)
        .next()
        .expect("self pair missing");

    let other_pv = result
        .pvalues
        .iter()
        .filter(|pv| {
            result
                .db_sources
                .get(&pv.db_sid)
                .map_or(false, |n| n == "other.bed")
        })
        .map(|pv| pv.p_value)
        .next()
        .expect("other pair missing");

    assert!(
        self_pv < other_pv,
        "self p={self_pv:.2e} should be lower than other p={other_pv:.2e}"
    );
    assert!(self_pv < 1e-10, "self-query p-value should be highly significant, got {self_pv:.2e}");
}
