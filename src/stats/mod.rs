use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::matrix::SparseMatrix;

/// Result of statistical analysis for a single query-source pair
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StatResult {
    /// Query identifier
    pub query_id: usize,
    /// Database source ID
    pub db_sid: u32,
    /// Number of overlapping intervals
    pub overlap_count: u32,
    /// P-value from Fisher's exact test
    pub p_value: f64,
    /// Odds ratio
    pub odds_ratio: f64,
}

/// Statistical output for all queries
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StatsOutput {
    pub results: Vec<StatResult>,
}

impl StatsOutput {
    /// Convert to JSON string
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }
}

/// Compute statistics for intersection results
///
/// # Arguments
/// * `counts` - Sparse matrix of intersection counts (rows = queries, cols = db sources)
/// * `query_sizes` - Size of each query in base pairs
/// * `db_sizes` - Size of each database source in base pairs
/// * `genome_size` - Total genome size for background calculation
pub fn compute_statistics(
    counts: &SparseMatrix,
    query_sizes: &[u64],
    db_sizes: &HashMap<u32, u64>,
    genome_size: u64,
) -> StatsOutput {
    let mut results = Vec::new();

    // Iterate over non-zero entries in the sparse matrix
    for (row_idx, row_vec) in counts.outer_iterator().enumerate() {
        for (col_idx, &count) in row_vec.iter() {
            if count == 0 {
                continue;
            }

            let db_sid = col_idx as u32;
            let query_size = query_sizes.get(row_idx).copied().unwrap_or(0);
            let db_size = db_sizes.get(&db_sid).copied().unwrap_or(0);

            // Compute Fisher's exact test
            let (p_value, odds_ratio) = fishers_exact_test(
                count as u64,
                query_size,
                db_size,
                genome_size,
            );

            results.push(StatResult {
                query_id: row_idx,
                db_sid,
                overlap_count: count,
                p_value,
                odds_ratio,
            });
        }
    }

    // Sort by p-value
    results.sort_by(|a, b| a.p_value.partial_cmp(&b.p_value).unwrap_or(std::cmp::Ordering::Equal));

    StatsOutput { results }
}

/// Compute Fisher's exact test for 2x2 contingency table
///
/// The contingency table is:
/// ```text
///                  | Overlaps DB | Doesn't overlap DB |
/// Query region     |     a       |         b          |
/// Non-query region |     c       |         d          |
/// ```
///
/// # Arguments
/// * `overlap_count` - Number of overlapping intervals (a)
/// * `query_size` - Total query size in base pairs
/// * `db_size` - Total database source size in base pairs
/// * `genome_size` - Total genome size
///
/// # Returns
/// (p_value, odds_ratio)
pub fn fishers_exact_test(
    overlap_count: u64,
    query_size: u64,
    db_size: u64,
    genome_size: u64,
) -> (f64, f64) {
    if genome_size == 0 || query_size == 0 || db_size == 0 {
        return (1.0, 1.0);
    }

    // Build contingency table
    // a = overlap_count (observed overlaps)
    // b = query_size - overlap_count (query that doesn't overlap)
    // c = db_size - overlap_count (db that doesn't overlap with query)
    // d = genome_size - query_size - db_size + overlap_count

    let a = overlap_count as f64;
    let b = (query_size.saturating_sub(overlap_count)) as f64;
    let c = (db_size.saturating_sub(overlap_count)) as f64;
    let d = (genome_size
        .saturating_sub(query_size)
        .saturating_sub(db_size)
        .saturating_add(overlap_count)) as f64;

    // Odds ratio = (a*d) / (b*c)
    let odds_ratio = if b * c > 0.0 {
        (a * d) / (b * c)
    } else if a > 0.0 {
        f64::INFINITY
    } else {
        1.0
    };

    // Compute p-value using log-space hypergeometric
    // P-value is the probability of observing this many or more overlaps by chance
    let p_value = hypergeometric_pvalue(a as u64, (a + b) as u64, (a + c) as u64, (a + b + c + d) as u64);

    (p_value, odds_ratio)
}

/// Compute hypergeometric p-value (one-tailed, right)
///
/// P(X >= k) where X ~ Hypergeometric(N, K, n)
/// - N = population size
/// - K = success states in population
/// - n = number of draws
/// - k = observed successes
fn hypergeometric_pvalue(k: u64, n: u64, big_k: u64, big_n: u64) -> f64 {
    if big_n == 0 || n == 0 || big_k == 0 {
        return 1.0;
    }

    // Compute cumulative probability P(X >= k)
    // Using log-space to avoid overflow
    let mut p_value = 0.0;
    let max_possible = n.min(big_k);

    for i in k..=max_possible {
        let log_prob = log_hypergeometric_pmf(i, n, big_k, big_n);
        p_value += log_prob.exp();

        // Early termination if probability becomes negligible
        if log_prob < -50.0 && i > k {
            break;
        }
    }

    p_value.min(1.0)
}

/// Compute log of hypergeometric PMF
/// P(X = k) = C(K, k) * C(N-K, n-k) / C(N, n)
fn log_hypergeometric_pmf(k: u64, n: u64, big_k: u64, big_n: u64) -> f64 {
    log_binomial(big_k, k) + log_binomial(big_n - big_k, n - k) - log_binomial(big_n, n)
}

/// Compute log of binomial coefficient C(n, k) = n! / (k! * (n-k)!)
fn log_binomial(n: u64, k: u64) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    if k == 0 || k == n {
        return 0.0;
    }

    // Use Stirling's approximation for large values
    // log(C(n,k)) = log(n!) - log(k!) - log((n-k)!)
    log_factorial(n) - log_factorial(k) - log_factorial(n - k)
}

/// Compute log factorial using Stirling's approximation for large values
fn log_factorial(n: u64) -> f64 {
    if n <= 1 {
        return 0.0;
    }

    // For small values, compute directly
    if n <= 20 {
        let mut result = 0.0;
        for i in 2..=n {
            result += (i as f64).ln();
        }
        return result;
    }

    // Stirling's approximation: ln(n!) ≈ n*ln(n) - n + 0.5*ln(2πn)
    let n_f = n as f64;
    n_f * n_f.ln() - n_f + 0.5 * (2.0 * std::f64::consts::PI * n_f).ln()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_log_factorial() {
        // Small values should be exact
        assert!((log_factorial(0) - 0.0).abs() < 1e-10);
        assert!((log_factorial(1) - 0.0).abs() < 1e-10);
        assert!((log_factorial(5) - (120.0_f64).ln()).abs() < 1e-10);

        // Larger values use approximation
        let log_20_factorial = log_factorial(20);
        let expected = (1..=20).map(|i| (i as f64).ln()).sum::<f64>();
        assert!((log_20_factorial - expected).abs() < 1e-10);
    }

    #[test]
    fn test_log_binomial() {
        // C(5, 2) = 10
        let result = log_binomial(5, 2).exp();
        assert!((result - 10.0).abs() < 1e-10);

        // C(10, 5) = 252
        let result = log_binomial(10, 5).exp();
        assert!((result - 252.0).abs() < 1e-8);

        // Edge cases
        assert_eq!(log_binomial(5, 0).exp(), 1.0);
        assert_eq!(log_binomial(5, 5).exp(), 1.0);
        assert!(log_binomial(5, 6).is_infinite() && log_binomial(5, 6).is_sign_negative());
    }

    #[test]
    fn test_fishers_exact_basic() {
        // Test with known values
        // If everything overlaps perfectly, p-value should be very low
        let (p_value, odds_ratio) = fishers_exact_test(100, 100, 100, 1000);
        assert!(p_value < 0.05);
        assert!(odds_ratio > 1.0);
    }

    #[test]
    fn test_fishers_exact_no_overlap() {
        // No overlap at all
        let (p_value, _odds_ratio) = fishers_exact_test(0, 100, 100, 1000);
        // P-value should be 1.0 or close to it (not enriched)
        assert!(p_value >= 0.5);
    }

    #[test]
    fn test_fishers_exact_edge_cases() {
        // Zero genome size
        let (p_value, odds_ratio) = fishers_exact_test(10, 100, 100, 0);
        assert_eq!(p_value, 1.0);
        assert_eq!(odds_ratio, 1.0);

        // Zero query size
        let (p_value, odds_ratio) = fishers_exact_test(0, 0, 100, 1000);
        assert_eq!(p_value, 1.0);
        assert_eq!(odds_ratio, 1.0);
    }

    #[test]
    fn test_compute_statistics() {
        use crate::matrix::{condense_to_sparse, BitwiseMask, DenseMatrix};

        let mut dense = DenseMatrix::new(2, 2);
        let mut mask = BitwiseMask::new(2, 2);

        dense.set(0, 0, 10);
        dense.set(1, 1, 5);
        mask.flag(0, 0);
        mask.flag(1, 1);

        let counts = condense_to_sparse(&dense, &mask);

        let query_sizes = vec![100, 200];
        let db_sizes: HashMap<u32, u64> = [(0, 150), (1, 250)].into_iter().collect();
        let genome_size = 10000;

        let stats = compute_statistics(&counts, &query_sizes, &db_sizes, genome_size);

        assert_eq!(stats.results.len(), 2);
        assert!(stats.results.iter().all(|r| r.p_value >= 0.0 && r.p_value <= 1.0));
    }

    #[test]
    fn test_stats_output_json() {
        let output = StatsOutput {
            results: vec![StatResult {
                query_id: 0,
                db_sid: 1,
                overlap_count: 10,
                p_value: 0.001,
                odds_ratio: 5.5,
            }],
        };

        let json = output.to_json().unwrap();
        assert!(json.contains("query_id"));
        assert!(json.contains("p_value"));
    }
}
