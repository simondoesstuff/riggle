//! Analytic NB/Poisson statistics for genomic interval overlap.
//!
//! # Null model
//! Under the free-relative-gap null, each query interval slides independently
//! over the DB with a uniform random position.  For query interval of length L
//! bins on a chromosome where the DB has pre-computed moments M[L]:
//!
//!   μ_i    = M[L].mean × scale     (expected overlap in bins; scale = l/l_int)
//!   σ²_i   = M[L].var  × scale²
//!   μ_raw_i = M[L].mean / l_int    (≈ coverage fraction; scale-free)
//!
//! Genome-wide null: μ = Σ μ_i, σ² = Σ σ²_i, μ_raw = Σ μ_raw_i.
//!
//! # Dispatch
//! When σ² > μ (overdispersed): NB null on bin-overlap statistic.
//! When σ² ≤ μ (underdispersed, e.g. 1 bp SNPs or sparse single-bin marks):
//!   Poisson null on count units.  The effective scale factor
//!   s = μ / μ_raw ≈ mean interval length rescales observed_bins → observed_count
//!   so the LLR has the correct magnitude (not deflated by 0.01 for SNPs).

use super::moments::QueryChromData;

/// Compute analytic statistics for a (query, DB) pair.
///
/// `lookup(chrom, l_bins)` returns `(mean, var)` for a sliding window of
/// `l_bins` bins on `chrom`.  Returns `(p_value, llr)`.
pub fn compute_analytic_stats(
    q_data: &[QueryChromData],
    lookup: impl Fn(&str, f64) -> Option<(f64, f64)>,
    observed: f64,
) -> Option<(f64, Option<f64>)> {
    let mut mu_total = 0.0f64;
    let mut var_total = 0.0f64;
    // μ_raw: Σ mean(l_int) / l_int ≈ Σ coverage_fraction_i (scale-free).
    let mut mu_raw = 0.0f64;
    let mut any_shared = false;

    for q_cd in q_data {
        for &l in &q_cd.interval_lengths {
            let l_int = (l.ceil() as usize).max(1);
            let scale = l / l_int as f64;
            if let Some((mean_l, var_l)) = lookup(&q_cd.chrom, l_int as f64) {
                mu_total += mean_l * scale;
                var_total += var_l * scale * scale;
                mu_raw += mean_l / l_int as f64;
                any_shared = true;
            }
        }
    }

    if !any_shared {
        return None;
    }

    let llr = if var_total > mu_total {
        nb_llr(mu_total, var_total, observed)
    } else {
        // Underdispersed: Poisson on count units.
        // effective_scale ≈ mean interval length in bins; rescaling removes the
        // factor that would otherwise shrink LLR to scale × LLR_count.
        let effective_scale = if mu_raw > 0.0 { mu_total / mu_raw } else { 1.0 };
        let observed_count = observed / effective_scale;
        poisson_llr(mu_raw, observed_count)
    };

    let p_value = match llr {
        Some(l) if l > 0.0 => erfc(l.sqrt()) * 0.5,
        _ => 1.0,
    };

    Some((p_value, llr))
}

/// Saddlepoint LLR for the Poisson null; fallback when NB requires underdispersion.
///
/// LLR = O·ln(O/μ) − (O − μ); positive for enrichment (O > μ).
fn poisson_llr(mu: f64, observed: f64) -> Option<f64> {
    if mu <= 0.0 || observed <= 0.0 {
        return None;
    }
    Some(observed * (observed / mu).ln() - (observed - mu))
}

/// Saddlepoint LLR for the NB null tilted to observation `o`.
fn nb_llr(mu: f64, var: f64, observed: f64) -> Option<f64> {
    if mu <= 0.0 || var <= mu {
        return None;
    }
    let r = mu * mu / (var - mu);
    let p = mu / var;
    let p_o = r / (r + observed);
    if !(0.0 < p && p < 1.0 && 0.0 < p_o && p_o < 1.0) {
        return None;
    }
    let theta = ((1.0 - p_o) / (1.0 - p)).ln();
    Some(observed * theta + r * (p_o / p).ln())
}

/// Complementary error function (max error < 1.5e-7 for x ≥ 0).
fn erfc(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc(-x);
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let p = t
        * (0.254_829_592
            + t * (-0.284_496_736
                + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
    p * (-x * x).exp()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fourier::{
        BedMap, ChromMoments, DepthMap, build_depth_moments, build_query_chrom_data,
    };

    fn make_lookup(moments: &[ChromMoments]) -> impl Fn(&str, f64) -> Option<(f64, f64)> + '_ {
        |chrom, l| moments.iter().find(|m| m.chrom == chrom)?.lookup(l)
    }

    // ── compute_analytic_stats ────────────────────────────────────────────────

    #[test]
    fn test_compute_analytic_stats_no_shared_chrom() {
        let mut q_bed = BedMap::new();
        q_bed.insert("chr21".to_string(), vec![(10_000_000, 10_001_000)]);
        let q_data = build_query_chrom_data(&q_bed);

        let mut db_bed = BedMap::new();
        db_bed.insert("chr22".to_string(), vec![(10_000_000, 10_001_000)]);
        let db_dm = DepthMap::build(&db_bed);
        let db_moments = build_depth_moments(&db_dm);

        assert!(compute_analytic_stats(&q_data, make_lookup(&db_moments), 0.0).is_none());
    }

    #[test]
    fn test_compute_analytic_stats_zero_observed() {
        let mut q_bed = BedMap::new();
        q_bed.insert("chr22".to_string(), vec![(0, 1_000_000)]);
        let q_data = build_query_chrom_data(&q_bed);

        let mut db_bed = BedMap::new();
        db_bed.insert("chr22".to_string(), vec![(49_000_000, 50_000_000)]);
        let db_dm = DepthMap::build(&db_bed);
        let db_moments = build_depth_moments(&db_dm);

        let (p_value, _) =
            compute_analytic_stats(&q_data, make_lookup(&db_moments), 0.0).unwrap();
        assert_eq!(p_value, 1.0);
    }

    #[test]
    fn test_compute_analytic_stats_enriched() {
        let ivs: Vec<(u32, u32)> =
            (0..500u32).map(|i| (i * 10_000, i * 10_000 + 5_000)).collect();

        let mut db_bed = BedMap::new();
        db_bed.insert("chr22".to_string(), ivs.clone());
        let db_dm = DepthMap::build(&db_bed);
        let db_moments = build_depth_moments(&db_dm);

        let mut q_bed = BedMap::new();
        q_bed.insert("chr22".to_string(), ivs);
        let q_data = build_query_chrom_data(&q_bed);

        let (_p_value, llr) =
            compute_analytic_stats(&q_data, make_lookup(&db_moments), 25_000.0).unwrap();
        assert!(llr.map_or(true, |l| l >= 0.0), "llr={llr:?}");
    }

    // ── nb_llr ────────────────────────────────────────────────────────────────

    #[test]
    fn test_nb_llr_zero_at_null_mean() {
        let llr = nb_llr(4000.0, 8000.0, 4000.0).unwrap();
        assert!(llr.abs() < 1e-6, "llr={llr}");
    }

    #[test]
    fn test_nb_llr_positive_for_enrichment() {
        let llr = nb_llr(4000.0, 8000.0, 20_000.0).unwrap();
        assert!(llr > 0.0, "llr={llr}");
    }

    #[test]
    fn test_nb_llr_none_underdispersed() {
        assert!(nb_llr(2000.0, 1000.0, 1.0).is_none());
    }

    // ── erfc ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_erfc_known_values() {
        assert!((erfc(0.0) - 1.0).abs() < 1e-6);
        assert!((erfc(1.0) - 0.157_299_2).abs() < 1e-5, "erfc(1)={}", erfc(1.0));
        assert!(erfc(5.0) < 1e-10);
    }
}
