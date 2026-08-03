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
/// `resolve_chrom(chrom)` is called once per query chromosome and returns an
/// opaque handle `H` (e.g., a pre-resolved chromosome index).  `lookup_l(h,
/// l_int)` is then called once per unique integer block size on that
/// chromosome and returns `(mean, var)`.  Separating the two calls lets
/// callers pay the O(N_chroms) chrom scan once per chrom rather than once
/// per unique block size.
///
/// Returns `(p_value, llr)`, or `None` when no (query, DB) chromosome pair shares
/// any moments.
pub fn compute_analytic_stats<H: Copy>(
    q_data: &[QueryChromData],
    resolve_chrom: impl Fn(&str) -> Option<H>,
    lookup_l: impl Fn(H, usize) -> Option<(f64, f64)>,
    observed: f64,
) -> Option<(f64, Option<f64>)> {
    let mut mu_total = 0.0f64;
    let mut var_total = 0.0f64;
    // μ_raw: Σ mean(l_int) / l_int ≈ Σ coverage_fraction_i (scale-free).
    let mut mu_raw = 0.0f64;
    let mut any_shared = false;

    for q_cd in q_data {
        let chrom_handle = match resolve_chrom(&q_cd.chrom) {
            Some(h) => h,
            None => continue,
        };
        for g in &q_cd.grouped {
            if let Some((mean_l, var_l)) = lookup_l(chrom_handle, g.l_int) {
                mu_total += mean_l * g.sum_scale;
                var_total += var_l * g.sum_scale2;
                mu_raw += mean_l / g.l_int as f64 * g.count as f64;
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

    // LLR is always ≥ 0 (KL divergence); direction comes from observed vs μ.
    // Enrichment (observed ≥ μ): p = erfc(√LLR)/2 ∈ (0, 0.5].
    // Depletion  (observed < μ): p = 1 − erfc(√LLR)/2 ∈ [0.5, 1).
    let p_value = match llr {
        Some(l) if l > 0.0 => {
            let tail = erfc(l.sqrt()) * 0.5;
            if observed >= mu_total { tail } else { 1.0 - tail }
        }
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

    fn make_resolve(moments: &[ChromMoments]) -> impl Fn(&str) -> Option<usize> + '_ {
        |chrom| moments.iter().position(|m| m.chrom == chrom)
    }

    fn make_lookup_l(moments: &[ChromMoments]) -> impl Fn(usize, usize) -> Option<(f64, f64)> + '_ {
        |idx, l_int| moments.get(idx)?.lookup(l_int as f64)
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

        assert!(compute_analytic_stats(&q_data, make_resolve(&db_moments), make_lookup_l(&db_moments), 0.0).is_none());
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
            compute_analytic_stats(&q_data, make_resolve(&db_moments), make_lookup_l(&db_moments), 0.0).unwrap();
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

        let (p_value, llr) =
            compute_analytic_stats(&q_data, make_resolve(&db_moments), make_lookup_l(&db_moments), 25_000.0).unwrap();
        assert!(llr.map_or(true, |l| l >= 0.0), "llr={llr:?}");
        assert!(p_value < 0.5, "enriched pair should have p < 0.5, got {p_value}");
    }

    #[test]
    fn test_compute_analytic_stats_depleted_gives_high_p() {
        // Use mock closures to control μ/σ² exactly: mean=20, var=100 (overdispersed).
        // Query: one interval of 1000bp = 10 bins on chr22 → scale=1, mu_total=20, var_total=100.
        let mut q_bed = BedMap::new();
        q_bed.insert("chr22".to_string(), vec![(0, 1_000)]);
        let q_data = build_query_chrom_data(&q_bed);

        let resolve_all = |_: &str| -> Option<()> { Some(()) };
        let fixed_moments = |_: (), _: usize| -> Option<(f64, f64)> { Some((20.0, 100.0)) };

        // observed = 0: NB p_o = 1 fails bounds → LLR=None → p=1.0.
        let (p_zero, _) = compute_analytic_stats(&q_data, resolve_all, fixed_moments, 0.0).unwrap();
        assert_eq!(p_zero, 1.0, "observed=0 → p=1.0, got {p_zero}");

        // observed = 10 < mu=20 (moderate depletion): p ∈ (0.5, 1.0).
        let (p_depleted, llr) = compute_analytic_stats(&q_data, resolve_all, fixed_moments, 10.0).unwrap();
        assert!(p_depleted > 0.5, "depletion (obs=10<mu=20) must have p>0.5, got {p_depleted} (llr={llr:?})");
        assert!(p_depleted < 1.0, "depletion must have p<1.0 (not saturated), got {p_depleted}");
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
