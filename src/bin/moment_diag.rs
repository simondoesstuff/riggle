/// Diagnostic: examine var/mean ratio across the moment store and for a GWAS query.
///
/// Usage:
///   cargo run --bin moment_diag --release -- <momentmap.rkyv> [gwas.bed ...]
///
/// Output:
///   1. For L = 1, 2, 5, 10, 50, 100, 500: fraction of (sid, chrom) pairs with var > mean,
///      and median var/mean ratio.
///   2. For each supplied GWAS BED file: aggregate mu_total vs var_total over all SIDs
///      that appear in the store, so we can see whether the Poisson fallback is always needed.

use std::collections::HashMap;
use std::env;
use std::path::Path;

use chuckle::fourier::{BedMap, QueryChromData, build_query_chrom_data, parse_bed_as_map};
use chuckle::io::MappedMomentStore;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: moment_diag <momentmap.rkyv> [gwas.bed ...]");
        std::process::exit(1);
    }

    let store_path = Path::new(&args[1]);
    let store = MappedMomentStore::open(store_path).expect("open momentmap");

    let probe_ls: &[usize] = &[1, 2, 5, 10, 25, 50, 100, 200, 500, 1000];

    // ── Part 1: sweep over SIDs and chromosomes ───────────────────────────────

    // Collect (mean, var) at each probe_l across all (sid, chrom) pairs.
    let mut by_l: HashMap<usize, Vec<(f64, f64)>> = HashMap::new();

    let max_sid = 2000u32; // sample up to this many SIDs
    let mut n_sids = 0usize;

    for sid in 0..max_sid {
        let Some(sid_m) = store.get_sid(sid) else { continue };
        n_sids += 1;
        let chroms = ["chr1","chr2","chr3","chr6","chr11","chr17","chr22"];
        for &chrom in &chroms {
            for &l in probe_ls {
                if let Some((mean, var)) = sid_m.lookup(chrom, l as f64) {
                    by_l.entry(l).or_default().push((mean, var));
                }
            }
        }
    }

    println!("=== Moment store overdispersion (over {} SIDs) ===", n_sids);
    println!("{:>6}  {:>8}  {:>12}  {:>12}  {:>10}",
             "L", "n_obs", "median_mean", "median_var/mean", "frac_overdispersed");

    for &l in probe_ls {
        let Some(pairs) = by_l.get(&l) else { continue };
        if pairs.is_empty() { continue }

        let mut ratios: Vec<f64> = pairs.iter()
            .filter(|&&(m, _)| m > 1e-10)
            .map(|&(m, v)| v / m)
            .collect();
        ratios.sort_by(f64::total_cmp);

        let mut means: Vec<f64> = pairs.iter().map(|&(m, _)| m).collect();
        means.sort_by(f64::total_cmp);

        let frac_over = ratios.iter().filter(|&&r| r > 1.0).count() as f64 / ratios.len() as f64;
        let med_ratio = if ratios.is_empty() { f64::NAN } else { ratios[ratios.len() / 2] };
        let med_mean = if means.is_empty() { f64::NAN } else { means[means.len() / 2] };

        println!("{:>6}  {:>8}  {:>12.4}  {:>12.4}  {:>10.1}%",
                 l, pairs.len(), med_mean, med_ratio, frac_over * 100.0);
    }

    // ── Part 2: aggregate GWAS query moments ─────────────────────────────────

    for gwas_path in &args[2..] {
        let path = Path::new(gwas_path);
        let bed: BedMap = match parse_bed_as_map(path) {
            Ok(b) => b,
            Err(e) => { eprintln!("skip {gwas_path}: {e}"); continue }
        };
        let q_data: Vec<QueryChromData> = build_query_chrom_data(&bed);
        let n_intervals: usize = q_data.iter().map(|c| c.interval_lengths.len()).sum();

        println!("\n=== GWAS query: {} ({} intervals) ===", gwas_path, n_intervals);
        println!("{:>8}  {:>12}  {:>12}  {:>10}  {:>10}",
                 "sid", "mu_total", "var_total", "var/mu", "overdispersed");

        let mut n_over = 0usize;
        let mut n_under = 0usize;
        let mut n_nb_fail = 0usize;
        let mut sample_rows: Vec<(u32, f64, f64)> = Vec::new();

        for sid in 0..max_sid {
            let Some(sid_m) = store.get_sid(sid) else { continue };
            let mut mu_total = 0.0f64;
            let mut var_total = 0.0f64;
            let mut any = false;

            for q_cd in &q_data {
                for &l in &q_cd.interval_lengths {
                    let l_int = (l.ceil() as usize).max(1);
                    let scale = l / l_int as f64;
                    if let Some((mean_l, var_l)) = sid_m.lookup(&q_cd.chrom, l_int as f64) {
                        mu_total += mean_l * scale;
                        var_total += var_l * scale * scale;
                        any = true;
                    }
                }
            }

            if !any { continue }

            if var_total > mu_total { n_over += 1; } else { n_under += 1; }
            if mu_total > 0.0 && var_total <= mu_total { n_nb_fail += 1; }

            if sample_rows.len() < 5 {
                sample_rows.push((sid, mu_total, var_total));
            }
        }

        for (sid, mu, var) in &sample_rows {
            let ratio = if *mu > 1e-12 { var / mu } else { f64::NAN };
            let od = if var > mu { "yes" } else { "no" };
            println!("{:>8}  {:>12.5e}  {:>12.5e}  {:>10.4}  {:>10}", sid, mu, var, ratio, od);
        }

        let total = n_over + n_under;
        println!("  overdispersed: {}/{} ({:.0}%)",
                 n_over, total, n_over as f64 / total as f64 * 100.0);
        println!("  NB fit fails (var<=mu): {}/{} ({:.0}%)",
                 n_nb_fail, total, n_nb_fail as f64 / total as f64 * 100.0);
    }
}
