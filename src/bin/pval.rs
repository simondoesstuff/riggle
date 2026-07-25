//! pval — analytic NB p-value for genomic interval overlap.
//!
//! Usage: pval <a.bed> <b.bed>

use std::env;
use std::path::Path;

use chuckle::fourier::{
    DEFAULT_MOMENTS_EPS, DepthMap, build_depth_moments, build_query_chrom_data,
    compute_analytic_stats, coverage_dot_product, parse_bed_as_map,
};

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: pval <a.bed> <b.bed>");
        std::process::exit(1);
    }

    let a_path = Path::new(&args[1]);
    let b_path = Path::new(&args[2]);

    let a_bed =
        parse_bed_as_map(a_path).unwrap_or_else(|e| panic!("Cannot parse {}: {e}", args[1]));
    let b_bed =
        parse_bed_as_map(b_path).unwrap_or_else(|e| panic!("Cannot parse {}: {e}", args[2]));

    let n_a: usize = a_bed.values().map(|v| v.len()).sum();
    let n_b: usize = b_bed.values().map(|v| v.len()).sum();
    eprintln!("A: {} intervals across {} chroms", n_a, a_bed.len());
    eprintln!("B: {} intervals across {} chroms", n_b, b_bed.len());

    let a_dm = DepthMap::build(&a_bed);
    let b_dm = DepthMap::build(&b_bed);

    let a_data = build_query_chrom_data(&a_bed, None);
    let b_moments = build_depth_moments(&b_dm, DEFAULT_MOMENTS_EPS);

    let observed = coverage_dot_product(&a_dm, &b_dm);

    match compute_analytic_stats(&a_data, &b_moments, observed) {
        None => eprintln!("No shared chromosomes."),
        Some((p_value, llr)) => {
            println!("---");
            println!(
                "observed overlap : {:.0} bins  ({:.0} bp)",
                observed,
                observed * 100.0
            );
            println!("p-value (right)  : {:.4e}", p_value);
            if let Some(llr) = llr {
                println!("LLR              : {:.4}", llr);
            } else {
                println!("LLR              : N/A (NB fit not feasible)");
            }
        }
    }
}
