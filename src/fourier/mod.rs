mod depth;
mod genome;
mod moments;
mod stats;

pub use depth::{
    ChromDepthMap, DepthMap, bed_map_v, coverage_dot_product, intervals_to_bed_map,
    parse_bed_as_map,
};
pub use genome::{BedMap, hg38_chrom_sizes};
pub use moments::{
    ChromMoments, EPS, QueryChromData, T_THRESHOLD, build_chrom_moments, build_depth_moments,
    build_query_chrom_data, compact_index, mean_interval_bins,
};
pub use stats::compute_analytic_stats;
