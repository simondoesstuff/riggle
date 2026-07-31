use serde::{Deserialize, Serialize};

/// One overlapping (query interval, DB interval) pair from `--intervals` mode.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IntervalHit {
    pub query_name: String,
    pub db_name: String,
    pub chrom: String,
    pub query_start: u32,
    pub query_end: u32,
    pub db_start: u32,
    pub db_end: u32,
    /// Exact intersection length in base pairs.
    pub intersection_bp: u32,
}

/// Full output for an `--intervals` query run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IntervalsOutput {
    pub hits: Vec<IntervalHit>,
}

impl IntervalsOutput {
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }
}

/// Statistical result for one (query, database source) pair.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StatResult {
    pub query_name: String,
    pub db_name: String,
    /// Number of overlapping interval pairs from the sweep phase.
    pub overlap_count: u32,
    /// Base-pair overlap in 100 bp bins (dot product of coverage arrays).
    /// Present only when `--stats` is enabled.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub observed_bins: Option<f64>,
    /// Right-tailed p-value under the rigid-body shift null model.
    /// Present only when `--stats` is enabled.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub p_value: Option<f64>,
    /// Saddlepoint log-likelihood ratio: NB null tilted to the observation.
    /// Present only when `--stats` is enabled and the NB fit succeeds.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub llr: Option<f64>,
}

/// Full statistical output for a query run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StatsOutput {
    pub results: Vec<StatResult>,
}

impl StatsOutput {
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stats_output_json_no_pvalue() {
        let out = StatsOutput {
            results: vec![StatResult {
                query_name: "q.bed".into(),
                db_name: "db.bed".into(),
                overlap_count: 42,
                observed_bins: None,
                p_value: None,
                llr: None,
            }],
        };
        let json = out.to_json().unwrap();
        assert!(json.contains("overlap_count"));
        assert!(!json.contains("p_value"));
        assert!(!json.contains("llr"));
    }

    #[test]
    fn test_stats_output_json_with_pvalue() {
        let out = StatsOutput {
            results: vec![StatResult {
                query_name: "q.bed".into(),
                db_name: "db.bed".into(),
                overlap_count: 42,
                observed_bins: Some(420.0),
                p_value: Some(1e-5),
                llr: Some(12.3),
            }],
        };
        let json = out.to_json().unwrap();
        assert!(json.contains("p_value"));
        assert!(json.contains("observed_bins"));
        assert!(json.contains("llr"));
    }
}
