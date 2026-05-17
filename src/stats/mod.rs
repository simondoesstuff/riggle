use serde::{Deserialize, Serialize};

/// Statistical result for one (query, database source) pair.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StatResult {
    pub query_name: String,
    pub db_name: String,
    /// Number of overlapping interval pairs from the sweep phase.
    pub overlap_count: u32,
    /// Base-pair overlap at shift 0 (in 100 bp bins).  Present only when
    /// `--stats` is enabled and the Fourier cache exists for the DB file.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub observed_bins: Option<f64>,
    /// Right-tailed p-value under the rigid-body shift null model.
    /// Present only when `--stats` is enabled.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub p_value: Option<f64>,
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
            }],
        };
        let json = out.to_json().unwrap();
        assert!(json.contains("overlap_count"));
        assert!(!json.contains("p_value"));
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
            }],
        };
        let json = out.to_json().unwrap();
        assert!(json.contains("p_value"));
        assert!(json.contains("observed_bins"));
    }
}
