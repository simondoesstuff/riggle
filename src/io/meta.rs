use std::collections::HashMap;
use std::fs;
use std::path::Path;

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Errors from meta.json operations
#[derive(Debug, Error)]
pub enum MetaError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("JSON error: {0}")]
    Json(#[from] serde_json::Error),
}

/// Exponential layer configuration.
///
/// Layer K holds all intervals whose size falls in the half-open range
/// `[min_size * gf^K, min_size * gf^(K+1))`.
/// Layer 0 holds intervals with size in `[0, min_size)`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayerConfig {
    pub min_size: u32,
    pub growth_factor: u32,
}

impl LayerConfig {
    /// Assign an interval of the given size to a layer index.
    ///
    /// Layer 0: size in [0, min_size)
    /// Layer K: size in [min_size * gf^(K-1), min_size * gf^K)  for K >= 1
    pub fn layer_for_size(&self, size: u32) -> usize {
        if size < self.min_size {
            return 0;
        }
        // Find K such that min_size * gf^(K-1) <= size < min_size * gf^K
        // i.e. K = ceil(log_{gf}(size / min_size)) but with integer arithmetic
        let mut threshold = self.min_size as u64;
        let gf = self.growth_factor as u64;
        let mut layer = 0usize;
        while threshold <= size as u64 {
            threshold = threshold.saturating_mul(gf);
            layer += 1;
        }
        layer
    }

    /// The exclusive upper bound on interval sizes in layer K.
    ///
    /// Any interval in layer K has size < layer_max_size(K).
    /// This is also the fast-forward distance used by the sweep query.
    pub fn layer_max_size(&self, layer: usize) -> u32 {
        // layer 0: [0, min_size)  → max = min_size
        // layer K: max = min_size * gf^K
        let mut val = self.min_size as u64;
        let gf = self.growth_factor as u64;
        for _ in 0..layer {
            val = val.saturating_mul(gf);
        }
        val.min(u32::MAX as u64) as u32
    }

    /// Default config suitable for genomic data (~15 layers covers u32 range).
    pub fn default_genomic() -> Self {
        Self {
            min_size: 128,
            growth_factor: 4,
        }
    }
}

/// Metadata for a single indexed source file (Database SID entry).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SidEntry {
    /// File name or human-readable identifier
    pub name: String,
}

/// Global database metadata stored as `{db}/meta.json`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Meta {
    /// Map from D_SID (u32) to source metadata
    pub sid_map: HashMap<u32, SidEntry>,
    /// Layer partitioning rules
    pub layer_config: LayerConfig,
    /// All shards (coordinate spaces) present in the database
    pub shards: Vec<String>,
    /// Number of layers written (highest layer index used + 1)
    pub num_layers: usize,
}

impl Meta {
    /// Create a new, empty `Meta` with the given layer configuration.
    pub fn new(layer_config: LayerConfig) -> Self {
        Self {
            sid_map: HashMap::new(),
            layer_config,
            shards: Vec::new(),
            num_layers: 0,
        }
    }

    /// Load `meta.json` from a database directory.
    pub fn load(db_path: &Path) -> Result<Self, MetaError> {
        let content = fs::read_to_string(db_path.join("meta.json"))?;
        Ok(serde_json::from_str(&content)?)
    }

    /// Save `meta.json` to a database directory (creates or overwrites).
    pub fn save(&self, db_path: &Path) -> Result<(), MetaError> {
        let json = serde_json::to_string_pretty(self)?;
        fs::write(db_path.join("meta.json"), json)?;
        Ok(())
    }

    /// The next D_SID to assign (one past the current maximum, or 0 if empty).
    pub fn next_sid(&self) -> u32 {
        self.sid_map.keys().copied().max().map(|k| k + 1).unwrap_or(0)
    }

    /// Number of distinct indexed source files.
    pub fn num_sources(&self) -> usize {
        self.sid_map.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_layer_for_size_defaults() {
        let cfg = LayerConfig::default_genomic(); // min=128, gf=4
        assert_eq!(cfg.layer_for_size(0), 0);
        assert_eq!(cfg.layer_for_size(127), 0);
        assert_eq!(cfg.layer_for_size(128), 1);
        assert_eq!(cfg.layer_for_size(511), 1);
        assert_eq!(cfg.layer_for_size(512), 2);
        assert_eq!(cfg.layer_for_size(2047), 2);
        assert_eq!(cfg.layer_for_size(2048), 3);
    }

    #[test]
    fn test_layer_max_size() {
        let cfg = LayerConfig::default_genomic();
        assert_eq!(cfg.layer_max_size(0), 128);  // [0, 128)
        assert_eq!(cfg.layer_max_size(1), 512);  // [128, 512)
        assert_eq!(cfg.layer_max_size(2), 2048); // [512, 2048)
    }

    #[test]
    fn test_meta_next_sid() {
        let mut meta = Meta::new(LayerConfig::default_genomic());
        assert_eq!(meta.next_sid(), 0);
        meta.sid_map.insert(0, SidEntry { name: "a.bed".into() });
        assert_eq!(meta.next_sid(), 1);
        meta.sid_map.insert(5, SidEntry { name: "b.bed".into() });
        assert_eq!(meta.next_sid(), 6);
    }

    #[test]
    fn test_meta_save_load_roundtrip() {
        let dir = TempDir::new().unwrap();
        let mut meta = Meta::new(LayerConfig { min_size: 256, growth_factor: 2 });
        meta.sid_map.insert(0, SidEntry { name: "test.bed".into() });
        meta.shards = vec!["chr1".into(), "chr2".into()];
        meta.num_layers = 3;

        meta.save(dir.path()).unwrap();

        let loaded = Meta::load(dir.path()).unwrap();
        assert_eq!(loaded.sid_map.len(), 1);
        assert_eq!(loaded.shards, vec!["chr1", "chr2"]);
        assert_eq!(loaded.num_layers, 3);
        assert_eq!(loaded.layer_config.min_size, 256);
    }
}
