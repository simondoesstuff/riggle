//! Single-file rkyv-backed store for per-source FFT moment tables.
//!
//! Parallel to `depthmap.rkyv`; written at index time, read at query time.
//! `maps[sid]` holds the per-chromosome compacted (mean, var) tables for that
//! source.  SIDs are dense, so the vec index is the SID — O(1) lookup.

use std::collections::HashMap;
use std::fs;
use std::io;
use std::path::Path;

use memmap2::Mmap;
use rkyv::{Archive, Deserialize, Serialize};

use crate::fourier::ChromMoments;

// ── Archived types ────────────────────────────────────────────────────────────

#[derive(Archive, Serialize, Deserialize, Clone)]
struct StoredChromMoments {
    chrom: String,
    n_bins: u32,
    /// Stored as f32 to save space; converted to f64 on load.
    eps: f32,
    dense_max: u32,
    /// Flattened (mean, var) pairs: [mean_0, var_0, mean_1, var_1, …]
    moments: Vec<f32>,
}

#[derive(Archive, Serialize, Deserialize, Clone)]
struct MomentStore {
    /// `maps[sid]` = per-chromosome moments for that source; empty = absent.
    maps: Vec<Vec<StoredChromMoments>>,
}

// ── Zerocopy handle for one SID ───────────────────────────────────────────────

/// Zerocopy moment lookup for a single DB source.
///
/// Backed directly by the mmap'd `momentmap.rkyv` — no f32→f64 conversion or
/// allocation beyond the one-time O(chroms) index built at access time.
pub struct MappedSidMoments<'a> {
    entry: &'a [ArchivedStoredChromMoments],
    /// chrom name → index into `entry`; built once, O(1) lookup.
    index: HashMap<&'a str, usize>,
}

impl<'a> MappedSidMoments<'a> {
    /// Return `(mean, var)` for a sliding window of `l_bins` bins on `chrom`.
    ///
    /// Reads directly from the archived f32 slice — no allocation.
    /// Returns `None` when `chrom` is absent, `l_bins` rounds to 0, or ≥ n_bins.
    pub fn lookup(&self, chrom: &str, l_bins: f64) -> Option<(f64, f64)> {
        let &i = self.index.get(chrom)?;
        let acm = &self.entry[i];
        let n_bins = u32::from(acm.n_bins) as usize;
        let dense_max = u32::from(acm.dense_max) as usize;
        let eps = f32::from(acm.eps) as f64;
        let l = l_bins.round() as usize;
        if l == 0 || l >= n_bins || acm.moments.is_empty() {
            return None;
        }
        let n_pairs = acm.moments.len() / 2;
        let idx = if l <= dense_max {
            l - 1
        } else {
            let m = ((l as f64 / dense_max as f64).ln() / (1.0 + eps).ln())
                .floor()
                .max(1.0) as usize;
            (dense_max + m - 1).min(n_pairs - 1)
        };
        let mean = f32::from(acm.moments[idx * 2]) as f64;
        let var = f32::from(acm.moments[idx * 2 + 1]) as f64;
        Some((mean, var))
    }
}

// ── Read-only memmap handle ───────────────────────────────────────────────────

pub struct MappedMomentStore {
    mmap: Mmap,
}

impl MappedMomentStore {
    pub fn open(path: &Path) -> io::Result<Self> {
        let file = fs::File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Ok(Self { mmap })
    }

    fn archived(&self) -> &ArchivedMomentStore {
        unsafe { rkyv::access_unchecked::<ArchivedMomentStore>(&self.mmap[..]) }
    }

    /// Return a zerocopy handle for `sid`, or `None` if not present.
    ///
    /// Builds a chrom→index map (O(chroms)) once; subsequent `lookup` calls are O(1).
    pub fn get_sid(&self, sid: u32) -> Option<MappedSidMoments<'_>> {
        let entry = self.archived().maps.get(sid as usize)?;
        if entry.is_empty() {
            return None;
        }
        let index: HashMap<&str, usize> = entry
            .iter()
            .enumerate()
            .map(|(i, c)| (c.chrom.as_str(), i))
            .collect();
        Some(MappedSidMoments { entry: &entry[..], index })
    }
}

// ── Mutable builder ───────────────────────────────────────────────────────────

pub struct MomentStoreBuilder {
    maps: Vec<Vec<StoredChromMoments>>,
}

impl Default for MomentStoreBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl MomentStoreBuilder {
    pub fn new() -> Self {
        Self { maps: Vec::new() }
    }

    pub fn load(path: &Path) -> io::Result<Self> {
        let data = fs::read(path)?;
        let archived = unsafe { rkyv::access_unchecked::<ArchivedMomentStore>(&data) };
        let store: MomentStore =
            rkyv::deserialize::<MomentStore, rkyv::rancor::Error>(archived)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
        Ok(Self { maps: store.maps })
    }

    pub fn insert(&mut self, sid: u32, moments: &[ChromMoments]) {
        let idx = sid as usize;
        if idx >= self.maps.len() {
            self.maps.resize_with(idx + 1, Vec::new);
        }
        self.maps[idx] = moments
            .iter()
            .map(|m| {
                let mut flat = Vec::with_capacity(m.moments.len() * 2);
                for &(mean, var) in &m.moments {
                    flat.push(mean as f32);
                    flat.push(var as f32);
                }
                StoredChromMoments {
                    chrom: m.chrom.clone(),
                    n_bins: m.n_bins as u32,
                    eps: m.eps as f32,
                    dense_max: m.dense_max as u32,
                    moments: flat,
                }
            })
            .collect();
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        let store = MomentStore {
            maps: self.maps.clone(),
        };
        let bytes = rkyv::to_bytes::<rkyv::rancor::Error>(&store)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
        fs::write(path, bytes.as_slice())
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fourier::{BedMap, DepthMap, build_depth_moments, DEFAULT_MOMENTS_EPS};
    use tempfile::TempDir;

    fn make_moments() -> Vec<ChromMoments> {
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(10_000_000, 20_000_000)]);
        bed.insert("chr1".to_string(), vec![(1_000_000, 2_000_000)]);
        let dm = DepthMap::build(&bed);
        build_depth_moments(&dm, DEFAULT_MOMENTS_EPS)
    }

    #[test]
    fn test_moment_store_roundtrip() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("momentmap.rkyv");

        let moments = make_moments();
        let mut builder = MomentStoreBuilder::new();
        builder.insert(0, &moments);
        builder.insert(1, &moments);
        builder.save(&path).unwrap();

        let store = MappedMomentStore::open(&path).unwrap();

        for sid in 0..2u32 {
            let sid_m = store.get_sid(sid).unwrap();
            for orig in &moments {
                // Verify dense and sparse lookups round-trip within f32 precision.
                for &l in &[50.0f64, 500.0] {
                    if let Some((exp_m, exp_v)) = orig.lookup(l) {
                        let (got_m, got_v) = sid_m.lookup(&orig.chrom, l).unwrap();
                        assert!(
                            (got_m - exp_m).abs() < exp_m.abs() * 1e-5 + 1e-10,
                            "mean mismatch L={}: {} vs {}", l, got_m, exp_m
                        );
                        assert!(
                            (got_v - exp_v).abs() < exp_v.abs() * 1e-5 + 1e-10,
                            "var mismatch L={}: {} vs {}", l, got_v, exp_v
                        );
                    }
                }
            }
        }

        assert!(store.get_sid(2).is_none());
    }

    #[test]
    fn test_moment_store_incremental() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("momentmap.rkyv");

        let moments = make_moments();

        let mut b1 = MomentStoreBuilder::new();
        b1.insert(0, &moments);
        b1.save(&path).unwrap();

        let mut b2 = MomentStoreBuilder::load(&path).unwrap();
        b2.insert(1, &moments);
        b2.save(&path).unwrap();

        let store = MappedMomentStore::open(&path).unwrap();
        assert!(store.get_sid(0).is_some());
        assert!(store.get_sid(1).is_some());
        assert!(store.get_sid(2).is_none());
    }

    #[test]
    fn test_lookup_absent_chrom_returns_none() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("momentmap.rkyv");

        let moments = make_moments();
        let mut builder = MomentStoreBuilder::new();
        builder.insert(0, &moments);
        builder.save(&path).unwrap();

        let store = MappedMomentStore::open(&path).unwrap();
        let sid_m = store.get_sid(0).unwrap();
        assert!(sid_m.lookup("chrXYZ", 100.0).is_none());
    }
}
