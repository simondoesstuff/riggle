//! Single-file rkyv-backed store for per-source FFT moment tables.
//!
//! Parallel to `depthmap.rkyv`; written at index time, read at query time.
//! `maps[sid]` holds the per-chromosome compacted (mean, var) tables for that
//! source.  SIDs are dense, so the vec index is the SID — O(1) lookup.

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

    /// Return the moment tables for `sid`, or `None` if not present.
    pub fn get(&self, sid: u32) -> Option<Vec<ChromMoments>> {
        let entry = self.archived().maps.get(sid as usize)?;
        if entry.is_empty() {
            return None;
        }
        let chroms = entry
            .iter()
            .map(|c| {
                let n_bins = u32::from(c.n_bins) as usize;
                let eps = f32::from(c.eps) as f64;
                let dense_max = u32::from(c.dense_max) as usize;
                let moments: Vec<(f64, f64)> = c
                    .moments
                    .chunks_exact(2)
                    .map(|pair| (f32::from(pair[0]) as f64, f32::from(pair[1]) as f64))
                    .collect();
                ChromMoments {
                    chrom: c.chrom.as_str().to_string(),
                    n_bins,
                    eps,
                    dense_max,
                    moments,
                }
            })
            .collect();
        Some(chroms)
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
            let got = store.get(sid).unwrap();
            assert_eq!(got.len(), moments.len());
            for (a, b) in moments.iter().zip(got.iter()) {
                assert_eq!(a.chrom, b.chrom);
                assert_eq!(a.n_bins, b.n_bins);
                assert_eq!(a.moments.len(), b.moments.len());
                // Check round-trip precision (f32 store, f64 in-memory)
                for (&(am, av), &(bm, bv)) in a.moments.iter().zip(b.moments.iter()) {
                    assert!((am - bm).abs() < am.abs() * 1e-5 + 1e-10,
                        "mean mismatch: {} vs {}", am, bm);
                    assert!((av - bv).abs() < av.abs() * 1e-5 + 1e-10,
                        "var mismatch: {} vs {}", av, bv);
                }
            }
        }

        assert!(store.get(2).is_none());
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
        assert!(store.get(0).is_some());
        assert!(store.get(1).is_some());
        assert!(store.get(2).is_none());
    }
}
