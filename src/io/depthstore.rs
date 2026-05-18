//! Single-file rkyv-backed store for all per-source depth maps.
//!
//! Replaces the `{db}/depthmap/{sid}.bin` per-file directory with a single
//! memmap-backed archive at `{db}/depthmap.rkyv`.
//!
//! # On-disk layout
//! A single rkyv archive: `maps[sid]` holds the chromosome impulse data for
//! that source.  SIDs are dense (0, 1, …, N−1), so the vec index *is* the SID
//! — O(1) lookup with no auxiliary index.
//!
//! # Read path (`MappedDepthStore`)
//! The file is opened once per query run and kept memory-mapped.  `get(sid)`
//! is a single slice index into the archived vec — no binary search, no file
//! seek.
//!
//! # Write path (`DepthStoreBuilder`)
//! At build time, a `DepthStoreBuilder` loads the existing archive (if any),
//! appends new entries, and re-serialises the whole store.

use std::fs;
use std::io;
use std::path::Path;

use memmap2::Mmap;
use rkyv::{Archive, Deserialize, Serialize};

use crate::fourier::{ChromDepthMap, DepthMap};

// ── Archived types ────────────────────────────────────────────────────────────

#[derive(Archive, Serialize, Deserialize, Clone)]
struct StoredChrom {
    chrom: String,
    /// Number of 100 bp bins on this chromosome (fits in u32 for hg38).
    n_bins: u32,
    pos_spikes: Vec<u32>,
    neg_spikes: Vec<u32>,
}

/// Root of the rkyv archive.  `maps[sid]` is the chrom data for that source;
/// an empty inner vec means the slot is unoccupied (gap in the SID space).
#[derive(Archive, Serialize, Deserialize, Clone)]
struct DepthStore {
    maps: Vec<Vec<StoredChrom>>,
}

// ── Read-only memmap handle ───────────────────────────────────────────────────

/// Memory-mapped, read-only view over all depth maps in the database.
pub struct MappedDepthStore {
    mmap: Mmap,
}

impl MappedDepthStore {
    /// Open and memory-map the store at `path`.
    pub fn open(path: &Path) -> io::Result<Self> {
        let file = fs::File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Ok(Self { mmap })
    }

    fn archived(&self) -> &ArchivedDepthStore {
        // SAFETY: bytes were written by `DepthStoreBuilder::save` via rkyv.
        unsafe { rkyv::access_unchecked::<ArchivedDepthStore>(&self.mmap[..]) }
    }

    /// Return the `DepthMap` for `sid`, or `None` if not present.
    pub fn get(&self, sid: u32) -> Option<DepthMap> {
        let entry = self.archived().maps.get(sid as usize)?;
        if entry.is_empty() {
            return None;
        }
        let chroms = entry
            .iter()
            .map(|c| ChromDepthMap {
                chrom: c.chrom.as_str().to_string(),
                n_bins: u32::from(c.n_bins) as usize,
                pos_spikes: c.pos_spikes.iter().map(|v| u32::from(*v)).collect(),
                neg_spikes: c.neg_spikes.iter().map(|v| u32::from(*v)).collect(),
            })
            .collect();
        Some(DepthMap { chroms })
    }
}

// ── Mutable builder ───────────────────────────────────────────────────────────

/// Accumulates depth maps (one per SID) and serialises them as a single rkyv
/// archive.
pub struct DepthStoreBuilder {
    maps: Vec<Vec<StoredChrom>>,
}

impl Default for DepthStoreBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl DepthStoreBuilder {
    pub fn new() -> Self {
        Self { maps: Vec::new() }
    }

    /// Load an existing store from `path` so new entries can be appended.
    pub fn load(path: &Path) -> io::Result<Self> {
        let data = fs::read(path)?;
        let archived = unsafe { rkyv::access_unchecked::<ArchivedDepthStore>(&data) };
        let store: DepthStore =
            rkyv::deserialize::<DepthStore, rkyv::rancor::Error>(archived)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
        Ok(Self { maps: store.maps })
    }

    /// Insert a depth map for `sid`.  Extends the vec if needed; any skipped
    /// slots are left as empty vecs (SIDs are dense so this should not occur).
    pub fn insert(&mut self, sid: u32, dm: &DepthMap) {
        let idx = sid as usize;
        if idx >= self.maps.len() {
            self.maps.resize_with(idx + 1, Vec::new);
        }
        self.maps[idx] = dm
            .chroms
            .iter()
            .map(|c| StoredChrom {
                chrom: c.chrom.clone(),
                n_bins: c.n_bins as u32,
                pos_spikes: c.pos_spikes.clone(),
                neg_spikes: c.neg_spikes.clone(),
            })
            .collect();
    }

    /// Serialise to `path`, overwriting any existing file.
    pub fn save(&self, path: &Path) -> io::Result<()> {
        let store = DepthStore {
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
    use crate::fourier::{BedMap, DepthMap};
    use tempfile::TempDir;

    fn make_dm() -> DepthMap {
        let mut bed = BedMap::new();
        bed.insert("chr22".to_string(), vec![(10_000_000, 10_001_000)]);
        bed.insert("chr1".to_string(), vec![(1_000_000, 2_000_000)]);
        DepthMap::build(&bed)
    }

    #[test]
    fn test_store_roundtrip() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("depthmap.rkyv");

        let dm = make_dm();

        let mut builder = DepthStoreBuilder::new();
        builder.insert(0, &dm);
        builder.insert(1, &dm);
        builder.insert(2, &dm);
        builder.save(&path).unwrap();

        let store = MappedDepthStore::open(&path).unwrap();

        // All three SIDs round-trip correctly.
        for sid in 0..3u32 {
            let got = store.get(sid).unwrap();
            assert_eq!(got.chroms.len(), dm.chroms.len());
            for (a, b) in dm.chroms.iter().zip(got.chroms.iter()) {
                assert_eq!(a.chrom, b.chrom);
                assert_eq!(a.n_bins, b.n_bins);
                assert_eq!(a.pos_spikes, b.pos_spikes);
                assert_eq!(a.neg_spikes, b.neg_spikes);
            }
        }

        assert!(store.get(3).is_none());
    }

    #[test]
    fn test_store_incremental_load() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("depthmap.rkyv");

        let dm = make_dm();

        // First batch: SIDs 0, 1.
        let mut b1 = DepthStoreBuilder::new();
        b1.insert(0, &dm);
        b1.insert(1, &dm);
        b1.save(&path).unwrap();

        // Second batch: load existing, add SID 2.
        let mut b2 = DepthStoreBuilder::load(&path).unwrap();
        b2.insert(2, &dm);
        b2.save(&path).unwrap();

        let store = MappedDepthStore::open(&path).unwrap();
        assert!(store.get(0).is_some());
        assert!(store.get(1).is_some());
        assert!(store.get(2).is_some());
        assert!(store.get(3).is_none());
    }
}
