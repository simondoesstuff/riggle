use std::fs::{File, OpenOptions};
use std::mem::size_of;
use std::path::Path;

use memmap2::{Mmap, MmapMut};
use thiserror::Error;

use crate::core::Interval;

/// Errors from layer file operations
#[derive(Debug, Error)]
pub enum LayerError {
    #[error("IO error: {0}")]
    Io(std::io::Error),
    #[error("Position file has invalid size: {0} bytes")]
    InvalidPositionSize(u64),
    #[error("Sid file has invalid size: {0} bytes")]
    InvalidSidSize(u64),
    #[error("Position count {pos} does not match sid count {sid}")]
    SidSizeMismatch { pos: usize, sid: usize },
    #[error("Index file has invalid size: {0} bytes (not a multiple of 4)")]
    InvalidIndexSize(u64),
}

impl From<std::io::Error> for LayerError {
    fn from(e: std::io::Error) -> Self {
        LayerError::Io(e)
    }
}

/// A (start, end) coordinate pair stored in the positions file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct Position {
    pub start: u32,
    pub end: u32,
}

/// A read-only memory-mapped layer, split into a positions file and a sids file.
///
/// The positions file is a flat sorted array of `Position` (8 bytes each).
/// The sids file is a parallel flat array of `u32` D_SIDs (4 bytes each).
/// The hot scan loop reads only positions; sids are fetched only on confirmed hits.
pub struct MappedLayer {
    pos_mmap: Mmap,
    sid_mmap: Mmap,
}

impl MappedLayer {
    /// Open a layer from its `.pos` and `.sid` files.
    pub fn open(pos_path: &Path, sid_path: &Path) -> Result<Self, LayerError> {
        let pos_file = File::open(pos_path)?;
        let sid_file = File::open(sid_path)?;

        let pos_size = pos_file.metadata()?.len();
        let sid_size = sid_file.metadata()?.len();

        if pos_size % 8 != 0 {
            return Err(LayerError::InvalidPositionSize(pos_size));
        }
        if sid_size % 4 != 0 {
            return Err(LayerError::InvalidSidSize(sid_size));
        }

        let pos_count = (pos_size / 8) as usize;
        let sid_count = (sid_size / 4) as usize;
        if pos_count != sid_count {
            return Err(LayerError::SidSizeMismatch { pos: pos_count, sid: sid_count });
        }

        let pos_mmap = unsafe { Mmap::map(&pos_file)? };
        let sid_mmap = unsafe { Mmap::map(&sid_file)? };

        // Positions are always scanned sequentially; sids are accessed sparsely.
        let _ = pos_mmap.advise(memmap2::Advice::Sequential);

        Ok(Self { pos_mmap, sid_mmap })
    }

    /// Sorted positions slice (zero-copy).
    #[inline]
    pub fn positions(&self) -> &[Position] {
        let count = self.pos_mmap.len() / size_of::<Position>();
        unsafe { std::slice::from_raw_parts(self.pos_mmap.as_ptr() as *const Position, count) }
    }

    /// Parallel sids slice (zero-copy).
    #[inline]
    pub fn sids(&self) -> &[u32] {
        let count = self.sid_mmap.len() / size_of::<u32>();
        unsafe { std::slice::from_raw_parts(self.sid_mmap.as_ptr() as *const u32, count) }
    }

    /// Number of intervals in this layer.
    #[inline]
    pub fn len(&self) -> usize {
        self.pos_mmap.len() / size_of::<Position>()
    }

    /// True if the layer contains no intervals.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.pos_mmap.is_empty()
    }
}

/// Write a sorted slice of intervals to new `.pos` and `.sid` files.
pub fn write_layer(pos_path: &Path, sid_path: &Path, intervals: &[Interval]) -> Result<(), LayerError> {
    let positions: Vec<Position> = intervals.iter().map(|iv| Position { start: iv.start, end: iv.end }).collect();
    let sids: Vec<u32> = intervals.iter().map(|iv| iv.sid).collect();

    let pos_bytes: &[u8] = unsafe {
        std::slice::from_raw_parts(
            positions.as_ptr() as *const u8,
            positions.len() * size_of::<Position>(),
        )
    };
    let sid_bytes: &[u8] = unsafe {
        std::slice::from_raw_parts(
            sids.as_ptr() as *const u8,
            sids.len() * size_of::<u32>(),
        )
    };

    std::fs::write(pos_path, pos_bytes)?;
    std::fs::write(sid_path, sid_bytes)?;
    Ok(())
}

/// Merge `new_sorted` into existing `.pos` and `.sid` files using the
/// Extend & Reverse-Shift algorithm applied simultaneously to both files.
///
/// The invariant `out_ptr = old_ptr + new_ptr ≥ old_ptr` holds on both arrays
/// independently, so the write index never aliases the old-read index in either file.
pub fn extend_layer(pos_path: &Path, sid_path: &Path, new_sorted: &[Interval]) -> Result<(), LayerError> {
    if new_sorted.is_empty() {
        return Ok(());
    }

    let pos_file = OpenOptions::new().read(true).write(true).open(pos_path)?;
    let sid_file = OpenOptions::new().read(true).write(true).open(sid_path)?;

    let old_pos_byte_len = pos_file.metadata()?.len();
    let old_sid_byte_len = sid_file.metadata()?.len();

    if old_pos_byte_len % 8 != 0 {
        return Err(LayerError::InvalidPositionSize(old_pos_byte_len));
    }
    if old_sid_byte_len % 4 != 0 {
        return Err(LayerError::InvalidSidSize(old_sid_byte_len));
    }

    let old_count = (old_pos_byte_len / 8) as usize;
    let new_count = new_sorted.len();
    let total_count = old_count + new_count;

    pos_file.set_len(total_count as u64 * 8)?;
    sid_file.set_len(total_count as u64 * 4)?;

    let mut pos_mmap = unsafe { MmapMut::map_mut(&pos_file)? };
    let mut sid_mmap = unsafe { MmapMut::map_mut(&sid_file)? };

    let pos_data: &mut [Position] = unsafe {
        std::slice::from_raw_parts_mut(pos_mmap.as_mut_ptr() as *mut Position, total_count)
    };
    let sid_data: &mut [u32] = unsafe {
        std::slice::from_raw_parts_mut(sid_mmap.as_mut_ptr() as *mut u32, total_count)
    };

    let mut old_ptr = old_count;
    let mut new_ptr = new_count;
    let mut out_ptr = total_count;

    while old_ptr > 0 && new_ptr > 0 {
        out_ptr -= 1;
        if pos_data[old_ptr - 1].start >= new_sorted[new_ptr - 1].start {
            pos_data[out_ptr] = pos_data[old_ptr - 1];
            sid_data[out_ptr] = sid_data[old_ptr - 1];
            old_ptr -= 1;
        } else {
            let niv = &new_sorted[new_ptr - 1];
            pos_data[out_ptr] = Position { start: niv.start, end: niv.end };
            sid_data[out_ptr] = niv.sid;
            new_ptr -= 1;
        }
    }

    while new_ptr > 0 {
        out_ptr -= 1;
        let niv = &new_sorted[new_ptr - 1];
        pos_data[out_ptr] = Position { start: niv.start, end: niv.end };
        sid_data[out_ptr] = niv.sid;
        new_ptr -= 1;
    }

    pos_mmap.flush()?;
    sid_mmap.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    fn iv(start: u32, end: u32, sid: u32) -> Interval {
        Interval::new(start, end, sid)
    }

    fn make_files() -> (NamedTempFile, NamedTempFile) {
        (NamedTempFile::new().unwrap(), NamedTempFile::new().unwrap())
    }

    fn open(pos: &NamedTempFile, sid: &NamedTempFile) -> MappedLayer {
        MappedLayer::open(pos.path(), sid.path()).unwrap()
    }

    #[test]
    fn test_write_and_read() {
        let (pf, sf) = make_files();
        let intervals = vec![iv(10, 20, 0), iv(30, 40, 1), iv(50, 60, 2)];
        write_layer(pf.path(), sf.path(), &intervals).unwrap();

        let mapped = open(&pf, &sf);
        assert_eq!(mapped.len(), 3);
        let pos = mapped.positions();
        let sids = mapped.sids();
        assert_eq!(pos[0], Position { start: 10, end: 20 });
        assert_eq!(pos[1], Position { start: 30, end: 40 });
        assert_eq!(pos[2], Position { start: 50, end: 60 });
        assert_eq!(sids, &[0u32, 1, 2]);
    }

    #[test]
    fn test_extend_merges_correctly() {
        let (pf, sf) = make_files();
        let existing = vec![iv(10, 20, 0), iv(30, 40, 1), iv(70, 80, 2)];
        write_layer(pf.path(), sf.path(), &existing).unwrap();

        let new_ivs = vec![iv(5, 15, 3), iv(25, 35, 4), iv(60, 65, 5)];
        extend_layer(pf.path(), sf.path(), &new_ivs).unwrap();

        let mapped = open(&pf, &sf);
        assert_eq!(mapped.len(), 6);

        let pos = mapped.positions();
        for w in pos.windows(2) {
            assert!(w[0].start <= w[1].start, "not sorted: {:?}", pos);
        }
        let starts: Vec<u32> = pos.iter().map(|p| p.start).collect();
        assert_eq!(starts, vec![5, 10, 25, 30, 60, 70]);
    }

    #[test]
    fn test_extend_into_empty() {
        let (pf, sf) = make_files();
        write_layer(pf.path(), sf.path(), &[]).unwrap();

        let new_ivs = vec![iv(1, 2, 0), iv(3, 4, 1)];
        extend_layer(pf.path(), sf.path(), &new_ivs).unwrap();

        let mapped = open(&pf, &sf);
        assert_eq!(mapped.len(), 2);
        assert_eq!(mapped.positions()[0].start, 1);
        assert_eq!(mapped.positions()[1].start, 3);
    }

    #[test]
    fn test_extend_with_empty_new() {
        let (pf, sf) = make_files();
        let existing = vec![iv(10, 20, 0), iv(30, 40, 1)];
        write_layer(pf.path(), sf.path(), &existing).unwrap();

        extend_layer(pf.path(), sf.path(), &[]).unwrap();

        let mapped = open(&pf, &sf);
        assert_eq!(mapped.len(), 2);
    }

    #[test]
    fn test_extend_all_new_before_existing() {
        let (pf, sf) = make_files();
        let existing = vec![iv(100, 200, 0), iv(300, 400, 1)];
        write_layer(pf.path(), sf.path(), &existing).unwrap();

        let new_ivs = vec![iv(1, 2, 2), iv(3, 4, 3)];
        extend_layer(pf.path(), sf.path(), &new_ivs).unwrap();

        let mapped = open(&pf, &sf);
        let starts: Vec<u32> = mapped.positions().iter().map(|p| p.start).collect();
        assert_eq!(starts, vec![1, 3, 100, 300]);
    }

    #[test]
    fn test_extend_all_new_after_existing() {
        let (pf, sf) = make_files();
        let existing = vec![iv(1, 2, 0), iv(3, 4, 1)];
        write_layer(pf.path(), sf.path(), &existing).unwrap();

        let new_ivs = vec![iv(100, 200, 2), iv(300, 400, 3)];
        extend_layer(pf.path(), sf.path(), &new_ivs).unwrap();

        let mapped = open(&pf, &sf);
        let starts: Vec<u32> = mapped.positions().iter().map(|p| p.start).collect();
        assert_eq!(starts, vec![1, 3, 100, 300]);
    }

    #[test]
    fn test_sids_travel_with_positions() {
        let (pf, sf) = make_files();
        write_layer(pf.path(), sf.path(), &[iv(10, 20, 7), iv(30, 40, 3)]).unwrap();
        extend_layer(pf.path(), sf.path(), &[iv(5, 15, 99), iv(25, 35, 42)]).unwrap();

        let mapped = open(&pf, &sf);
        let pos = mapped.positions();
        let sids = mapped.sids();
        // sorted order: 5, 10, 25, 30
        assert_eq!(pos[0].start, 5);  assert_eq!(sids[0], 99);
        assert_eq!(pos[1].start, 10); assert_eq!(sids[1], 7);
        assert_eq!(pos[2].start, 25); assert_eq!(sids[2], 42);
        assert_eq!(pos[3].start, 30); assert_eq!(sids[3], 3);
    }
}
