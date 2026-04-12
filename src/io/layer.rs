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
    #[error("Layer file has invalid size: {0} bytes (not a multiple of interval stride)")]
    InvalidSize(u64),
    #[error("Index file has invalid size: {0} bytes (not a multiple of 4)")]
    InvalidIndexSize(u64),
}

impl From<std::io::Error> for LayerError {
    fn from(e: std::io::Error) -> Self {
        LayerError::Io(e)
    }
}

/// A read-only memory-mapped layer file.
///
/// The file is a flat, tightly packed, sorted array of `Interval` tuples.
/// Zero-copy access via `intervals()` casts the mmap bytes directly.
pub struct MappedLayer {
    mmap: Mmap,
}

impl MappedLayer {
    /// Open a layer `.bin` file for zero-copy reading.
    pub fn open(path: &Path) -> Result<Self, LayerError> {
        let file = File::open(path)?;
        let size = file.metadata()?.len();
        let stride = size_of::<Interval>() as u64;
        if size % stride != 0 {
            return Err(LayerError::InvalidSize(size));
        }
        let mmap = unsafe { Mmap::map(&file)? };
        Ok(Self { mmap })
    }

    /// View the layer as a sorted slice of intervals (zero-copy).
    #[inline]
    pub fn intervals(&self) -> &[Interval] {
        let stride = size_of::<Interval>();
        let count = self.mmap.len() / stride;
        // SAFETY: the file was written by `write_layer`/`extend_layer` as
        // packed `repr(C)` Interval structs, so the bytes are valid.
        unsafe { std::slice::from_raw_parts(self.mmap.as_ptr() as *const Interval, count) }
    }

    /// Number of intervals in this layer.
    #[inline]
    pub fn len(&self) -> usize {
        self.mmap.len() / size_of::<Interval>()
    }

    /// True if the layer contains no intervals.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.mmap.is_empty()
    }
}

/// Write a sorted slice of intervals to a new layer file.
///
/// The file is created (or truncated) and filled with the raw bytes of
/// the intervals. The caller must ensure `intervals` is sorted by `start`.
pub fn write_layer(path: &Path, intervals: &[Interval]) -> Result<(), LayerError> {
    let bytes: &[u8] = unsafe {
        std::slice::from_raw_parts(
            intervals.as_ptr() as *const u8,
            intervals.len() * size_of::<Interval>(),
        )
    };
    std::fs::write(path, bytes)?;
    Ok(())
}

/// Merge `new_sorted` into an existing layer file using the **Extend & Reverse-Shift**
/// algorithm so the result is sorted by `start` without loading the full file into RAM.
///
/// ### Algorithm
///
/// 1. Extend the file on disk by `new_sorted.len() * sizeof(Interval)` bytes.
/// 2. Map the whole extended file as a mutable slice of `Interval`.
/// 3. Perform a two-pointer reverse merge starting from the back:
///    - `old_ptr` scans backwards through the original `[0..old_count]` region.
///    - `new_ptr` scans backwards through `new_sorted`.
///    - The larger interval at each step is written to `out_ptr` (scanning back
///      from the absolute end of the file).
///    - Because `out_ptr >= old_ptr` at all times (the gap equals the number of
///      `new_ptr` steps taken so far), no unread data is ever overwritten.
/// 4. Flush and unmap.
pub fn extend_layer(path: &Path, new_sorted: &[Interval]) -> Result<(), LayerError> {
    if new_sorted.is_empty() {
        return Ok(());
    }

    let file = OpenOptions::new().read(true).write(true).open(path)?;
    let old_byte_len = file.metadata()?.len();
    let stride = size_of::<Interval>() as u64;

    if old_byte_len % stride != 0 {
        return Err(LayerError::InvalidSize(old_byte_len));
    }

    let old_count = (old_byte_len / stride) as usize;
    let new_count = new_sorted.len();
    let total_count = old_count + new_count;

    // Extend the file to make room for the incoming intervals.
    file.set_len(total_count as u64 * stride)?;

    let mut mmap = unsafe { MmapMut::map_mut(&file)? };
    let data: &mut [Interval] = unsafe {
        std::slice::from_raw_parts_mut(mmap.as_mut_ptr() as *mut Interval, total_count)
    };

    // Two-pointer reverse merge in-place.
    let mut old_ptr = old_count; // exclusive end of the original region
    let mut new_ptr = new_count; // exclusive end of new_sorted
    let mut out_ptr = total_count; // exclusive end of the output

    while old_ptr > 0 && new_ptr > 0 {
        out_ptr -= 1;
        // Compare by start; ties keep the existing element first (stable).
        if data[old_ptr - 1].start >= new_sorted[new_ptr - 1].start {
            data[out_ptr] = data[old_ptr - 1];
            old_ptr -= 1;
        } else {
            data[out_ptr] = new_sorted[new_ptr - 1];
            new_ptr -= 1;
        }
    }

    // Drain any remaining new_sorted elements into the front.
    while new_ptr > 0 {
        out_ptr -= 1;
        data[out_ptr] = new_sorted[new_ptr - 1];
        new_ptr -= 1;
    }
    // If old_ptr > 0 the remaining existing elements are already in place.

    mmap.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    fn iv(start: u32, end: u32, sid: u32) -> Interval {
        Interval::new(start, end, sid)
    }

    #[test]
    fn test_write_and_read() {
        let file = NamedTempFile::new().unwrap();
        let intervals = vec![iv(10, 20, 0), iv(30, 40, 1), iv(50, 60, 2)];
        write_layer(file.path(), &intervals).unwrap();

        let mapped = MappedLayer::open(file.path()).unwrap();
        assert_eq!(mapped.len(), 3);
        assert_eq!(mapped.intervals(), intervals.as_slice());
    }

    #[test]
    fn test_extend_merges_correctly() {
        let file = NamedTempFile::new().unwrap();
        let existing = vec![iv(10, 20, 0), iv(30, 40, 1), iv(70, 80, 2)];
        write_layer(file.path(), &existing).unwrap();

        let new_ivs = vec![iv(5, 15, 3), iv(25, 35, 4), iv(60, 65, 5)];
        extend_layer(file.path(), &new_ivs).unwrap();

        let mapped = MappedLayer::open(file.path()).unwrap();
        let result = mapped.intervals();
        assert_eq!(result.len(), 6);

        // Verify sorted by start
        for w in result.windows(2) {
            assert!(w[0].start <= w[1].start, "not sorted: {:?}", result);
        }
        // All original starts present
        let starts: Vec<u32> = result.iter().map(|i| i.start).collect();
        assert_eq!(starts, vec![5, 10, 25, 30, 60, 70]);
    }

    #[test]
    fn test_extend_into_empty() {
        let file = NamedTempFile::new().unwrap();
        write_layer(file.path(), &[]).unwrap();

        let new_ivs = vec![iv(1, 2, 0), iv(3, 4, 1)];
        extend_layer(file.path(), &new_ivs).unwrap();

        let mapped = MappedLayer::open(file.path()).unwrap();
        assert_eq!(mapped.len(), 2);
        assert_eq!(mapped.intervals()[0].start, 1);
        assert_eq!(mapped.intervals()[1].start, 3);
    }

    #[test]
    fn test_extend_with_empty_new() {
        let file = NamedTempFile::new().unwrap();
        let existing = vec![iv(10, 20, 0), iv(30, 40, 1)];
        write_layer(file.path(), &existing).unwrap();

        extend_layer(file.path(), &[]).unwrap();

        let mapped = MappedLayer::open(file.path()).unwrap();
        assert_eq!(mapped.len(), 2);
    }

    #[test]
    fn test_extend_all_new_before_existing() {
        let file = NamedTempFile::new().unwrap();
        let existing = vec![iv(100, 200, 0), iv(300, 400, 1)];
        write_layer(file.path(), &existing).unwrap();

        let new_ivs = vec![iv(1, 2, 2), iv(3, 4, 3)];
        extend_layer(file.path(), &new_ivs).unwrap();

        let mapped = MappedLayer::open(file.path()).unwrap();
        let starts: Vec<u32> = mapped.intervals().iter().map(|i| i.start).collect();
        assert_eq!(starts, vec![1, 3, 100, 300]);
    }

    #[test]
    fn test_extend_all_new_after_existing() {
        let file = NamedTempFile::new().unwrap();
        let existing = vec![iv(1, 2, 0), iv(3, 4, 1)];
        write_layer(file.path(), &existing).unwrap();

        let new_ivs = vec![iv(100, 200, 2), iv(300, 400, 3)];
        extend_layer(file.path(), &new_ivs).unwrap();

        let mapped = MappedLayer::open(file.path()).unwrap();
        let starts: Vec<u32> = mapped.intervals().iter().map(|i| i.start).collect();
        assert_eq!(starts, vec![1, 3, 100, 300]);
    }
}
