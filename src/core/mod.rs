use voracious_radix_sort::Radixable;

/// The universal interval primitive used everywhere in Riggle.
///
/// At index time, `sid` is the Database SID (D_SID) identifying the source file.
/// At query time, `sid` is the Query SID (Q_SID) identifying the query file.
///
/// `repr(C)` guarantees the field layout matches what is written to / read from
/// memory-mapped layer files, enabling zero-copy casting between `&[u8]` and
/// `&[Interval]`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
#[repr(C)]
pub struct Interval {
    pub start: u32,
    pub end: u32,
    pub sid: u32,
}

impl Interval {
    #[inline]
    pub fn new(start: u32, end: u32, sid: u32) -> Self {
        Self { start, end, sid }
    }

    #[inline]
    pub fn size(&self) -> u32 {
        self.end - self.start
    }
}

// ---------------------------------------------------------------------------
// voracious_radix_sort integration
// Sorts by `start` coordinate, which is the canonical order for layer files.
// ---------------------------------------------------------------------------

impl Radixable<u32> for Interval {
    type Key = u32;

    #[inline]
    fn key(&self) -> u32 {
        self.start
    }
}
