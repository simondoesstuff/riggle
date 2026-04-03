use rkyv::{Archive, Deserialize, Serialize};

/// Basic interval representation with half-open coordinates [start, end)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Archive, Serialize, Deserialize)]
pub struct Interval {
    pub start: u32,
    pub end: u32,
}

impl Interval {
    pub fn new(start: u32, end: u32) -> Self {
        debug_assert!(start <= end, "Invalid interval: start > end");
        Self { start, end }
    }

    /// Check if this interval overlaps with another (half-open intervals)
    #[inline]
    pub fn overlaps(&self, other: &Interval) -> bool {
        self.start < other.end && other.start < self.end
    }

    /// Check if this interval fully contains another
    #[inline]
    pub fn contains(&self, other: &Interval) -> bool {
        self.start <= other.start && other.end <= self.end
    }

    /// Get the length of the interval
    #[inline]
    pub fn len(&self) -> u32 {
        self.end - self.start
    }

    /// Check if interval is empty
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.start == self.end
    }
}

/// Interval tagged with a source ID
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Archive, Serialize, Deserialize)]
pub struct TaggedInterval {
    pub iv: Interval,
    pub sid: u32,
}

impl TaggedInterval {
    pub fn new(start: u32, end: u32, sid: u32) -> Self {
        Self {
            iv: Interval::new(start, end),
            sid,
        }
    }
}

/// Exponential layer identifier (log2 of interval size range)
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Archive, Serialize, Deserialize,
)]
pub struct LayerID(pub u8);

impl LayerID {
    /// Compute layer ID from interval length
    pub fn from_length(len: u32) -> Self {
        if len == 0 {
            LayerID(0)
        } else {
            // floor(log2(x))
            LayerID((32 - len.leading_zeros() - 1) as u8)
        }
    }
}

/// Fixed genomic territory identifier
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Archive, Serialize, Deserialize,
)]
pub struct ChunkID(pub u32);

impl ChunkID {
    /// Compute chunk ID from coordinate and chunk size
    pub fn from_coord(coord: u32, chunk_size: u32) -> Self {
        ChunkID(coord / chunk_size)
    }
}

/// Tile subdivision within a chunk
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Archive, Serialize, Deserialize,
)]
pub struct TileID(pub u32);

impl TileID {
    /// Compute tile ID from coordinate and tile size
    pub fn from_coord(coord: u32, tile_size: u32) -> Self {
        TileID(coord / tile_size)
    }
}

/// Tile structure - the leaf node in the index
/// Contains intervals that touch this tile, organized for efficient querying.
///
/// INVARIANT: `start_ivs` and `end_ivs` are sorted by offset (first tuple element)
/// to enable binary search during query. This sorting is performed by `index_sweep`.
#[derive(Debug, Clone, Archive, Serialize, Deserialize)]
pub struct Tile {
    /// Start coordinate of this tile
    pub start_coord: u32,
    /// Running counts: intervals that span the entire tile (sid, count)
    pub running_counts: Vec<(u32, u32)>,
    /// Intervals starting within this tile (offset from tile start, sid)
    /// SORTED by offset for binary search
    /// Note: u32 offset required since tile_size can exceed 65535 for larger layers
    pub start_ivs: Vec<(u32, u32)>,
    /// Intervals ending within this tile (offset from tile start, sid)
    /// SORTED by offset for binary search
    /// Note: u32 offset required since tile_size can exceed 65535 for larger layers
    pub end_ivs: Vec<(u32, u32)>,
}

impl Tile {
    pub fn new(start_coord: u32) -> Self {
        Self {
            start_coord,
            running_counts: Vec::new(),
            start_ivs: Vec::new(),
            end_ivs: Vec::new(),
        }
    }

    /// Check if the tile is empty (no intervals touch it)
    pub fn is_empty(&self) -> bool {
        self.running_counts.is_empty() && self.start_ivs.is_empty() && self.end_ivs.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interval_overlap() {
        let a = Interval::new(10, 20);
        let b = Interval::new(15, 25);
        let c = Interval::new(20, 30);
        let d = Interval::new(5, 15);

        assert!(a.overlaps(&b)); // Partial overlap
        assert!(b.overlaps(&a)); // Symmetric
        assert!(!a.overlaps(&c)); // Adjacent, no overlap (half-open)
        assert!(a.overlaps(&d)); // Partial overlap
    }

    #[test]
    fn test_interval_contains() {
        let outer = Interval::new(10, 30);
        let inner = Interval::new(15, 25);
        let partial = Interval::new(5, 20);

        assert!(outer.contains(&inner));
        assert!(!inner.contains(&outer));
        assert!(!outer.contains(&partial));
    }

    #[test]
    fn test_interval_length() {
        assert_eq!(Interval::new(10, 20).len(), 10);
        assert_eq!(Interval::new(0, 100).len(), 100);
        assert_eq!(Interval::new(50, 50).len(), 0);
    }

    #[test]
    fn test_tagged_interval_creation() {
        let ti = TaggedInterval::new(100, 200, 42);
        assert_eq!(ti.iv.start, 100);
        assert_eq!(ti.iv.end, 200);
        assert_eq!(ti.sid, 42);
    }

    #[test]
    fn test_layer_id_from_length() {
        assert_eq!(LayerID::from_length(1), LayerID(0));
        assert_eq!(LayerID::from_length(2), LayerID(1));
        assert_eq!(LayerID::from_length(3), LayerID(1));
        assert_eq!(LayerID::from_length(4), LayerID(2));
        assert_eq!(LayerID::from_length(1024), LayerID(10));
    }

    #[test]
    fn test_chunk_id_from_coord() {
        let chunk_size = 1000;
        assert_eq!(ChunkID::from_coord(0, chunk_size), ChunkID(0));
        assert_eq!(ChunkID::from_coord(500, chunk_size), ChunkID(0));
        assert_eq!(ChunkID::from_coord(1000, chunk_size), ChunkID(1));
        assert_eq!(ChunkID::from_coord(2500, chunk_size), ChunkID(2));
    }

    #[test]
    fn test_tile_id_from_coord() {
        let tile_size = 100;
        assert_eq!(TileID::from_coord(0, tile_size), TileID(0));
        assert_eq!(TileID::from_coord(50, tile_size), TileID(0));
        assert_eq!(TileID::from_coord(100, tile_size), TileID(1));
    }

    #[test]
    fn test_tile_empty() {
        let empty = Tile::new(0);
        assert!(empty.is_empty());

        let mut non_empty = Tile::new(0);
        non_empty.running_counts.push((1, 5));
        assert!(!non_empty.is_empty());
    }

    #[test]
    fn test_tile_serialization_roundtrip() {
        use rkyv::rancor::Error;

        let mut tile = Tile::new(1000);
        tile.running_counts.push((1, 10));
        tile.running_counts.push((2, 5));
        tile.start_ivs.push((50, 3));
        tile.end_ivs.push((75, 4));

        // Serialize
        let bytes = rkyv::to_bytes::<Error>(&tile).unwrap();

        // Deserialize
        let archived = rkyv::access::<ArchivedTile, Error>(&bytes).unwrap();

        assert_eq!(archived.start_coord, 1000);
        assert_eq!(archived.running_counts.len(), 2);
        assert_eq!(archived.start_ivs.len(), 1);
        assert_eq!(archived.end_ivs.len(), 1);
    }
}
