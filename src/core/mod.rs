use rkyv::{Archive, Deserialize, Serialize};
use voracious_radix_sort::Radixable;

// =============================================================================
// Centralized Layer Configuration Constants
// =============================================================================
// These define the hierarchical tiling structure. All layer computations
// should derive from these constants to ensure consistency.

/// Base tile size for layer 0: 2^14 = 16,384 bp (common convention)
pub const BASE_TILE_SIZE_LOG2: u32 = 14;

/// Ratio of tile_size to max_interval_size: 2^3 = 8
/// Ensures no interval can span an entire tile, eliminating running_counts
pub const TILE_MAX_INTERVAL_RATIO_LOG2: u32 = 3;

/// Number of layers (0..NUM_LAYERS)
pub const NUM_LAYERS: u8 = 16;

/// Chunk size (coordinate span): 2^21 = 2,097,152 bp (~2MB)
/// Sized for good parallelism across 128 cores
pub const CHUNK_SIZE_LOG2: u32 = 21;

/// Derived: log2 of max interval size for layer 0
/// max_size[0] = tile_size[0] / ratio = 2^14 / 2^3 = 2^11 = 2048
pub const LAYER_0_MAX_SIZE_LOG2: u32 = BASE_TILE_SIZE_LOG2 - TILE_MAX_INTERVAL_RATIO_LOG2;

/// Derived: shift value for layer ID computation from interval length
/// For length `len`, layer = max(0, floor(log2(len)) - LAYER_SHIFT)
pub const LAYER_SHIFT: u32 = LAYER_0_MAX_SIZE_LOG2 - 1;

/// Basic interval representation with half-open coordinates [start, end)
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Archive, Serialize, Deserialize,
)]
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
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Archive, Serialize, Deserialize,
)]
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

// Implement Radixable to enable O(n) radix sorting by start coordinate
impl Radixable<u32> for TaggedInterval {
    type Key = u32;

    #[inline]
    fn key(&self) -> Self::Key {
        self.iv.start
    }
}

/// Exponential layer identifier (log2 of interval size range)
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Archive, Serialize, Deserialize,
)]
pub struct LayerID(pub u8);

impl LayerID {
    /// Compute layer ID from interval length
    ///
    /// Maps interval length to the appropriate layer based on centralized constants.
    /// Layer k handles intervals with length in [max_size[k-1], max_size[k])
    /// where max_size[k] = 2^(LAYER_0_MAX_SIZE_LOG2 + k).
    ///
    /// Formula: layer = clamp(floor(log2(len)) - LAYER_SHIFT, 0, NUM_LAYERS-1)
    pub fn from_length(len: u32) -> Self {
        if len == 0 {
            LayerID(0)
        } else {
            // floor(log2(len)) = 31 - leading_zeros for non-zero values
            let log2_len = 31 - len.leading_zeros();
            let layer = log2_len.saturating_sub(LAYER_SHIFT);
            LayerID(layer.min(NUM_LAYERS as u32 - 1) as u8)
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
/// INVARIANT: `intervals` are sorted by start coordinate
/// to enable binary search during query. This sorting is performed by `index_sweep`.
///
/// NOTE: tile_size is set larger than the max interval size for each layer. This guarantees
/// no interval can span an entire tile, eliminating the need for running_counts.
#[derive(Debug, Clone, Archive, Serialize, Deserialize)]
pub struct Tile {
    /// Start coordinate of this tile
    pub start_coord: u32,
    /// Intervals that overlap with this tile
    /// SORTED by start coordinate for binary search via radix sort
    pub intervals: Vec<TaggedInterval>,
}

impl Tile {
    pub fn new(start_coord: u32) -> Self {
        Self {
            start_coord,
            intervals: Vec::new(),
        }
    }

    /// Check if the tile is empty (no intervals touch it)
    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
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
        // With LAYER_SHIFT = 10, layer = max(0, floor(log2(len)) - 10)
        // Layer 0 handles intervals with length < 2048 (2^11)
        assert_eq!(LayerID::from_length(0), LayerID(0));
        assert_eq!(LayerID::from_length(1), LayerID(0));
        assert_eq!(LayerID::from_length(1024), LayerID(0)); // log2(1024)=10, 10-10=0
        assert_eq!(LayerID::from_length(2047), LayerID(0)); // log2(2047)=10, 10-10=0

        // Layer 1 handles [2048, 4096)
        assert_eq!(LayerID::from_length(2048), LayerID(1)); // log2(2048)=11, 11-10=1
        assert_eq!(LayerID::from_length(4095), LayerID(1)); // log2(4095)=11, 11-10=1

        // Layer 2 handles [4096, 8192)
        assert_eq!(LayerID::from_length(4096), LayerID(2)); // log2(4096)=12, 12-10=2
        assert_eq!(LayerID::from_length(8191), LayerID(2));

        // Higher layers
        assert_eq!(LayerID::from_length(1 << 20), LayerID(10)); // 1M -> layer 10
        assert_eq!(LayerID::from_length(1 << 25), LayerID(15)); // 32M -> layer 15 (max)

        // Values beyond max layer are clamped to NUM_LAYERS-1
        assert_eq!(LayerID::from_length(u32::MAX), LayerID(NUM_LAYERS - 1));
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
        non_empty.intervals.push(TaggedInterval::new(50, 60, 1));
        assert!(!non_empty.is_empty());
    }

    #[test]
    fn test_tile_serialization_roundtrip() {
        use rkyv::rancor::Error;

        let mut tile = Tile::new(1000);
        tile.intervals.push(TaggedInterval::new(1050, 1100, 3));
        tile.intervals.push(TaggedInterval::new(1075, 1150, 4));

        // Serialize
        let bytes = rkyv::to_bytes::<Error>(&tile).unwrap();

        // Deserialize
        let archived = rkyv::access::<ArchivedTile, Error>(&bytes).unwrap();

        assert_eq!(archived.start_coord, 1000);
        assert_eq!(archived.intervals.len(), 2);
    }
}
