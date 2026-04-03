# Riggle - Claude Code Context

## Project Overview
Riggle is a statistical interval intersection engine for genomic data. It indexes BED files into a hierarchical structure (layers → chunks → tiles) and supports fast parallel queries with Fisher's exact test statistics.

## Architecture

### Hierarchical Structure
- **Layers**: Partition intervals by size (log2 of length). Each interval goes to exactly ONE layer.
- **Chunks**: Fixed-size genomic territories within a layer. O(1) lookup via `chunk_id = coord / chunk_size`.
- **Tiles**: Subdivisions within chunks. Intervals are classified as:
  - `running_counts`: spans entire tile (sid, count)
  - `start_ivs`: starts in tile - Vec<OffsetSid>, SORTED by offset via radix sort
  - `end_ivs`: ends in tile - Vec<OffsetSid>, SORTED by offset via radix sort

### Key Files
```
src/
├── core/mod.rs      # Interval, TaggedInterval, Tile, LayerID, ChunkID, TileID
├── matrix/
│   ├── dense.rs     # DenseMatrix, BitwiseMask, allocate_dense_accumulator
│   └── sparse.rs    # condense_to_sparse, merge_sparse
├── io/
│   ├── header.rs    # MasterHeader, ChunkHeader, LayerConfig
│   ├── mmap.rs      # MappedChunk (zero-copy via memmap2)
│   └── parse.rs     # BED parsing
├── sweep/
│   ├── index.rs     # index_sweep - builds tiles from intervals
│   └── query.rs     # query_sweep - counts intersections
├── tasks/
│   ├── build.rs     # build_database pipeline
│   └── query.rs     # query_database pipeline
├── stats/mod.rs     # Fisher's exact test
└── main.rs          # CLI (build, query, info commands)
```

## Critical Implementation Details

### Offset Storage (FIXED)
- `start_ivs` and `end_ivs` use `OffsetSid` struct (offset: u32, sid: u32)
- `OffsetSid` implements `Radixable<u32>` for O(n) sorting by offset
- Previously used u16 for offset which caused silent truncation for large tiles (>65535)
- Tile sizes can exceed 65535 for layers > 11

### Binary Search in Tiles
- `start_ivs` and `end_ivs` are sorted by offset at index time (in `index_sweep`)
- Query uses `partition_point()` for O(log n) lookup to initiate iteration

### O(1) Chunk Lookup
- Chunk paths are predictable: `layer_{layer_id}/chunk_{chunk_id}.bin`
- `compute_chunk_tasks()` in query.rs computes which chunks to access from coordinates
- No filesystem scanning needed

### Deduplication (FIXED)
- Uses tile position logic: `is_first_overlap_tile = query.iv.start >= tile_start`
- Only counts running_counts and end_ivs in the first overlap tile
- start_ivs always counted (unique to their tile)
- O(1) determination instead of O(N) seen.contains()

### Thread-Local Buffers (FIXED)
- Uses `thread_local!` in tasks/query.rs to reuse DenseMatrix/BitwiseMask
- Clears only flagged regions between uses via `zero_flagged_regions()`
- Buffers resize only when dimensions increase

### Sparse Matrix Construction (FIXED)
- Builds CsMat directly from row-major ordered data in O(nnz)
- No TriMat intermediate, no O(N log N) sort needed
- Data is already sorted: rows in order, cols in order within rows

## Build & Test
```bash
nix develop --command cargo test
nix develop --command cargo build --release
```

Note: Requires nix shell for libiconv on macOS.

## CLI Usage
```bash
riggle build --input <dir> --output <db_path>
riggle add --input <dir> --db <db_path>           # Add files to existing database
riggle query --db <db_path> --query <bed_or_dir> --output <json>
riggle info --db <db_path>
```

### Query Batch Mode
Query accepts either a single BED file or a directory of BED files:
- Single file: `--query file.bed`
- Directory: `--query query_dir/` (all .bed files processed)

## Completed Optimizations

### From docs/_critiques.md
1. ✅ u16 truncation bug - Changed to u32 for tile offsets
2. ✅ O(N) deduplication - Uses tile position logic (`is_first_overlap_tile`)
3. ✅ Thread-local buffer reuse - Uses `thread_local!` for DenseMatrix/BitwiseMask
4. ✅ O(N log N) sparse construction - Builds CsMat directly in O(nnz)

### From Code TODOs
5. ✅ Radix sort for tile intervals - Uses voracious_radix_sort with OffsetSid wrapper
6. ✅ Radix sort for build pipeline - TaggedInterval implements Radixable for O(n) layer sorting
7. ✅ Pre-allocation in build pipeline - Exact capacity via two-pass
8. ✅ Add command - Incremental database updates
9. ✅ Batch query support - Query with directory of BED files

## rkyv Notes
- Uses rkyv 0.8 for zero-copy deserialization
- Archived tuples accessed via `entry.0.into()`, `entry.1.into()` (not destructuring)
- ArchivedTile, ArchivedInterval etc. for zero-copy access from mmap
