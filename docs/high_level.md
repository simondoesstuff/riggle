# Riggle Technical Specification

## 1. Core Paradigm

Riggle is a statistical interval intersection engine, an evolution of IGD that incorporates elements from Giggle. It prioritizes memory-bandwidth and zero-copy access. Riggle does not return individual overlapping intervals; instead, it rapidly computes intersection matrices (`Query Set ID` $\times$ `Database Set ID` $\rightarrow$ `Overlap Count`) to power downstream statistical analysis.

## 2. Data Hierarchy & Memory Layout

The database is sharded by chromosome/contig (detected from the first non-integer column in BED files).
Each shard is further distributed across layers (exponentially bounded interval sizes) and chunks (memory-mapped files representing fixed genomic territories).

### A. Master Header (Global Metadata)

- **SID Map:** `database_sid -> metadata` (name, filepath, etc.)
- **Layer Configurations:** Defines the exponential tiling bounds and chunk sizes per layer.
- **Shards:** List of shard names (e.g., `["chr1", "chr2", "chrX"]`)
- **Shard Max Coords:** Maximum coordinate seen per shard.

### A.1 Shard Detection (BED Parsing)

Columns are scanned left-to-right:
- First non-integer column → shard name (e.g., "chr1", "chrX")
- First two integer columns → (start, end) coordinates
- If no non-integer column exists → uses `"default"` as shard name

This flexible parsing supports standard BED (chrom, start, end) as well as coordinate-only formats.

### B. Chunk Layout (The Mmap)

- **Chunk Header:** `tile_id -> byte_offset`. Provides $O(1)$ routing to a specific tile.
- **Tile Data:** The tightly packed, serialized byte arrays of the individual tiles.

### C. Tile Composition (The Leaf Node)

Tiles track the _state_ of intervals relative to the tile boundaries. Data is serialized using a zero-copy framework (e.g., `rkyv`).

- **`start_coord`:** The absolute genomic coordinate where this tile begins.
- **`intervals`:** `[(start, end, sid)]`. Intervals that intersect with this tile, sorted by their start coordinate.

**Key Invariant:** Tile size is set to 4× the maximum interval size for each layer. This guarantees no interval can completely span a tile, eliminating the need for `running_counts` and simplifying both indexing and query logic.

## 3. Indexing Pipeline (`ChunkWrite`)

**Phase 0: Parse & Shard**

Intervals are parsed from BED files and grouped by shard (chromosome). Each shard is processed independently into `{db}/{shard}/layer_{id}/chunk_{id}.bin`.

**Phase 1: Map & Sort (per shard)**

1. **Flatten & Tag:** Intervals are parsed and flattened into `(interval, sid)`.
2. **Layer Splitting:** Intervals are partitioned into their appropriate exponential layers based on their size.
3. **Sort:** The flattened array undergoes a Radix sort by the `start` coordinate.
4. **Group:** Intervals are grouped by target `Chunk`. Boundary-crossing intervals are duplicated into all necessary chunks.

**Phase 2: Reduce & Write (The Sweep)**

- A worker thread takes a `Chunk` and its sorted `[(interval, sid)]` list.
- It performs the **Riggle Scan** (detailed below) to construct the `intervals` array for each tile.
- The tile is serialized to disk at `{shard}/layer_{id}/chunk_{id}.bin`, updating the Chunk Header.

## 4. Query Pipeline (`ChunkQuery`)

**Phase 0: Shard Grouping**

Query intervals are parsed and grouped by shard. Only shards present in both the query AND database are processed (missing shards are silently skipped—no overlaps possible).

**Phase 1: Query Preparation (per shard)**

1. **Flatten & Tag:** Query intervals are flattened into `(global_index, interval)` preserving indices for aggregation.
2. **Layer-Wise Splitting:** _Crucial Step._ A query interval must be tested against all layers because overlaps can happen irrespective of size.
3. **Results Initialization:** The thread initializes a dense results matrix of dimensions `(batch_query_sids, database_sids)`. Because we expect fewer than 100k database SIDs and significantly fewer query SIDs per batch, a dense matrix is highly memory-efficient for continuous writes. Concurrently, a bit mask is initialized to track non-zero regions during the sweep.

**Phase 2: The Sweep & Aggregate (per shard)**

- The thread executes the **Riggle Scan** (detailed below) across the mapped tiles.
- Results are mutated in-place within the thread's dense matrix, while the bit mask flags intersecting regions.
- Following the sweep, the dense matrix is efficiently condensed into a sparse matrix by utilizing the bit mask to locate non-zero data.

**Phase 3: Cross-Shard Aggregation**

- Thread-local sparse matrices are merged (across chunks, across layers, across shards) via tree-reduce.
- Statistics (e.g., Fisher's Exact) are computed on the aggregated totals.
- JSON results are returned.

---

# The Riggle Scan Algorithm

The core of Riggle is a one-dimensional sweep-line algorithm. Because both the tiles in a chunk and the input intervals (query or database batches) are strictly sorted by their start coordinates, we can resolve all spatial relationships in a single forward pass without secondary buffer allocations.

### The Universal State

Regardless of whether Riggle is indexing or querying, the thread operates directly on the sorted batch of input intervals for a specific chunk. For any given tile $T$, the thread tracks the active intervals entirely via slice/pointer manipulation:

1. **The Head Pointer:** A pointer (or index) into the batch representing the _first_ interval that intersects the current tile $T$.
2. **Evaluating Overlaps:**

   - Everything before the head has started earlier and may be overlapping the tile.

   * Everything `tile_size` past the head starts after the head and is considered strictly inactive for this tile.

3. **In-Place Pruning:** As the sweep traverses forward across tiles, if it hits an interval's end, that interval is logically removed from the active subset. This guarantees that at any given time, there are no intervals that both start and end before the head.
4. **Trivial Calculation:** Because expired intervals are continuously dropped, the algorithm can trivially determine if the remaining active intervals within the pointer bounds started before the tile or within the tile.

### Variant A: Indexing (`ChunkWrite` Scan)

_Goal: Translate a stream of intervals into tile-specific state arrays._

For the current Tile $T$, and for each active interval $I$ identified relative to the head pointer:

- **Check Intersection:** If $I.start < T.end$ AND $I.end > T.start$:
  - Push `(I.start, I.end, I.sid)` to $T$`.intervals`.

_Note: Since tile_size ≥ 4× max_interval_size, no interval can span an entire tile. Each interval will have at least one endpoint (start or end) within a tile, never both outside._

### Variant B: Querying (`ChunkQuery` Scan)

_Goal: Calculate exact intersection counts while strictly preventing double-counting of duplicated database intervals._

We can guarantee no duplication without hashes or bit masks with a simple logical rule. The overlap between a query $Q$ and database interval $D$ begins exactly at `max(Q.start, D.start)`. Since tiles perfectly partition the coordinate space, this overlap start coordinate falls into exactly one tile. We simply only count the overlap if it "begins" in the current tile.

For the current Tile $T$, we evaluate the active query intervals. For a given active query $Q$:

For each $D \in T$`.intervals`:

- If $Q$ overlaps $D$ (i.e., $Q.start < D.end$ AND $Q.end > D.start$):
  - _Deduplication Check:_ We only count this if the intersection _began_ in this tile. The intersection begins at $S = \max(Q.start, D.start)$.
  - If $S \ge T.start$ AND $S < T.end$, increment `results[Q.sid][D.sid]`.
