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
- **`running_counts`:** `{database_sid -> count}`. Tracks intervals that completely span the tile.
- **`start_ivs`:** `[(start_offset, database_sid)]`. Intervals that begin within this tile.
- **`end_ivs`:** `[(end_offset, database_sid)]`. Intervals that terminate within this tile.

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
- It performs the **Riggle Scan** (detailed below) to construct the `running_counts`, `start_ivs`, and `end_ivs` for each tile.
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

1. **Check Start:** If $I.start \ge T.start$ AND $I.start < T.end$:

   - Push `(I.start - T.start, I.sid)` to $T$`.start_ivs`.

2. **Check End:** If $I.end > T.start$ AND $I.end \le T.end$:

   - Push `(I.end - T.start, I.sid)` to $T$`.end_ivs`.

3. **Check Span:** If $I.start < T.start$ AND $I.end > T.end$:

   - Increment the count for $I.sid$ in $T$`.running_counts`.

_Note: Because $I$ is duplicated across chunk boundaries prior to this step, an interval that crosses from Chunk A to Chunk B will correctly register an `end_ivs` event in Chunk A and a `start_ivs` event in Chunk B if it terminates/begins within those boundary tiles._

### Variant B: Querying (`ChunkQuery` Scan)

_Goal: Calculate exact intersection counts while strictly preventing double-counting of duplicated database intervals._

We can guarantee no duplication without hashes or a bit masks with a simple logical rule, borrowed from IGD. Only take an interval as overlapping if either the query or database interval start in the current chunk, that is, at most one of them can intersect the left boundary. We encode this as follows...

For the current Tile $T$, we evaluate the active query intervals. For a given active query $Q$:

**1. Apply Running Counts (The Fast Path)**

If $Q$ completely spans $T$ ($Q.start \le T.start$ AND $Q.end \ge T.end$):

- $Q$ intersects _everything_ in this tile.
- Add $T$`.running_counts` directly to `results[Q.sid]`.
- Add all `database_sid`s present in $T$`.start_ivs` and $T$`.end_ivs` to `results[Q.sid]`.
- _IGD Deduplication Check:_ Only apply the above if $T$ is the _first_ tile where this intersection occurs (i.e., $T$ is the tile containing $\max(Q.start, \text{Tile.start})$).

**2. Process Sub-Tile Starts ($T$.`start_ivs`)**

For each $D_{start} \in T$`.start_ivs`:

- If $Q$ overlaps $D_{start}$ (i.e., $Q.start < D.start$ AND $Q.end > D.start$):
  - Increment `results[Q.sid][D.sid]`.
  - _No deduplication check needed:_ Because we are intersecting exactly on the database interval's start coordinate, it is impossible for this specific intersection to have been counted in a previous tile.

**3. Process Sub-Tile Ends ($T$.`end_ivs`)**

For each $D_{end} \in T$`.end_ivs`:

- If $Q$ overlaps $D_{end}$ (i.e., $Q.start < D.end$ AND $Q.end > D.end$):
  - _IGD Deduplication Check:_ We only count this if the intersection _began_ in this tile. If $Q$ crosses the left boundary of $T$ ($Q.start < T.start$), this intersection was already counted by a `running_counts` or `start_ivs` evaluation in a previous tile.
  - If $Q.start \ge T.start$, increment `results[Q.sid][D.sid]`.
