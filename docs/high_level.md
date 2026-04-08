## Riggle Technical Specification

## 1. Core Paradigm

Riggle is a statistical interval intersection engine focused on **set-level analysis**. It prioritizes memory-bandwidth and zero-copy access. Riggle does not return individual overlapping intervals; instead, it rapidly computes intersection counts across sets to power downstream statistical analysis.

In Riggle, everything is a **SetID (SID)**:

- **Database SID:** Represents a set of intervals (e.g., a specific ChIP-seq experiment).
- **Query SID:** Represents a set of query intervals.
- **The Goal:** Compute an intersection matrix where each cell is the count of overlaps between a specific Query SID and a Database SID.

---

## 2. Data Hierarchy & Memory Layout

The database is sharded by chromosome/contig. Each shard is distributed across **layers** (exponentially bounded interval sizes) and **chunks** (memory-mapped files representing fixed genomic territories).

### A. Master Header (Global Metadata)

- **SID Map:** `database_sid -> metadata`.
- **Layer Configurations:** Defines exponential tiling bounds and chunk sizes per layer.
- **Shards:** List of shard names and their maximum coordinates.

### B. Chunk Layout (The Memmap)

A chunk is a memory-mapped file consisting of **rkyv blocks** (referred to as **tiles**).

- **Chunk Header:** `tile_id -> byte_offset`. Provides $O(1)$ routing.
- **Tile Data:** Tightly packed, serialized byte arrays.

### C. Tile Composition (The Leaf Node)

Tiles track the state of intervals relative to tile boundaries.

- **`start_coord`:** The absolute genomic coordinate where this tile begins.
- **`intervals`:** `[(start, end, sid)]`. Sorted by start coordinate.

---

## 3. Indexing & Chunk Updates

### Phase 1: Map & Sort (per shard)

Intervals are partitioned into exponential layers based on size and sorted by `start` coordinate. Intervals are then grouped by target `Chunk`.

### Phase 2: The Update Sweep (`ChunkWrite`)

To update an existing chunk without the memory spike of loading the entire file, Riggle uses a **Shift-and-Sync** strategy:

1.  **In-Memory Delta:** Form an in-memory chunk containing only the _new_ intervals.
2.  **Offset Calculation:** Calculate the new byte offsets for every tile in the chunk (existing + new).
3.  **Memmap Extension:** Extend the existing memory-mapped file on disk to accommodate the new size.
4.  **Back-to-Front Shift:** Move existing tiles to their new offsets starting from the end of the file. This prevents overwriting data before it is moved.
5.  **Synchronous Write:** Write the new tiles into the vacated gaps and update the header.

---

## 4. Query Pipeline

### Phase 1: High-Level Set Preparation

Queries are handled strictly at the **set level**.

1.  **Flattening:** A batch of query sets is flattened into a list of **intervals**.
2.  **Tagging:** Each interval is paired with its original **Query SID**. This allows the engine to operate on raw intervals for speed while preserving the ability to remap results back to the original sets for statistical accumulation.

### Phase 2: The General Form

The query algorithm follows a **parallel tree map-reduce** strategy over chunks.

1.  **Parallel Map:** Queries are dispatched across layers and chunks.
2.  **Deferred Statistics:** The map phase omits the generation of dense matrices for interval densities. Instead, it yields raw intersection counts per SID pair.
3.  **Tree Reduce:** Results are aggregated across tiles, chunks, layers, and shards. Statistical tests (e.g., Fisher’s Exact) are deferred until all layers have yielded their final set-level counts.

---

## The Riggle Scan Algorithm (Index Time)

The core of Riggle's indexing is a one-dimensional sweep-line algorithm. Because both the tiles and input intervals are strictly sorted, we resolve spatial relationships in a single forward pass.

### The Universal State

The indexing thread tracks active intervals via slice manipulation:

1.  **The Head Pointer:** Points to the first interval intersecting current tile $T$.
2.  **In-Place Pruning:** As the sweep moves forward, intervals that end before the current tile starts are logically dropped.
3.  **Tile Intersection:** For each tile $T$, Riggle identifies intervals $I$ where $I.start < T.end$ and $I.end > T.start$.
4.  **Serialization:** These intervals are pushed to $T.intervals$ and serialized using `rkyv`.

_**Note**: Because tile_size $\ge$ max_interval_size for a given layer, no interval can span an entire tile without having at least one endpoint within it._
