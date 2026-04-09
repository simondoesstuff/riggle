## Riggle Technical Specification: System Overview

## 1. Core Paradigm

Riggle is a domain-agnostic, statistical interval intersection engine focused strictly on **set-level analysis**. While capable of handling genomic data, Riggle is generalized for any 1D coordinate system. It prioritizes memory-bandwidth and zero-copy access to rapidly compute intersection counts across sets, powering downstream statistical analysis. Riggle is not concerned with returning individual overlapping intervals or calculating internal interval densities.

### The Universal Primitive

Everywhere in the Riggle codebase, an interval is universally represented as a simple, flat tuple: **`(start, end, sid)`**.

The meaning of the SetID (`sid`) depends entirely on the operational context:

- **Index Time (Database SID / D_SID):** The `sid` represents a distinct indexed item or target feature set in the database.
- **Query Time (Query SID / Q_SID):** The `sid` represents a specific input query set being compared against the database.

**Input & Output:**

- **Inputs:** Both index and query operations ingest a **batch** of interval files. An interval file is a TSV containing `start` and `end` coordinates, and optionally a `shard` column (a discrete coordinate space or partition identifier).
- **The Goal:** The ultimate output of a query is a statistical intersection matrix. The matrix represents Query SIDs by Database SIDs, populating each cell with overlap counts and computed statistics.

---

## 2. Data Hierarchy & Memory Layout

The database is heavily simplified for contiguous memory access. It is partitioned first by **shard** (a distinct coordinate space), and then divided into **layers** based on exponentially bounded interval sizes.

### A. The Layer (The Memmap)

Chunks and tiles do not exist in this architecture. A layer is a **single, giant memory-mapped file**.

- **Structure:** It contains a flat, tightly packed array of intervals.
- **Interval Tuple:** Every entry is strictly a `(start, end, sid)` tuple.
- **Ordering:** Intervals within the memmap are strictly sorted by their `start` coordinate.

### B. Global Metadata (`meta.json`)

Instead of per-layer headers, the entire index is governed by a single `meta.json` file. This file contains only the essential global context:

- **SID Map:** A mapping of `D_SID -> metadata`, storing the system context or identity of each indexed set.
- **Layer Configuration:** The exponential scaling rules for the index, defined simply by a `min_size` and a `growth_factor` (represented as a `u32`).

---

## 3. High-Level Indexing Pipeline

The indexer processes a batch of interval TSV files into the memmap structure through a highly parallelized pipeline:

1.  **Parse & Split:** The batched TSVs are parsed in parallel. Intervals are partitioned into their appropriate shards and layers based on their size and coordinate space.
2.  **Sort:** Intervals within each partition are sorted by their `start` coordinate.
3.  **Memmap Write:** The sorted `(start, end, sid)` tuples are written directly to their corresponding layer's giant memmap. During this index phase, the `sid` attached to each interval corresponds to its original **Database SID**.

---

## 4. High-Level Query Pipeline

Queries operate strictly at the set level to generate the final statistical matrix:

1.  **Batch Preparation (Flatten & Tag):** A batch of query interval files is ingested and flattened. Each query interval is tagged with its originating **Query SID**, forming the standard `(start, end, sid)` tuple. This uniform structure allows Riggle to execute raw interval mathematics at maximum speed while preserving the mapping necessary to accumulate set-level statistics.
2.  **Parallel Map-Reduce:**
    - **Map:** The `(start, end, sid)` query intervals are dispatched in parallel across the relevant shards and layer memmaps to compute raw intersection counts.
    - **Reduce:** Results are aggregated across the engine.
3.  **Matrix Generation:** The pipeline culminates in a matrix mapping Query SIDs to Database SIDs, containing the final aggregated overlap counts and the requested statistical tests.
