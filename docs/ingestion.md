## Riggle Technical Specification: Ingestion & Indexing Algorithms

## 1. Unified Interval Ingestion Pipeline

Both the indexing and query pipelines must resolve raw input files into a structured, layer-partitioned memory layout before further processing. To guarantee consistency and maximize code reuse, Riggle utilizes a unified ingestion routine for all incoming interval data.

Before processing begins, the entrypoint parses the global `meta.json` to load the layer configuration (specifically `min_size` and `growth_factor`), which dictates how intervals will be spatially partitioned.

### The Ingestion Flow

1.  **Parallel Flat-Mapping:** A batch of input files (TSVs) is read and parsed in parallel. Every valid record is transformed into the universal Riggle primitive: `(start, end, sid)`. The parser pairs each interval with its corresponding `shardName` (e.g., chromosome or coordinate space identifier).
2.  **Layer Assignment:** Using the parsed layer configuration, each interval's size (`end - start`) is evaluated to determine its target `layerID`.
3.  **Partitioning & Grouping:** The engine groups the flat-mapped intervals by their target destination.
4.  **The Output:** The pipeline yields a mapping for the entire batch:
    `{ (shard, layerID) -> vec<Interval> }`

---

## 2. Parallel Sort-Merge Strategy

Riggle requires strictly sorted intervals to function efficiently. Before the mapped data can be merged into the index or traversed during a query, the vectors within the `{ (shard, layerID) -> vec<Interval> }` map must be sorted by their `start` coordinate.

To achieve maximum throughput, Riggle employs a nested parallel sort-merge strategy:

- **Outer Parallelism:** The engine parallelizes across the map entries (the distinct `(shard, layerID)` pairs) using nested `iter_par()` calls.
- **Inner Parallelism (Voracious Radix):** Within each vector, chunks are sorted utilizing **`voracious_radix_sort`**. This provides a significant performance boost over standard comparison sorts for flat coordinate data, rapidly ordering the `(start, end, sid)` tuples.

---

## 3. Indexing Architecture & Thread Delegation

At index time, the system architecture separates the CPU-bound parsing/sorting from the I/O-bound disk operations to prevent bottlenecking.

1.  **The Main Thread:** Responsible for batching the input file paths and executing the unified ingestion and sorting pipelines.
2.  **Handoff:** Once a batch is fully processed and sorted, the main thread transmits the resulting slice of data—structured as `&[(shard, layer, sorted_ivs)]`—to a dedicated **Index Writer Thread**.
3.  **The Index Writer Thread:** Exclusively handles the memory-mapped file mutations, ensuring thread-safe, sequential writes to disk without blocking the main thread's ingestion loop.

---

## 4. The Zero-Spike Merge Algorithm (Index Insertion)

Because both the existing memmap and the incoming batch vectors are strictly sorted, inserting new data is logically a merge operation. However, loading an entire multi-gigabyte layer memmap into RAM to perform a standard merge would cause catastrophic memory spikes.

To solve this, the Index Writer thread uses an in-place **Extend & Reverse-Shift** strategy to merge data directly on disk.

### Step-by-Step Merge Process

1.  **Memmap Extension:** The existing memory-mapped file for the target `(shard, layerID)` is extended on disk to accommodate the exact byte size of the incoming `sorted_ivs`. The new space at the end of the file is initially zero-filled.
2.  **Chunked Reading:** Instead of loading the whole file, Riggle reads the existing data in discrete, memory-safe chunks.
3.  **Reverse Merge & Overwrite (Down-Shift):** To prevent overwriting unread existing data, the merge is executed from **bottom-to-top** (end-to-start).
    - The engine simultaneously iterates the layer memmap and `sorted_ivs` in reverse order.
    - Using the merge routine (as used in merge sort), a sorted intermediate buffer is accumulated.
    - Once this buffer is fully populated (or a data source is exhausted), the buffer can be written to the end of the memmap.
    - This reverse-merge safely populates the file from back to front. As the pointers move backwards, the algorithm naturally creates the necessary gaps for the new intervals without ever clobbering data that hasn't been read yet.
