### 1. High-Level Module Architecture

- **`riggle::core`**: Fundamental domain types (`Interval`, `TaggedInterval`, `TileID`, `LayerID`).
- **`riggle::io`**: Memory-mapped file wrappers, zero-copy deserialization (`rkyv` or `zerocopy`), and chunk boundary routing.
- **`riggle::sweep`**: The isolated, bufferless line-sweep algorithms (Indexing and Querying variants).
- **`riggle::matrix`**: Dense matrix allocation, bitmask tracking, and sparse condensation logic.
- **`riggle::tasks`**: Thread pool execution units and concurrency primitives.
- **`riggle::stats`**: Downstream statistical aggregation (Fisher's Exact, Jaccard, etc.).

---

### 2. Thread Pool Task Typology

Riggle operates as a Directed Acyclic Graph (DAG) of tasks submitted to a central work-stealing thread pool (e.g., Rayon).

#### A. Indexing Tasks (Build Phase)

- **`ParseAndFlatten`**: Reads a raw database file, validates coordinates, and outputs a `Vec<TaggedInterval>`.
- **`LayerPartition`**: Takes flattened batches and routes them into layer-specific buckets based on $\log_2(\text{size})$.
- **`RadixSort`**: Sorts a partition locally by `start` coordinate. **Defer to fast (non-comparative sort) library.**
- **`ChunkWrite`**: Takes a sorted slice of intervals and a target Chunk ID. Executes the `IndexSweep` and serializes the resulting `Tile` structures to disk.

#### B. Query Tasks (Query Phase)

- **`QueryFlattenAndSort`**: Parses a query set, creates `TaggedInterval`s, and Radix sorts them.
- **`ChunkQuery`**: The primary workhorse. Takes a sorted query batch and a memory-mapped chunk. Allocates a thread-local dense matrix, runs the `QuerySweep`, and condenses it into a `SparseMatrix`.
- **`MatrixMerge`**: A reduction task. Merges multiple `SparseMatrix` objects. Executed in a tree-reduction pattern across chunks, and then across layers.
- **`StatCompute`**: Takes a finalized `SparseMatrix` row (one Query SID against all Database SIDs) and computes the exact statistical significance of the overlaps.

---

### 3. Key Function Signatures

#### Module: `riggle::core`

```rust
pub struct Interval {
    pub start: u32,
    pub end: u32,
}

pub struct TaggedInterval {
    pub iv: Interval,
    pub sid: u32,
}
```

#### Module: `riggle::sweep`

This module encapsulates the zero-buffer, slice-based traversal logic.

```rust
/// INDEXING SWEEP
/// Translates a sorted batch of intervals into serialized tile states.
///
/// * `chunk_bounds`: The absolute genomic coordinates of this chunk.
/// * `tile_size`: The fixed size of tiles in this layer.
/// * `active_batch`: The slice of ALL intervals routed to this chunk, sorted by start.
///
/// Returns: A vector of serialized Tile byte arrays ready for disk I/O.
pub fn index_sweep(
    chunk_bounds: Interval,
    tile_size: u32,
    active_batch: &[TaggedInterval]
) -> Vec<SerializedTile> { ... }


/// QUERY SWEEP
/// Executes the read-only sweep against a memory-mapped chunk.
///
/// * `mapped_chunk`: Zero-copy reference to the mmap'd chunk data.
/// * `query_batch`: Sorted slice of query intervals.
/// * `results`: Mutable reference to the thread-local dense matrix.
/// * `mask`: Mutable reference to the bitmask tracking non-zero coordinates.
pub fn query_sweep(
    mapped_chunk: &MappedChunk,
    query_batch: &[TaggedInterval],
    results: &mut DenseMatrix,
    mask: &mut BitwiseMask
) { ... }
```

#### Module: `riggle::matrix`

Handles the memory-efficient aggregation discussed in the technical spec.

```rust
/// Allocates a dense matrix optimized for continuous cache-friendly writes.
/// Dimensions: (batch_query_sids) x (total_database_sids)
pub fn allocate_dense_accumulator(
    num_queries: usize,
    num_db_sids: usize
) -> (DenseMatrix, BitwiseMask) { ... }

/// Condenses the highly-mutated dense matrix into a sparse representation
/// immediately after a ChunkQuery finishes, using the bitmask for $O(K)$ traversal.
pub fn condense_to_sparse(
    dense: &DenseMatrix,
    mask: &BitwiseMask
) -> SparseMatrix { ... }

/// Merges two sparse matrices. Used by the `MatrixMerge` reduction tasks.
pub fn merge_sparse(
    left: SparseMatrix,
    right: SparseMatrix
) -> SparseMatrix { ... }
```

#### Module: `riggle::tasks` (The Executor Interface)

Defines how the pipeline is orchestrated.

```rust
pub enum RiggleTask {
    // Build Tasks
    ParseAndFlatten { filepath: PathBuf, db_sid: u32 },
    LayerPartition { batch: Vec<TaggedInterval> },
    ChunkWrite { chunk_id: ChunkID, data: Vec<TaggedInterval> },

    // Query Tasks
    QueryFlattenAndSort { queries: Vec<Interval>, query_sid: u32 },
    ChunkQuery {
        chunk_id: ChunkID,
        query_batch: Arc<Vec<TaggedInterval>> // Arc allows multiple chunk tasks to share the query batch
    },
    MatrixMerge { a: SparseMatrix, b: SparseMatrix },
    StatCompute { final_matrix: SparseMatrix }
}

/// The core scheduler loop
pub fn execute_pipeline(tasks: Vec<RiggleTask>, thread_pool: &mut ThreadPool) { ... }
```

### 4. Implementation Details for the Query Thread Execution

When a thread pulls a `ChunkQuery` task off the queue, its internal lifecycle looks exactly like this:

1.  **Acquire Memory:** Thread checks its local cache for a reusable `DenseMatrix` and `BitwiseMask` buffer. If none exists, it calls `allocate_dense_accumulator`.
2.  **Zero-Out:** It zeroes out only the regions of the dense matrix that were flagged by the bitmask in the _previous_ run (preventing a full $O(N \times M)$ wipe).
3.  **Sweep:** It calls `query_sweep()`, passing the `mapped_chunk` and the `query_batch`. The algorithm moves the `head` pointer along the slice, updating the dense matrix and flipping bits in the mask in-place.
4.  **Condense:** It calls `condense_to_sparse()`.
5.  **Return:** The `SparseMatrix` is returned to the thread pool channel, triggering a `MatrixMerge` task, while the thread retains its Dense buffers for the next `ChunkQuery` task.
