## Riggle Technical Specification: Query Architecture

### System Architecture Overview

Riggle's query engine computes dense intersection statistics across entire sets without allocating intermediate overlapping intervals. It leverages massive bounded exponential layer memmaps for optimal disk I/O, a single-pass dual-active-set sweep for cache efficiency, dense matrix SIMD accumulation for zero-allocation scaling, and a binary-search fast-forward heuristic to aggressively bypass dead space.

---

### 1. Data Structures & Memory Layout

**1.1. Storage Engine (Disk / OS Page Cache)**

- **Database ($D$):** Partitioned into distinct coordinate spaces (shards) and divided into **Exponential Layers**. There are no chunks or tiles. Each layer is a single, giant memory-mapped file.
  - _Constraint:_ Every interval in Layer $K$ has length $L \le \text{layerSize}$ (the maximum length boundary for that layer).
  - _Layout:_ Flat, contiguous array of universally formatted structs: `(start: u32, end: u32, D_sid: u32)`. Strictly sorted by `start`.
- **Queries ($Q$):** Fully loaded into RAM as a flat array.
  - _Layout:_ `(start: u32, end: u32, Q_sid: u32)`. Strictly sorted by `start`.

**1.2. Execution Engine (Thread-Local RAM)**

- **`Active_D` (The Shield):** An Implicit Min-Heap ordered by `end` coordinate. Guards against density spikes and tracks when database intervals expire.
- **`Active_Q` (The Ring):** An Implicit Min-Heap or ring buffer ordered by `end` coordinate. Tracks active query intervals (low depth expected).
- **`overlapsFrame` (The Delta State):** A dense 1D array `[D_sids]u32` representing the current frequency of all active database tracks at the exact sweep line coordinate.
- **`results` (The Accumulator):** A dense 2D matrix `[Q_sids][D_sids]u32` pre-allocated per thread.

---

### 2. The Core Execution Loop (Single-Threaded Block Processing)

This loop is executed independently by each thread over a designated block of query intervals against a specific exponential layer memmap.

#### 2.1. The Fast-Forward Escapement (Sparsity Heuristic)

Because threads map over query blocks rather than DB chunks, a thread's assigned queries may start deep within the chromosome. To avoid scanning millions of irrelevant database intervals, the engine utilizes a distance-based fast-forward.

```text
LOOP START (Check Dead Space):
    IF D_cursor.start < (next_Q.start - layerSize):
        // We are in a dead zone. The current DB interval cannot possibly overlap the next Query.
        1. CLEAR Active_D heap.
        2. ZERO out overlapsFrame.
        3. BINARY SEARCH the layer memmap for the first D index where:
           D.start >= (next_Q.start - layerSize)
        4. D_cursor = found_index
```

#### 2.2. The Sweep Line State Machine

Once aligned, the sweep line advances to time $T$, defined as the minimum coordinate among the next four potential events: `D_cursor.start`, `Q_cursor.start`, `peek(Active_D).end`, `peek(Active_Q).end`.

Events are processed strictly in this order to handle boundary conditions (Ends process before Starts at the exact same coordinate):

**Event 1: $D$ Ends ($T = \text{peek}(\text{Active\_D}).\text{end}$)**

1. `expired_D = Active_D.pop()`
2. `overlapsFrame[expired_D.sid] -= 1`

**Event 2: $Q$ Ends ($T = \text{peek}(\text{Active\_Q}).\text{end}$)**

1. `Active_Q.remove()`

**Event 3: $D$ Starts ($T = \text{D\_cursor.start}$)**

1. `Active_D.push(D_cursor)`
2. `overlapsFrame[D_cursor.sid] += 1`
3. _Forward Match (Sparse Update):_
   `FOR q IN Active_Q:`
   `results[q.sid][D_cursor.sid] += 1`
4. Advance `D_cursor`

**Event 4: $Q$ Starts ($T = \text{Q\_cursor.start}$)**

1. `Active_Q.push(Q_cursor)`
2. _Catch-Up Match (SIMD Vector Addition):_
   `results[Q_cursor.sid] += overlapsFrame` // CPU executes AVX-512 block add
3. Advance `Q_cursor`

---

### 3. Parallel Aggregation & Output

To scale across massive core counts without memory bandwidth starvation or lock contention, the engine uses a query-partitioned map-reduce pipeline.

**Phase 1: Map (Query Block Processing)**

- The global flat array of query intervals ($Q$) is partitioned into evenly sized continuous blocks.
- A work-stealing scheduler (e.g., Rayon) assigns these query blocks to idle threads.
- Each thread runs the Core Execution Loop against the layer memmap. Due to the fast-forward heuristic, the thread instantly skips irrelevant memmap space and seeks to the region relevant to its assigned query block.
- Each thread generates a local, fully populated `results` dense matrix. No locks, no atomics, zero false sharing.

**Phase 2: Reduce (Parallel Tree Merge)**

- As threads finish their blocks, the scheduler dynamically pairs them.
- Idle threads take two thread-local `results` matrices and perform a SIMD vectorized matrix addition: `Matrix_A += Matrix_B`.
- This folds $N$ thread-local matrices into 1 final global results matrix in $O(\log(\text{threads}))$ steps.

**Phase 3: Sparse Export**

- The final global dense matrix is converted to a sparse format (e.g., COUP/Coordinate List or directly to a TSV/Parquet file) for the end-user, filtering out the millions of zeros. Only the requested statistical overlaps remain.
