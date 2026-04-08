### System Architecture Overview

It leverages bounded exponential layers for optimal disk I/O, a single-pass dual-active-set sweep for cache efficiency, dense matrix SIMD accumulation for zero-allocation scaling, and a fast-forward heuristic to bypass dead space.

---

### 1. Data Structures & Memory Layout

**1.1. Storage Engine (Disk / OS Page Cache)**

- **Database ($D$):** Partitioned into massive Memory-Mapped (memmap) chunks. Within each chunk, data is split into **Exponential Layers**.
  - _Constraint:_ Every interval in Layer $K$ has length $L \le T_k$ (the max length for that layer).
  - _Layout:_ Flat, contiguous array of structs: `(start: u32, end: u32, D_sid: u32)`. Strictly sorted by `start`.
- **Queries ($Q$):** Fully loaded into RAM as a flat array.
  - _Layout:_ `(start: u32, end: u32, Q_sid: u32)`. Strictly sorted by `start`.

**1.2. Execution Engine (Thread-Local RAM)**

- **`Active_D` (The Shield):** An Implicit Min-Heap ordered by `end` coordinate. Guards against density spikes and tracks when database intervals expire.
- **`Active_Q` (The Ring):** An Implicit Min-Heap or ring buffer (if sufficiently few queries), ordered by `end` coordinate. Tracks active query intervals. (Low depth expected).
- **`overlapsFrame` (The Delta State):** A dense 1D array `[D_sids]u32` representing the current frequency of all active database tracks at the exact sweep line coordinate.
- **`results` (The Accumulator):** A dense 2D matrix `[Q_sids][D_sids]u32` pre-allocated per thread.

---

### 2. The Core Execution Loop (Single-Threaded Chunk Processing)

This loop is executed independently by each thread over a specific memmap chunk and a specific exponential layer.

#### 2.1. The Fast-Forward Escapement (Sparsity Heuristic)

Before advancing the sweep line, the engine checks for batch congestion.

```text
LOOP START:
    IF not nextTile.exists() or nextTile().len():
        1. CLEAR Active_D heap.
        2. ZERO out overlapsFrame.
        3. Jump to next active tile. BINARY SEARCH in-chunk tiles for first index where tileEnd < nextQueryStart.
        5. current_tile = nextActiveTile
```

#### 2.2. The Sweep Line State Machine

The sweep line advances to time $T$, defined as the minimum coordinate among the next four potential events: `D_cursor.start`, `Q_cursor.start`, `peek(Active_D).end`, `peek(Active_Q).end`.

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

To scale to 128+ cores without memory bandwidth starvation, the engine uses a map-reduce pipeline.

**Phase 1: Map (Chunk Processing)**

- A work-stealing scheduler (e.g., Rayon) assigns memmap chunks to idle threads.
- Each thread runs the Core Execution Loop (Section 2) using entirely thread-local data structures, generating a local, fully populated `results` dense matrix. No locks, no atomics, zero false sharing.

**Phase 2: Reduce (Parallel Tree Merge)**

- As threads finish their chunks, the scheduler dynamically pairs them.
- Idle threads take two thread-local `results` matrices and perform a SIMD vectorized matrix addition: `Matrix_A += Matrix_B`.
- This folds $N$ matrices into 1 final global results matrix in $O(\log(\text{threads}))$ steps.

**Phase 3: Sparse Export**

- The final global dense matrix is converted to a sparse format (e.g., COUP/Coordinate List or directly to a TSV/Parquet file) for the end-user, filtering out the millions of zeros. Only the requested statistical overlaps remain.
