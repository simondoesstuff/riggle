## Chuckle Technical Specification: Analytic Statistics

### Overview

After the sweep phase produces an overlap count matrix, `--stats` mode computes
per-pair significance using an analytic negative-binomial (NB) null model. No
random trials are run at query time. The full pipeline — FFT moment computation,
NB parameter fitting, LLR, and p-value — runs in O(N log N) at index time and
O(1) per query interval at query time.

---

## 1. Null Model

The null is a **free-relative-gap shift**: each query interval slides
independently over the DB with a uniform random position on each chromosome.
The DB file is fixed with its full spatial structure.

The test statistic is the **total base-pair overlap** between all query
intervals and the DB file, measured in 100 bp bins and summed over all
shared chromosomes.

---

## 2. Coverage Representation

At index time, each source file is stored as a **dense signed impulse train** —
the derivative of its depth-of-coverage array.  Per chromosome, each interval
`[s, e)` contributes:

```
d[s/100]   += (1 − frac_s)      (fractional start spike, left bin)
d[s/100+1] += frac_s             (fractional start spike, right bin)
d[e/100]   -= (1 − frac_e)      (fractional end spike, left bin)
d[e/100+1] -= frac_e             (fractional end spike, right bin)
```

where `frac_x = (x mod 100) / 100` distributes each spike proportionally
across the two adjacent 100 bp bins.

The coverage array `g[j]` is recovered by a SIMD-accelerated prefix sum:

```
g[j] = Σ_{k ≤ j} d[k]
```

The depth signal g[] is transient — only the moment tables (below) are
persisted to `momentmap.rkyv`.

---

## 3. FFT Moment Computation

For each DB chromosome of length N bins, the null distribution moments for
all block sizes L are computed from three compacted arrays and one FFT:

```
d*[j]  = Σ_{k<j} g[k]              (prefix sum of coverage, size N+1)
d**[j] = Σ_{k<j} d*[k]             (double prefix sum, size N+2)
d2[j]  = Σ_{k<j} (d*[k])²         (prefix sum of squared d*, size N+2)
                                     computed with SIMD squaring (f64x4)
R[L]   = IFFT(|FFT(d*)|²)[L]       (autocorrelation of d* via zero-padded FFT)
```

The autocorrelation R is computed via zero-padded FFT (O(N log N)); all other
arrays are O(N) prefix sums.

For any block size L:

```
sum_f(L)  = d**[N+1] − d**[L] − d**[N−L+1]      (sum of window sums)
sum_f2(L) = (d2[N+1] − d2[L]) + d2[N−L+1] − 2·R[L]
n_w       = N − L + 1                              (number of windows)
mean(L)   = sum_f(L)  / n_w
var(L)    = sum_f2(L) / n_w  −  mean(L)²
```

`mean(L)` is the expected overlap (bins) when a block of L bins is placed
uniformly at random; `var(L)` captures the true clustering structure of the DB.

---

## 4. Compact Moment Storage

Only a subset of block sizes L are stored, tolerating at most 1% relative
error in L (eps = 0.01, T = floor(1/eps) = 100):

```
L = 1..T        : stored exactly (every block size, no approximation)
L > T           : only the first L in each log_{1+eps}-spaced slot is stored
                  slot k covers L in [T·(1+eps)^k, T·(1+eps)^{k+1})
```

**Lookup formula** — given a query block size L, the 0-indexed position in the
compact store is:

```
compact_index(L) = L − 1                                 if L ≤ T
compact_index(L) = T + floor(log_{1.01}(L/T)) − 1       if L > T
```

O(1) lookup via a single `ln` and integer cast; no stored L values needed.

**Storage per chromosome**: approximately T + ceil(log_{1.01}(N/T)) ≈ 1100
entries for chr1 (N ≈ 2.49 M), versus 2.49 M for the dense scheme.
Storage ≈ 2 × 1100 × 4 bytes ≈ 9 KB per chromosome, ≈ 220 KB total for hg38
per DB source (vs 250 MB dense).

The maximum approximation error is eps × L in block size, which translates to
a negligible shift in the returned `(mean, var)` for smoothly varying moment
functions.

---

## 5. NB Null Distribution

### 5.1 Moment Accumulation

For each query interval i of length l_i bins on chromosome c, look up the
stored moments for block size l_i on that chromosome:

```
μ_i  = mean(l_i)   (from momentmap.rkyv)
σ²_i = var(l_i)    (from momentmap.rkyv)
```

Genome-wide null moments (independent intervals):

```
μ_null  = Σ_i μ_i
σ²_null = Σ_i σ²_i
```

### 5.2 NB Parameter Fitting (Method of Moments)

A negative binomial NB(r, p) has mean `r(1−p)/p` and variance `r(1−p)/p²`.
Solving from `μ_null` and `σ²_null`:

```
p = μ_null / σ²_null
r = μ_null² / (σ²_null − μ_null)
```

Requires `σ²_null > μ_null > 0` (overdispersion).

---

## 6. Log-Likelihood Ratio (Saddlepoint / Exponential Tilt)

To score observation `O` against NB(r, p):

```
p_O = r / (r + O)
θ   = ln((1 − p_O) / (1 − p))
LLR = O · θ + r · ln(p_O / p)
```

`LLR = 0` when `O = μ_null`; `LLR > 0` for enrichment.

---

## 7. Analytic P-Value

Via Wilks' theorem: `2·LLR ~ χ²(1)`, so:

```
p = erfc(√LLR) / 2
```

---

## 8. Pipeline Summary

### Index time (per batch of BED files, parallel across files)

1. Parse BED → assemble dense signed impulse train per chromosome.
2. SIMD prefix sum (f32x8 two-pass scan) → depth signal g[].
3. Build d*, d**, d2 (SIMD f64x4 for squaring), autocorrelation R via FFT.
4. Compute mean(L) and var(L) for ~T + log_{1.01}(N/T) ≈ 1100 sampled L values.
5. Store compact (mean, var) pairs as f32 in `momentmap.rkyv`.

Total per chromosome: O(N log N) for the FFT; O(N) for prefix sums; O(log N)
for moment sampling and storage.

### Query time (per (query, DB) pair)

1. Sweep phase: count interval overlaps (unchanged).
2. Stats phase:
   - Load pre-stored moments from `momentmap.rkyv` (memmap'd, O(1) per lookup).
   - For each query interval: O(1) moment lookup by direct array index.
   - Accumulate μ_null, σ²_null; fit NB; compute LLR and p-value.

---

## 9. Output Fields

| Field           | Description                                                   |
| --------------- | ------------------------------------------------------------- |
| `overlap_count` | Number of overlapping interval pairs from the sweep.          |
| `observed_bins` | Heuristic base-pair overlap: `sweep_count × mean_query_bins`. |
| `p_value`       | Right-tailed analytic p-value under the NB null.              |
| `llr`           | Saddlepoint LLR; `null` when NB fit fails.                    |
