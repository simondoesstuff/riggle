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
intervals and the DB file (A), measured in 100 bp bins and summed over all
shared chromosomes.

---

## 2. Coverage Representation

At index time, each source file is stored as a **sparse impulse train** — the
derivative of its depth-of-coverage array.  Per chromosome, each interval
`[s, e)` contributes:

```
d[s] += (1 − frac_s),   d[s+1] += frac_s        (start spike, fractional)
d[e] -= (1 − frac_e),   d[e+1] -= frac_e        (end spike, fractional)
```

where `frac_x = (x mod 100) / 100` distributes each spike across the two
adjacent 100 bp bins proportional to the sub-bin offset.

The coverage array `g[j]` is recovered by a single prefix sum:

```
g[j] = Σ_{k ≤ j} d[k]
```

---

## 3. FFT Moment Computation

For each DB chromosome of length N bins, the null distribution moments for
any block size L are computed from three compacted arrays and one FFT:

```
d*[j]  = Σ_{k<j} g[k]              (prefix sum of coverage, size N+1)
d**[j] = Σ_{k<j} d*[k]             (double prefix sum, size N+2)
d2[j]  = Σ_{k<j} (d*[k])²         (prefix sum of squared d*, size N+2)
R[L]   = IFFT(|FFT(d*)|²)[L]       (autocorrelation of d* via FFT)
```

The autocorrelation R is computed via zero-padded FFT (O(N log N)); all other
arrays are O(N) prefix sums.

For any block size L:

```
sum_f(L)  = d**[N+1] − d**[L] − d**[N−L+1]      (sum of window sums)
sum_f2(L) = (d2[N+1] − d2[L]) + d2[N−L+1] − 2·R[L]  (sum of sq. window sums)
n_w       = N − L + 1                              (number of windows)
mean(L)   = sum_f(L)  / n_w
var(L)    = sum_f2(L) / n_w  −  mean(L)²
```

`mean(L)` is the expected overlap (bins) when a block of L bins is placed
uniformly at random; `var(L)` captures the true clustering structure of the DB.

---

## 4. Compacted Moment Storage

All moments are computed once at index time and stored in `momentmap.rkyv`.
Per chromosome, moments for a compacted set of representative block sizes are
stored with a relative error bound ε (default 1%):

- **Dense region** L ∈ [1, ⌈1/ε⌉]: one `(mean, var)` pair per integer L.
- **Sparse region** L > ⌈1/ε⌉: one pair per representative
  L_m = ⌊dense_max · (1+ε)^m⌋ for m = 1, 2, …

A query interval of length L maps to the nearest stored representative within
relative error ε in block size.

Storage per chromosome: O((1/ε) · log_{1+ε}(N/dense_max)) pairs ≈ 1,100
pairs at ε = 0.01 for N = 2.5 M bins. Total: roughly 215 KB per DB file at
single precision.

---

## 5. NB Null Distribution

### 5.1 Moment Accumulation

For each query interval i of length l_i bins on chromosome c, look up the
stored moments for block size l_i on that chromosome:

```
μ_i  = mean(l_i)   (from M[l_i])
σ²_i = var(l_i)    (from M[l_i])
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

## 8. Positional Filter (Whitelist / Blacklist)

An optional positional filter restricts analysis to accessible regions.
When a filter is active, per-chromosome moments are computed on-the-fly at
query time (via the same FFT approach) after zeroing non-accessible bins.
This is done once per DB file (not per query file).

---

## 9. Pipeline Summary

### Index time (per batch of BED files, parallel across files)

1. Parse BED → build sparse impulse train (`depthmap.rkyv`).
2. Build depth signal g → compute d*, d**, d2, R via FFT.
3. Enumerate representative block sizes; store `(mean, var)` pairs in
   `momentmap.rkyv`.

Total index-time cost per chromosome: O(N log N) for the FFT.

### Query time (per (query, DB) pair)

1. Sweep phase: count interval overlaps (unchanged).
2. Stats phase:
   - Load pre-stored moments from `momentmap.rkyv` (no-filter) or compute
     moments on-the-fly with filter applied (O(N log N) per DB file).
   - For each query interval: O(1) moment lookup.
   - Accumulate μ_null, σ²_null; fit NB; compute LLR and p-value.

---

## 10. Output Fields

| Field           | Description                                                   |
| --------------- | ------------------------------------------------------------- |
| `overlap_count` | Number of overlapping interval pairs from the sweep.          |
| `observed_bins` | Heuristic base-pair overlap: `sweep_count × mean_query_bins`. |
| `p_value`       | Right-tailed analytic p-value under the NB null.              |
| `llr`           | Saddlepoint LLR; `null` when NB fit fails.                    |
