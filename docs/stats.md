## Chuckle Technical Specification: Analytic Statistics

### Overview

After the sweep phase produces an overlap count matrix, `--stats` mode computes
per-pair significance using an analytic negative-binomial (NB) null model. No
random trials are run and no FFT is performed. The full pipeline — coverage
construction, moment estimation, LLR, and p-value — is closed-form.

---

## 1. Null Model

The null is a **free-relative-gap shift**: each query interval slides
independently over the DB, with no constraint on the relative gaps between
intervals. Under this null, each query interval's position is uniform and
independent on each chromosome; the DB file is fixed with its full spatial
structure.

The test statistic is the **total base-pair overlap** between all query
intervals and the DB file (A), measured in 100 bp bins and summed over all
shared chromosomes.

We do not enumerate positions explicitly. Instead we derive the null distribution
analytically from the statistics of the DB coverage array at the scale of each
query interval type. This captures the true clustering structure of the DB
(the key source of overdispersion) rather than treating each bin independently.

---

## 2. Coverage Representation

At index time, each source file is stored as a **sparse impulse train** — the
derivative of its depth-of-coverage array. Per chromosome, each interval
`[s, e)` contributes:

```
d[s] += (1 − frac_s),   d[s+1] += frac_s        (start spike, fractional)
d[e] -= (1 − frac_e),   d[e+1] -= frac_e        (end spike, fractional)
```

where `frac_x = (x mod 100) / 100` distributes each spike across the two
adjacent 100 bp bins proportional to the sub-bin offset. This anti-aliasing
step ensures sub-bin intervals contribute correctly rather than cancelling to
zero.

The coverage array `g[j]` is recovered at query time by a single prefix sum:

```
g[j] = Σ_{k ≤ j} d[k]
```

`g[j]` is the depth of coverage (number of overlapping intervals) in bin `j`.
Building `g` from the impulse train is O(N_bins) with no FFT.

---

## 3. Per-File Coverage Statistics

For each chromosome, the coverage array `g` of length `N_bins` yields:

```
μ_file = (1/N_bins) Σ_j g[j]               (mean per-bin depth)
σ²_file = (1/N_bins) Σ_j (g[j] − μ_file)²  (variance of per-bin depth)
```

These two numbers are the only per-chromosome statistics that need to be
retained. At query time they are aggregated across shared chromosomes to
parameterise the NB null.

---

## 4. Observed Overlap

The observed base-pair overlap between query A and database file B (in 100 bp
bins) is the **dot product** of their coverage arrays, summed over shared
chromosomes:

```
O = Σ_chrom Σ_j g_A[j] · g_B[j]
```

Equivalently, `O` is the total length (in 100 bp bins) of the intersection of
all pairs of overlapping intervals, accounting for depth multiplicity. It is
computed in O(N_bins) per shared chromosome.

---

## 5. NB Null Distribution

### 5.1 Moment Accumulation

Under the free-relative-gap null, each query interval i of length l_i (bins)
slides independently over the DB.  Its expected overlap with the DB is:

```
μ_i = l_i · μ_B
```

where μ_B is the mean depth per bin of the DB (after masking).

The variance of the overlap for interval i depends on the scale of the interval
relative to the DB's cluster structure.  Let s_i = 2^round(log2(l_i)) be the
nearest power of 2 to l_i in log-space (clamped to [1, 8192] bins).  Define the
**sliding-window variance** at scale s:

```
Var_s(g_B) = (1 / N_windows) · Σ_k [f_s(k) − μ_B·s]²
```

where f_s(k) = Σ_{j=k}^{k+s−1} g_B[j] is the s-bin sliding sum of the DB
coverage array and N_windows = N_bins − s + 1.

Then:

```
σ²_i = Var_{s_i}(g_B)
```

This variance captures the true cluster structure of the DB: when DB intervals
are longer than s, f_s(k) is bimodal (near 0 outside clusters, near s inside),
giving large variance (overdispersion) even for binary coverage.  The NB fit
succeeds whenever DB intervals significantly exceed the query interval scale.

Intervals are summed as independent contributions; chromosomes are also
independent:

```
μ_null  = Σ_i μ_i  = Σ_{chrom,i} l_i · μ_{B,chrom}
σ²_null = Σ_i σ²_i = Σ_{chrom,i} Var_{s_i}(g_{B,chrom})
```

### 5.2 NB Parameter Fitting (Method of Moments)

A negative binomial NB(r, p) has mean `r(1−p)/p` and variance `r(1−p)/p²`.
Solving for `r` and `p` from `μ_null` and `σ²_null`:

```
p = μ_null / σ²_null
r = μ_null² / (σ²_null − μ_null)
```

The NB fit requires `σ²_null > μ_null > 0` (overdispersion). If this
condition is not met — which can occur for highly uniform or very sparse signals
— the NB path returns `None` and no LLR or analytic p-value is reported.

---

## 6. Log-Likelihood Ratio (Saddlepoint / Exponential Tilt)

To score the observation `O` against the NB(r, p) null we compute the
**saddlepoint LLR** — the log-ratio of the tilted distribution (parameterised
to have mean `O`) to the null.

The tilted distribution is NB(r, p_O) where:

```
p_O = r / (r + O)
```

The LLR is then:

```
θ   = ln((1 − p_O) / (1 − p))
LLR = O · θ + r · ln(p_O / p)
```

When `O = μ_null`, `p_O = p` and `θ = 0`, so `LLR = 0`.
When `O > μ_null` (enrichment), `LLR > 0`.
When `O < μ_null` (depletion), `LLR < 0`.

A positive LLR represents evidence for enrichment over the rigid-body-shift null.

---

## 7. Analytic P-Value

The p-value is derived from the LLR via **Wilks' theorem**: under `H₀`, the
statistic `2 · LLR` converges in distribution to χ²(1), i.e., the square of a
standard normal. For the right-tailed enrichment test:

```
p = P(N(0,1) ≥ √(2 · LLR))  =  erfc(√LLR) / 2
```

where `erfc` is the complementary error function. When `LLR ≤ 0` (no
enrichment), the p-value is set to 1.0.

This approximation is tight for large effective sample sizes — i.e., when the
total number of contributing bins is large — and requires no special functions
beyond `erfc`, which is approximated via the Horner polynomial of Abramowitz &
Stegun 7.1.26 (max absolute error < 1.5 × 10⁻⁷).

---

## 8. Positional Filter (Whitelist / Blacklist)

An optional positional filter confines the statistics to a defined accessible
region of the genome, independent of the overlap sweep.

After the prefix sum builds `g[j]`, any bin `j` excluded by the filter is
zeroed:

```
g[j] = 0   for all j not in accessible region
```

This zeroing is applied identically to both the query and DB coverage arrays
before computing moments and observed overlap. The effect is that both the null
distribution and the test statistic are estimated over the same restricted
coordinate space, preventing confounding by inaccessible regions (e.g., ENCODE
blacklist regions, or non-immune tissue for immune-trait GWAS).

- **Whitelist:** Only bins covered by the filter BED are accessible. Chromosomes
  absent from the filter BED are excluded entirely.
- **Blacklist:** Bins covered by the filter BED are zeroed. Chromosomes absent
  from the filter BED remain fully accessible.

---

## 9. Query-Time Pipeline

The stats computation runs as a second phase after the sweep matrix is built:

**Phase A — Load DB depth maps.** The consolidated `depthmap.rkyv` store is
memory-mapped once. For each DB source that has at least one non-zero overlap
with any query, its sparse impulse train is deserialised.

**Phase B — Build query interval data.** For each query file with at least one
non-zero overlap, the BED file is parsed and `build_query_chrom_data` stores the
raw interval lengths (in 100 bp bins) per chromosome.  The mean interval size is
also computed here for use in the observed-overlap heuristic.  Each query's
interval data is cached for reuse across all DB partners.

**Phase C — Compute per-pair statistics (parallel over DB sources).** For each
DB source, `build_chrom_cov_data_with_filter` builds the DB coverage array and
computes `μ_B` and `sliding_vars[0..14]` (sliding-window variance at 14
power-of-2 scales) for each chromosome.  Then, for each query overlapping this
DB source, `compute_analytic_stats` computes:

1. `O` — the heuristic observed overlap: `sweep_count × mean_interval_bins`.
2. `μ_null`, `σ²_null` — independent-interval moment accumulation (Section 5.1).
3. `r`, `p` — NB fit.
4. `LLR` — exponential tilt.
5. `p_value` — Wilks approximation.

DB coverage data is built once per DB source and reused for all overlapping
queries. Total work is O(14 · N_bins · (N_queries + N_db)) instead of the
naive O(N_bins × N_pairs).

**Zero-overlap pairs** receive `p_value = 1.0`, `llr = Some(0.0)`, and
`observed_bins = 0.0` without invoking the NB path. This ensures every
(query, DB) pair appears in the output regardless of whether an overlap was
detected.

---

## 10. Output Fields

Each result record carries:

| Field           | Description                                                            |
| --------------- | ---------------------------------------------------------------------- |
| `overlap_count` | Number of overlapping interval pairs from the sweep.                   |
| `observed_bins` | Base-pair overlap `O` in 100 bp bins (dot product of coverage arrays). |
| `p_value`       | Right-tailed analytic p-value under the NB null.                       |
| `llr`           | Saddlepoint log-likelihood ratio; `null` when NB fit fails.            |

Results are sorted by `p_value` ascending (or by `overlap_count` descending
when `--stats` is not requested).

---
