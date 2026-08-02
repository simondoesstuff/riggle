## Chuckle Technical Specification: Analytic Statistics

### Overview

After the sweep phase produces an overlap count matrix, `--stats` mode computes
per-pair significance using an analytic null model — negative-binomial (NB) when the
null is overdispersed, Poisson on count units when it is not. No random trials are run
at query time. The full pipeline runs in $O(N \log N)$ at index time and $O(1)$ per
query interval at query time.

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
the derivative of its depth-of-coverage array. Per chromosome, each interval
$[s, e)$ contributes:

$$
\begin{aligned}
d\!\left[\lfloor s/100 \rfloor\right]     &\mathrel{+}= 1 - f_s \\
d\!\left[\lfloor s/100 \rfloor + 1\right] &\mathrel{+}= f_s \\
d\!\left[\lfloor e/100 \rfloor\right]     &\mathrel{-}= 1 - f_e \\
d\!\left[\lfloor e/100 \rfloor + 1\right] &\mathrel{-}= f_e
\end{aligned}
$$

where $f_x = (x \bmod 100) / 100$ distributes each spike proportionally
across the two adjacent 100 bp bins.

The coverage array $g[j]$ is recovered by a SIMD-accelerated prefix sum:

$$g[j] = \sum_{k \le j} d[k]$$

The depth signal $g[\cdot]$ is transient — only the moment tables (below) are
persisted to `momentmap.rkyv`.

---

## 3. FFT Moment Computation

For each DB chromosome of length $N$ bins, the null distribution moments for
all block sizes $L$ are computed from three compacted arrays and one FFT:

$$
\begin{aligned}
d^{*}[j]   &= \sum_{k < j} g[k]               & &(N+1 \text{ entries; prefix sum of coverage}) \\
d^{**}[j]  &= \sum_{k < j} d^{*}[k]           & &(N+2 \text{ entries; double prefix sum}) \\
d^{2}[j]   &= \sum_{k < j} \bigl(d^{*}[k]\bigr)^2 & &(N+2 \text{ entries; prefix sum of squared } d^{*}) \\
R[L]       &= \operatorname{IFFT}\!\bigl(|\operatorname{FFT}(d^{*})|^2\bigr)[L] & &(\text{autocorrelation of } d^{*} \text{ via zero-padded FFT})
\end{aligned}
$$

The autocorrelation $R$ is computed via zero-padded FFT in $O(N \log N)$; all other
arrays are $O(N)$ prefix sums.

For any block size $L$, defining $n_w = N - L + 1$ (number of windows):

$$
\begin{aligned}
S_1(L)      &= d^{**}[N+1] - d^{**}[L] - d^{**}[N-L+1] \\
S_2(L)      &= \bigl(d^{2}[N+1] - d^{2}[L]\bigr) + d^{2}[N-L+1] - 2R[L] \\
\mu(L)      &= S_1(L) \;/\; n_w \\
\sigma^2(L) &= S_2(L) \;/\; n_w \;-\; \mu(L)^2
\end{aligned}
$$

$\mu(L)$ is the expected overlap (bins) when a block of $L$ bins is placed
uniformly at random; $\sigma^2(L)$ captures the true clustering structure of the DB.

---

## 4. Compact Moment Storage

Only a subset of block sizes $L$ are stored, tolerating at most 1% relative
error in $L$ ($\varepsilon = 0.01$, $T = \lfloor 1/\varepsilon \rfloor = 100$):

- $L = 1, \ldots, T$: stored exactly, no approximation.
- $L > T$: only the first $L$ in each $\log_{1+\varepsilon}$-spaced slot is stored;
  slot $k$ covers $L \in \bigl[T(1+\varepsilon)^k,\; T(1+\varepsilon)^{k+1}\bigr)$.

**Lookup formula** — given a query block size $L$, the 0-indexed position in the
compact store is:

$$
\text{compact\_index}(L) = \begin{cases}
L - 1 & L \le T \\[4pt]
T + \left\lfloor \log_{1.01}(L/T) \right\rfloor - 1 & L > T
\end{cases}
$$

$O(1)$ lookup via a single $\ln$ and integer cast; no stored $L$ values needed.

**Storage per chromosome**: approximately $T + \lceil \log_{1.01}(N/T) \rceil \approx 1100$
entries for chr1 ($N \approx 2.49\,\text{M}$), versus 2.49 M for the dense scheme.
Storage $\approx 2 \times 1100 \times 4\ \text{bytes} \approx 9\ \text{KB}$ per chromosome,
$\approx 220\ \text{KB}$ total for hg38 per DB source (vs 250 MB dense).

The maximum approximation error is $\varepsilon L$ in block size, which translates to
a negligible shift in the returned $(\mu, \sigma^2)$ for smoothly varying moment functions.

---

## 5. Null Distribution and Dispatch

### 5.1 Moment Accumulation

For each query interval $i$ of length $l_i$ bins, let $l_\text{int} = \lceil l_i \rceil$ and
$s_i = l_i / l_\text{int}$ (scale factor; $s_i = 1$ for full-bin intervals, $s_i < 1$ for sub-bin).
Look up stored moments at $l_\text{int}$:

$$
\mu_i = \mu(l_\text{int}) \cdot s_i, \qquad
\sigma^2_i = \sigma^2(l_\text{int}) \cdot s_i^2, \qquad
\mu^\text{raw}_i = \mu(l_\text{int}) / l_\text{int}
$$

$\mu^\text{raw}_i$ approximates the per-interval overlap probability (coverage fraction);
it is scale-free and equals $\mu(1)$ for $l_i \le 1$.

Genome-wide sums:

$$
\mu = \sum_i \mu_i, \qquad
\sigma^2 = \sum_i \sigma^2_i, \qquad
\mu^\text{raw} = \sum_i \mu^\text{raw}_i
$$

### 5.2 Dispatch on Dispersion

**Overdispersed** ($\sigma^2 > \mu$): NB null on the bin-overlap statistic $O$.
This holds for $L \ge 2$ in the RME index (95–100% of pairs).

**Underdispersed** ($\sigma^2 \le \mu$): Poisson null on count units.
The effective scale $s = \mu / \mu^\text{raw}$ (≈ mean interval length) maps $O$ to
count units, recovering the full-power count-based LLR:

$$
\hat{O} = O / s, \qquad \mu^\text{count} = \mu^\text{raw}
$$

This is exact for uniform-length queries (e.g. all 1 bp SNPs, where $s = 0.01$).

### 5.3 NB Parameter Fitting (Method of Moments)

A negative binomial $\operatorname{NB}(r, p)$ has mean $r(1-p)/p$ and variance $r(1-p)/p^2$.
Solving from $\mu$ and $\sigma^2$:

$$
p = \frac{\mu}{\sigma^2}, \qquad r = \frac{\mu^2}{\sigma^2 - \mu}
$$

---

## 6. Log-Likelihood Ratios

**NB saddlepoint** (overdispersed), observation $O$ against $\operatorname{NB}(r, p)$:

$$
p_O = \frac{r}{r + O}, \qquad
\theta = \ln\frac{1 - p_O}{1 - p}, \qquad
\text{LLR}_\text{NB} = O\,\theta + r\ln\frac{p_O}{p}
$$

**Poisson saddlepoint** (underdispersed), observation $\hat{O}$ against $\operatorname{Poisson}(\mu^\text{count})$:

$$\text{LLR}_\text{Poisson} = \hat{O}\ln(\hat{O}/\mu^\text{count}) - (\hat{O} - \mu^\text{count})$$

Both LLRs equal 0 at the null mean and are positive for enrichment.

---

## 7. Analytic P-Value

Via Wilks' theorem, $2\cdot\text{LLR} \sim \chi^2(1)$, so:

$$p = \frac{\operatorname{erfc}\!\left(\sqrt{\text{LLR}}\right)}{2}$$

---

## 8. Pipeline Summary

### Index time (per batch of BED files, parallel across files)

1. Parse BED → assemble dense signed impulse train per chromosome.
2. SIMD prefix sum (`f32x8` two-pass scan) → depth signal $g[\cdot]$.
3. Build $d^{*}$, $d^{**}$, $d^{2}$ (SIMD `f64x4` for squaring), autocorrelation $R$ via FFT.
4. Compute $\mu(L)$ and $\sigma^2(L)$ for ${\sim}T + \log_{1.01}(N/T) \approx 1100$ sampled $L$ values.
5. Store compact $(\mu, \sigma^2)$ pairs as `f32` in `momentmap.rkyv`.

Total per chromosome: $O(N \log N)$ for the FFT; $O(N)$ for prefix sums; $O(\log N)$
for moment sampling and storage.

### Query time (per (query, DB) pair)

1. Sweep phase: count interval overlaps (unchanged).
2. Stats phase:
   - Load pre-stored moments from `momentmap.rkyv` (memmap'd, $O(1)$ per lookup).
   - Group query intervals by unique block size $L$ per chromosome; one moment lookup per
     unique $(L, \text{chrom})$ pair — $O(N_\text{unique})$ rather than $O(N_\text{intervals})$.
   - Accumulate $\mu$, $\sigma^2$, $\mu^\text{raw}$; dispatch on dispersion (§5.2); compute LLR and p-value.
   - **Zero-overlap pairs** (no shared overlapping intervals) always satisfy $O = 0$, giving
     $\text{LLR} = 0$ and $p = 1$ trivially.  These pairs are **not emitted** in the output;
     callers treat any absent (query, DB) entry as $p = 1$.

---

## 9. Output Fields

| Field           | Description                                                                                                                                                                                          |
| --------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `overlap_count` | Number of overlapping interval pairs from the sweep.                                                                                                                                                 |
| `observed_bins` | Exact base-pair overlap in 100 bp bins: $\sum \bigl(\min(q_e, d_e) - \max(q_s, d_s)\bigr) / 100$ summed over all overlapping (query, DB) pairs. Computed during the sweep phase alongside the count. |
| `p_value`       | Right-tailed analytic p-value: NB saddlepoint when overdispersed, Poisson saddlepoint on count units when underdispersed.  Only present for pairs with at least one overlapping interval; absent pairs have $p = 1$ implicitly. |
| `llr`           | Saddlepoint LLR (NB or Poisson); `null` only when no query intervals share a chromosome with the DB.                                                                                                 |
