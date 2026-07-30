## Chuckle Technical Specification: Analytic Statistics

### Overview

After the sweep phase produces an overlap count matrix, `--stats` mode computes
per-pair significance using an analytic negative-binomial (NB) null model. No
random trials are run at query time. The full pipeline — FFT moment computation,
NB parameter fitting, LLR, and p-value — runs in $O(N \log N)$ at index time and
$O(1)$ per query interval at query time.

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

## 5. NB Null Distribution

### 5.1 Moment Accumulation

For each query interval $i$ of length $l_i$ bins on chromosome $c$, look up the
stored moments for block size $l_i$:

$$\mu_i = \mu(l_i), \qquad \sigma^2_i = \sigma^2(l_i)$$

Genome-wide null moments (independent intervals):

$$\mu_\text{null} = \sum_i \mu_i, \qquad \sigma^2_\text{null} = \sum_i \sigma^2_i$$

### 5.2 NB Parameter Fitting (Method of Moments)

A negative binomial $\operatorname{NB}(r, p)$ has mean $r(1-p)/p$ and variance $r(1-p)/p^2$.
Solving from $\mu_\text{null}$ and $\sigma^2_\text{null}$:

$$
p = \frac{\mu_\text{null}}{\sigma^2_\text{null}}, \qquad
r = \frac{\mu_\text{null}^2}{\sigma^2_\text{null} - \mu_\text{null}}
$$

Requires $\sigma^2_\text{null} > \mu_\text{null} > 0$ (overdispersion).

---

## 6. Log-Likelihood Ratio (Saddlepoint / Exponential Tilt)

To score observation $O$ against $\operatorname{NB}(r, p)$:

$$
p_O = \frac{r}{r + O}, \qquad
\theta = \ln\frac{1 - p_O}{1 - p}, \qquad
\text{LLR} = O\,\theta + r\ln\frac{p_O}{p}
$$

$\text{LLR} = 0$ when $O = \mu_\text{null}$; $\text{LLR} > 0$ for enrichment.

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
   - For each query interval: $O(1)$ moment lookup by direct array index.
   - Accumulate $\mu_\text{null}$, $\sigma^2_\text{null}$; fit NB; compute LLR and p-value.

---

## 9. Output Fields

| Field           | Description                                                                                                                                                                                          |
| --------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `overlap_count` | Number of overlapping interval pairs from the sweep.                                                                                                                                                 |
| `observed_bins` | Exact base-pair overlap in 100 bp bins: $\sum \bigl(\min(q_e, d_e) - \max(q_s, d_s)\bigr) / 100$ summed over all overlapping (query, DB) pairs. Computed during the sweep phase alongside the count. |
| `p_value`       | Right-tailed analytic p-value under the NB null.                                                                                                                                                     |
| `llr`           | Saddlepoint LLR; `null` when NB fit fails.                                                                                                                                                           |
