# Deferred Work

## f32 FFT via g-autocorrelation reformulation

**Motivation.** The FFT dominates index-time cost (~58% of worker samples). Switching the
autocorrelation from f64 to f32 would halve memory traffic for the 8M-point chr1 transform
(128 MB → 64 MB per thread) and double NEON arithmetic throughput (4 f32 lanes vs 2 f64),
giving an estimated **~1.4× overall speedup** (~29% wall-time reduction).

**Why a type-swap doesn't work.** The current formula autocorrelates `d*` — the prefix sum
of the coverage signal `g`. For chr1 at typical coverage (~5%), `d*[N] ≈ 125 000`, so the
variance formula involves near-cancellation of terms ~10¹⁶:

```
sum_f2(L) = (d2[N+1] − d2[L]) + d2[N−L+1] − 2·R[L]
```

Condition number ~10⁸–10⁹; f32 machine epsilon (~6×10⁻⁸) leaves zero correct digits in
σ²(L). The f64 widening is load-bearing.

**The path forward.** Autocorrelate the raw signal `g` (values in [0, 1]) instead.
Let `A[m] = Σ_k g[k]·g[k+m]`. Then:

```
sum_f2(L) = Σ_{m=−(L−1)}^{L−1} (L − |m|) · A[m]
```

Since `g` is in [0, 1], `A[m] ≤ Σ g ≈ 125 000`, condition number drops to O(10) — f32-safe.
The per-L O(1) recovery from `A` needs two prefix sums: one over `A[m]` and one over
`m·A[m]`.

**Before implementing:**

1. Prove the identity algebraically, including exact boundary terms (finite valid-window
   range differs from the infinite-signal convolution).
2. Unit-test the reformulated `sum_f2` against the current `d*`-based path on at least one
   real chromosome with edge cases (sparse file, dense file, sub-bin intervals at both ends).
3. Verify that `g` can be passed as f32 directly to `realfft::RealFftPlanner::<f32>` without
   additional preprocessing.
