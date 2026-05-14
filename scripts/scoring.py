"""
Interval scoring utilities.

Build NCLS indices from BED files and score query sets against database sets
using multiple scoring methods.

Foundational API
----------------
    build_index(bed_path)                    -> Index
    query_overlaps(index, query_path)        -> list[np.ndarray]  (db ids per query interval)
    score_mean_overlaps(query, db)           -> float
    score_bm25_merge(query, db, k1, b)       -> float
    score_bm25_sum(query, db, k1, b)         -> float
    score_fix_size(query, db, n_trials)        -> float

TF in BM25 variants uses Jaccard similarity — size-normalised and independent
of interval length ranges.  bm25-merge collapses stacking db intervals to
coverage before computing Jaccard; bm25-sum sums pairwise Jaccards.

score_mc_fix_size runs a Monte Carlo permutation test: randomise query starts
(preserving per-chrom counts and interval sizes) and return the fraction of
trials whose total overlap count (pairwise query×db pairs) is >= the observed count.

score_mc_fix_size_distinct is the same test but counts distinct db intervals hit
rather than total pairs — lower variance when db intervals stack deeply.

calibrate_mc_stability plots the number of MC trials required for a stable
p-value estimate against the effective genome size.  A neighbourhood radius d
is swept from near-zero to the full chromosome span.  Valid placements at each
d are the union of per-db-interval windows, so the method works correctly even
when hit intervals are scattered (e.g. GWAS SNPs across a chromosome).
Because p ∝ 1/L, required_trials ∝ L, revealing when MC becomes intractable.

CLI
---
    python scoring.py [--method mean|bm25-merge|bm25-sum|mc-fix-size|calibrate] query.bed db1.bed ...
"""

from __future__ import annotations

import gzip
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Generator, Iterator

import numpy as np
from ncls import NCLS

# ---------------------------------------------------------------------------
# BED I/O
# ---------------------------------------------------------------------------


def _bed_lines(path: Path) -> Iterator[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                yield line


def parse_bed(path: str | Path) -> dict[str, list[tuple[int, int]]]:
    """Parse a BED file into chrom -> [(start, end), ...]."""
    result: dict[str, list[tuple[int, int]]] = {}
    for line in _bed_lines(Path(path)):
        parts = line.split("\t") if "\t" in line else line.split()
        chrom, start, end = parts[0], int(parts[1]), int(parts[2])
        result.setdefault(chrom, []).append((start, end))
    return result


# ---------------------------------------------------------------------------
# Interval geometry helpers
# ---------------------------------------------------------------------------


def _merge_length(starts: np.ndarray, ends: np.ndarray) -> int:
    """Total base-pair length of the union of a set of intervals."""
    if len(starts) == 0:
        return 0
    order = np.argsort(starts)
    s, e = starts[order], ends[order]
    total = 0
    cur_s, cur_e = int(s[0]), int(e[0])
    for i in range(1, len(s)):
        if s[i] <= cur_e:
            cur_e = max(cur_e, int(e[i]))
        else:
            total += cur_e - cur_s
            cur_s, cur_e = int(s[i]), int(e[i])
    total += cur_e - cur_s
    return total


def _merge_intervals(starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """Merge overlapping intervals; return (N, 2) int64 array of [start, end) pairs."""
    if len(starts) == 0:
        return np.empty((0, 2), dtype=np.int64)
    order = np.argsort(starts)
    s = starts[order].astype(np.int64)
    e = ends[order].astype(np.int64)
    merged: list[tuple[int, int]] = []
    cur_s, cur_e = int(s[0]), int(e[0])
    for i in range(1, len(s)):
        if s[i] <= cur_e:
            cur_e = max(cur_e, int(e[i]))
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = int(s[i]), int(e[i])
    merged.append((cur_s, cur_e))
    return np.array(merged, dtype=np.int64)


def _sample_from_intervals(merged: np.ndarray, rng: np.random.Generator) -> int:
    """Sample a uniform random integer from the union of merged [start, end) intervals."""
    lengths = merged[:, 1] - merged[:, 0]
    pos = int(rng.integers(0, int(lengths.sum())))
    cum = 0
    for i, ln in enumerate(lengths):
        cum += int(ln)
        if pos < cum:
            return int(merged[i, 0]) + pos - (cum - int(ln))
    return int(merged[-1, 0])  # unreachable in practice


def _pairwise_jaccard(
    q_start: int, q_end: int, db_s: np.ndarray, db_e: np.ndarray
) -> np.ndarray:
    """Jaccard similarity between query interval and each db interval individually."""
    intersection = np.minimum(db_e, q_end) - np.maximum(db_s, q_start)
    union = np.maximum(db_e, q_end) - np.minimum(db_s, q_start)
    return intersection / union


def _jaccard(q_start: int, q_end: int, db_s: np.ndarray, db_e: np.ndarray) -> float:
    """
    Jaccard similarity between query interval [q_start, q_end) and the union
    of a set of overlapping database intervals.

    intersection = bp of q covered by ≥1 db interval
    union        = bp in q or ≥1 db interval
    """
    clipped_s = np.maximum(db_s, q_start)
    clipped_e = np.minimum(db_e, q_end)
    intersection = _merge_length(clipped_s, clipped_e)

    all_s = np.concatenate([[q_start], db_s])
    all_e = np.concatenate([[q_end], db_e])
    union = _merge_length(all_s, all_e)

    return intersection / union if union > 0 else 0.0


# ---------------------------------------------------------------------------
# Index
# ---------------------------------------------------------------------------


@dataclass
class Index:
    """Per-chromosome NCLS indices with precomputed database statistics."""

    ncls: dict[str, NCLS]
    starts: dict[str, np.ndarray]  # db interval starts per chrom
    ends: dict[str, np.ndarray]  # db interval ends per chrom
    total_count: int  # |D|: total intervals across all chroms
    avg_depth: float  # avgdl: sum(lengths) / total_span


def build_index(bed_path: str | Path) -> Index:
    """Build an Index from a BED file."""
    intervals = parse_bed(bed_path)

    ncls_map: dict[str, NCLS] = {}
    starts_map: dict[str, np.ndarray] = {}
    ends_map: dict[str, np.ndarray] = {}

    total_count = 0
    total_length = 0  # sum of (end - start) across all intervals
    total_span = 0  # sum of (max_end - min_start) per chrom

    for chrom, ivs in intervals.items():
        starts = np.array([s for s, _ in ivs], dtype=np.int64)
        ends = np.array([e for _, e in ivs], dtype=np.int64)
        ids = np.arange(len(ivs), dtype=np.int64)

        ncls_map[chrom] = NCLS(starts, ends, ids)
        starts_map[chrom] = starts
        ends_map[chrom] = ends

        total_count += len(ivs)
        total_length += int((ends - starts).sum())
        total_span += int(ends.max() - starts.min())

    avg_depth = total_length / total_span if total_span > 0 else 1.0

    return Index(
        ncls=ncls_map,
        starts=starts_map,
        ends=ends_map,
        total_count=total_count,
        avg_depth=avg_depth,
    )


# ---------------------------------------------------------------------------
# Batch query
# ---------------------------------------------------------------------------


def _grouped_hits(
    ncls: NCLS,
    ivs: list[tuple[int, int]],
) -> Generator[tuple[int, int, int, np.ndarray], None, None]:
    """
    Batch-query *ncls* with all intervals in *ivs* in one call.

    Yields (q_idx, q_start, q_end, db_ids) for each query interval that has
    at least one hit, in ascending q_idx order.
    """
    if not ivs:
        return
    qs = np.array([s for s, _ in ivs], dtype=np.int64)
    qe = np.array([e for _, e in ivs], dtype=np.int64)
    qi = np.arange(len(ivs), dtype=np.int64)

    q_hit, db_hit = ncls.all_overlaps_both(qs, qe, qi)
    if len(q_hit) == 0:
        return

    order = np.argsort(q_hit, kind="stable")
    sq = q_hit[order]
    sd = db_hit[order]

    boundaries = np.flatnonzero(np.diff(sq)) + 1
    for lo, hi in zip(
        np.concatenate([[0], boundaries]),
        np.concatenate([boundaries, [len(sq)]]),
    ):
        i = int(sq[lo])
        yield i, ivs[i][0], ivs[i][1], sd[lo:hi]


# ---------------------------------------------------------------------------
# Query
# ---------------------------------------------------------------------------


def query_overlaps(
    index: Index,
    query_path: str | Path,
) -> list[np.ndarray]:
    """
    For each interval in *query_path*, return an array of overlapping database
    interval IDs (0-based within chrom).  Result is parallel to the query
    intervals (chrom-dict order, then position order within each chrom).
    """
    results: list[np.ndarray] = []
    for chrom, ivs in parse_bed(query_path).items():
        per_q = [np.empty(0, dtype=np.int64)] * len(ivs)
        if chrom in index.ncls:
            for i, _, _, db_ids in _grouped_hits(index.ncls[chrom], ivs):
                per_q[i] = db_ids
        results.extend(per_q)
    return results


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


def score_mean_overlaps(query_path: str | Path, db_path: str | Path) -> float:
    """Mean number of overlapping database intervals per query interval."""
    index = build_index(db_path)
    total_hits = 0
    total_queries = 0
    for chrom, ivs in parse_bed(query_path).items():
        total_queries += len(ivs)
        if chrom not in index.ncls:
            continue
        qs = np.array([s for s, _ in ivs], dtype=np.int64)
        qe = np.array([e for _, e in ivs], dtype=np.int64)
        qi = np.arange(len(ivs), dtype=np.int64)
        _, db_hit = index.ncls[chrom].all_overlaps_both(qs, qe, qi)
        total_hits += len(db_hit)
    return total_hits / total_queries if total_queries > 0 else 0.0


def score_bm25_merge(
    query_path: str | Path,
    db_path: str | Path,
    k1: float = 0.001,
    b: float = 0.1,
) -> float:
    """
    BM25-inspired score of a query BED file against a database BED file.

    Each query interval is a "term" t.  Summed over all query intervals:

        Score = Σ_t  IDF_t · TF_t·(k1+1) / (TF_t + k1·(1 - b + b·|D|/avgdl))

    Where:
        TF_t    = Jaccard(t, union of overlapping db intervals) ∈ [0, 1]
                  size-normalised: intersection / union in base-pair space
        IDF_t   = log((N - df + 0.5) / (df + 0.5)) + 1
                  N = total db interval count (|D|), df = # db intervals hitting t
        |D|     = index.total_count  (total intervals in the db BED file)
        avgdl   = index.avg_depth    (mean genomic depth of the db BED file)
    """
    index = build_index(db_path)
    N = index.total_count
    norm = 1.0 - b + b * (N / index.avg_depth)

    total = 0.0
    for chrom, ivs in parse_bed(query_path).items():
        if chrom not in index.ncls:
            continue
        db_s = index.starts[chrom]
        db_e = index.ends[chrom]
        for _, q_start, q_end, db_ids in _grouped_hits(index.ncls[chrom], ivs):
            df = len(db_ids)
            tf = _jaccard(q_start, q_end, db_s[db_ids], db_e[db_ids])
            idf = math.log((N - df + 0.5) / (df + 0.5)) + 1.0
            total += idf * tf * (k1 + 1.0) / (tf + k1 * norm)

    return total


def score_bm25_sum(
    query_path: str | Path,
    db_path: str | Path,
    k1: float = 0.001,
    b: float = 0.1,
) -> float:
    """
    BM25-inspired score using the sum of pairwise Jaccard similarities as TF.

    Same formula as score_bm25_merge, but TF is unbounded:

        TF_t = Σ_i  Jaccard(t, db_i)   for each overlapping db interval i

    Each pairwise Jaccard is computed without merging, so stacking db intervals
    accumulate additively rather than being collapsed to coverage.
    """
    index = build_index(db_path)
    N = index.total_count
    norm = 1.0 - b + b * (N / index.avg_depth)

    total = 0.0
    for chrom, ivs in parse_bed(query_path).items():
        if chrom not in index.ncls:
            continue
        db_s = index.starts[chrom]
        db_e = index.ends[chrom]
        for _, q_start, q_end, db_ids in _grouped_hits(index.ncls[chrom], ivs):
            df = len(db_ids)
            tf = float(
                _pairwise_jaccard(q_start, q_end, db_s[db_ids], db_e[db_ids]).sum()
            )
            idf = math.log((N - df + 0.5) / (df + 0.5)) + 1.0
            total += idf * tf * (k1 + 1.0) / (tf + k1 * norm)

    return total


def score_mc_fix_size(
    query_path: str | Path,
    db_path: str | Path,
    n_trials: int = 10000,
    rng: np.random.Generator | None = None,
) -> float:
    """
    Monte Carlo p-value: P(overlaps >= observed | random intervals of same sizes).

    Generates *n_trials* randomised query sets — same per-chrom interval sizes
    and counts as the real query, but with uniformly random start positions
    within the database span for each chromosome.  Returns the fraction of
    trials whose total overlap count is >= the observed count.
    """
    if rng is None:
        rng = np.random.default_rng()

    index = build_index(db_path)
    query = parse_bed(query_path)

    # --- observed total overlap count ---
    observed = 0
    for chrom, ivs in query.items():
        if chrom not in index.ncls:
            continue
        qs = np.array([s for s, _ in ivs], dtype=np.int64)
        qe = np.array([e for _, e in ivs], dtype=np.int64)
        qi = np.arange(len(ivs), dtype=np.int64)
        _, db_hit = index.ncls[chrom].all_overlaps_both(qs, qe, qi)
        observed += len(db_hit)

    # --- per-chrom: interval sizes and valid start range [db_min, db_max - size] ---
    chrom_data: dict[str, tuple[np.ndarray, int, int]] = {}
    for chrom, ivs in query.items():
        if chrom not in index.ncls:
            continue
        sizes = np.array([e - s for s, e in ivs], dtype=np.int64)
        db_min = int(index.starts[chrom].min())
        db_max = int(index.ends[chrom].max())
        chrom_data[chrom] = (sizes, db_min, db_max)

    # --- monte carlo ---
    at_least = 0
    for _ in range(n_trials):
        trial_count = 0
        for chrom, (sizes, db_min, db_max) in chrom_data.items():
            # upper bound for each start: clamp so interval stays within db span
            hi = np.maximum(db_max - sizes, db_min) + 1  # +1: rng.integers is exclusive
            rand_starts = rng.integers(db_min, hi)
            rand_ends = rand_starts + sizes
            qi = np.arange(len(sizes), dtype=np.int64)
            _, db_hit = index.ncls[chrom].all_overlaps_both(
                rand_starts.astype(np.int64),
                rand_ends.astype(np.int64),
                qi,
            )
            trial_count += len(db_hit)
        if trial_count >= observed:
            at_least += 1

    return at_least / n_trials


def score_mc_fix_size_distinct(
    query_path: str | Path,
    db_path: str | Path,
    n_trials: int = 10000,
    rng: np.random.Generator | None = None,
) -> float:
    """
    Monte Carlo p-value using distinct db intervals hit as the test statistic.

    Same permutation scheme as score_mc_fix_size, but instead of counting total
    (query, db) overlap pairs, counts the number of distinct db interval IDs
    touched by any query interval.  This reduces trial-to-trial variance when db
    intervals stack deeply (e.g. SNP queries against high-coverage annotation
    tracks), where a single lucky placement can add hundreds of pairs.
    """
    if rng is None:
        rng = np.random.default_rng()

    index = build_index(db_path)
    query = parse_bed(query_path)

    # --- observed: distinct db intervals hit ---
    observed = 0
    for chrom, ivs in query.items():
        if chrom not in index.ncls:
            continue
        qs = np.array([s for s, _ in ivs], dtype=np.int64)
        qe = np.array([e for _, e in ivs], dtype=np.int64)
        qi = np.arange(len(ivs), dtype=np.int64)
        _, db_hit = index.ncls[chrom].all_overlaps_both(qs, qe, qi)
        observed += len(np.unique(db_hit))

    # --- per-chrom: interval sizes and valid start range ---
    chrom_data: dict[str, tuple[np.ndarray, int, int]] = {}
    for chrom, ivs in query.items():
        if chrom not in index.ncls:
            continue
        sizes = np.array([e - s for s, e in ivs], dtype=np.int64)
        db_min = int(index.starts[chrom].min())
        db_max = int(index.ends[chrom].max())
        chrom_data[chrom] = (sizes, db_min, db_max)

    # --- monte carlo ---
    at_least = 0
    for _ in range(n_trials):
        trial_count = 0
        for chrom, (sizes, db_min, db_max) in chrom_data.items():
            hi = np.maximum(db_max - sizes, db_min) + 1
            rand_starts = rng.integers(db_min, hi)
            rand_ends = rand_starts + sizes
            qi = np.arange(len(sizes), dtype=np.int64)
            _, db_hit = index.ncls[chrom].all_overlaps_both(
                rand_starts.astype(np.int64),
                rand_ends.astype(np.int64),
                qi,
            )
            trial_count += len(np.unique(db_hit))
        if trial_count >= observed:
            at_least += 1

    return at_least / n_trials


# ---------------------------------------------------------------------------
# MC stability calibration
# ---------------------------------------------------------------------------


def calibrate_mc_stability(
    query_path: str | Path,
    db_path: str | Path,
    *,
    n_points: int = 10,
    n_trials: int = 5000,
    target_cv: float = 0.2,
    output: str | Path | None = None,
    rng: np.random.Generator | None = None,
) -> None:
    """
    Plot estimated MC trials required for stable p-values vs effective genome size.

    Tests the hypothesis that MC noise in SNP-based analyses arises because the
    overlap probability p ∝ 1/L shrinks as the valid placement space L grows, so
    the required number of trials ∝ L, quickly becoming intractable.

    Instead of a bounding-box expansion (which collapses when hits are spread
    across a chromosome), a neighbourhood radius d is swept.  For each d, a
    query interval of size s is valid at any position within d bp of *any*
    individual db interval — i.e. the union of per-db-interval windows:

        valid starts for size s = ⋃_i  [db_s_i − d − s + 1,  db_e_i + d)

    d is swept log-uniformly from 0 (only overlapping placements) to the span
    of the largest db chromosome (effectively the whole genome).  The effective
    genome size plotted on the x-axis is the total bp in valid start zones
    (averaged over query sizes per chromosome).

    At each d, n_trials random placements are drawn and p̂ estimated; then:

        n_required = (1 − p̂) / (p̂ · target_cv²)

    Parameters
    ----------
    query_path, db_path
        Paths to BED files.
    n_points
        Number of neighbourhood radii to sample (log-spaced, default 10).
    n_trials
        MC trials per radius used to estimate p (default 5000).
    target_cv
        Target coefficient of variation for p-value stability (default 0.2).
    output
        Path for the saved plot (PNG/PDF); show interactively if None.
    rng
        Numpy random Generator (optional).
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise ImportError("matplotlib is required for calibration: pip install matplotlib")

    if rng is None:
        rng = np.random.default_rng()

    index = build_index(db_path)
    query = parse_bed(query_path)

    # -- per-chrom data: query sizes, db starts/ends -----------------------------
    chrom_info: dict[str, np.ndarray] = {}   # chrom -> query sizes
    for chrom, ivs in query.items():
        if chrom not in index.ncls:
            continue
        qs = np.array([s for s, _ in ivs], dtype=np.int64)
        qe = np.array([e for _, e in ivs], dtype=np.int64)
        qi = np.arange(len(ivs), dtype=np.int64)
        _, db_hit = index.ncls[chrom].all_overlaps_both(qs, qe, qi)
        if len(db_hit) == 0:
            continue
        chrom_info[chrom] = np.array([e - s for s, e in ivs], dtype=np.int64)

    if not chrom_info:
        print("No overlapping chromosomes found; nothing to calibrate.")
        return

    # -- observed distinct db intervals hit --------------------------------------
    observed = 0
    for chrom in chrom_info:
        ivs = query[chrom]
        qs = np.array([s for s, _ in ivs], dtype=np.int64)
        qe = np.array([e for _, e in ivs], dtype=np.int64)
        qi = np.arange(len(ivs), dtype=np.int64)
        _, db_hit = index.ncls[chrom].all_overlaps_both(qs, qe, qi)
        observed += len(np.unique(db_hit))

    # -- d sweep range: 0 to span of largest db chromosome ----------------------
    max_d = max(
        int(index.ends[c].max() - index.starts[c].min()) for c in chrom_info
    )
    # d=0 means query must overlap db interval; start at 1 to avoid empty zones
    d_values = np.unique(
        np.logspace(0, math.log10(max_d), n_points).astype(np.int64)
    )

    # -- full-db effective size (for reference line on plot) ---------------------
    full_eff_size = sum(
        int(index.ends[c].max() - index.starts[c].min()) for c in chrom_info
    )

    span_bp: list[float] = []
    p_hats: list[float] = []
    n_required: list[float] = []
    resolved_flags: list[bool] = []

    for d in d_values:
        # -- build per-chrom valid-start zone maps (keyed by unique query size) --
        # For size s, valid starts = ⋃_i [db_s_i − d − s + 1, db_e_i + d)
        per_chrom: list[tuple[str, np.ndarray, dict[int, np.ndarray]]] = []
        total_eff_size = 0

        for chrom, sizes in chrom_info.items():
            db_s = index.starts[chrom]
            db_e = index.ends[chrom]
            size_zones: dict[int, np.ndarray] = {}
            for s in np.unique(sizes):
                zone_starts = db_s - int(d) - int(s) + 1
                zone_ends = db_e + int(d)
                size_zones[int(s)] = _merge_intervals(zone_starts, zone_ends)
            # effective size: use mean query size as representative
            mean_s = int(sizes.mean())
            rep_zones = size_zones[min(size_zones, key=lambda k: abs(k - mean_s))]
            total_eff_size += int((rep_zones[:, 1] - rep_zones[:, 0]).sum())
            per_chrom.append((chrom, sizes, size_zones))

        # -- MC estimate of p for this d -----------------------------------------
        at_least = 0
        for _ in range(n_trials):
            trial_count = 0
            for chrom, sizes, size_zones in per_chrom:
                rand_starts = np.empty(len(sizes), dtype=np.int64)
                for j, s in enumerate(sizes):
                    rand_starts[j] = _sample_from_intervals(size_zones[int(s)], rng)
                rand_ends = rand_starts + sizes
                qi = np.arange(len(sizes), dtype=np.int64)
                _, db_hit = index.ncls[chrom].all_overlaps_both(rand_starts, rand_ends, qi)
                trial_count += len(np.unique(db_hit))
            if trial_count >= observed:
                at_least += 1

        # resolved = at least one trial succeeded; unresolved = zero hits (p unknown)
        resolved = at_least > 0
        if resolved:
            p_hat = at_least / n_trials
            n_req = (1.0 - p_hat) / (p_hat * target_cv**2)
            status = f"p̂={p_hat:.4f}  n_req={n_req:,.0f}"
        else:
            # p < 1/n_trials; n_req is a strict lower bound of (n_trials-1)/target_cv²
            p_hat = 1.0 / n_trials          # upper bound on p
            n_req = (1.0 - p_hat) / (p_hat * target_cv**2)  # lower bound on n_req
            status = f"p̂<{1/n_trials:.2e}  n_req>{n_req:,.0f}  [unresolved]"

        span_bp.append(float(total_eff_size))
        p_hats.append(p_hat)
        n_required.append(n_req)
        resolved_flags.append(resolved)

        print(f"  d={d:>10,} bp  eff={total_eff_size / 1e6:>8.2f} Mbp  {status}")

    # -- plot --------------------------------------------------------------------
    span_mbp = [s / 1e6 for s in span_bp]
    full_mbp = full_eff_size / 1e6

    res = [i for i, r in enumerate(resolved_flags) if r]
    unres = [i for i, r in enumerate(resolved_flags) if not r]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # required-trials panel
    if res:
        ax1.loglog(
            [span_mbp[i] for i in res], [n_required[i] for i in res],
            "o-", color="steelblue", label="resolved",
        )
    if unres:
        ax1.scatter(
            [span_mbp[i] for i in unres], [n_required[i] for i in unres],
            marker=6,  # caretup = upward triangle = lower bound
            color="steelblue", s=80, zorder=5, label="lower bound (unresolved)",
        )
        # dashed line connecting last resolved point (if any) to first lower bound
        if res:
            join_x = [span_mbp[res[-1]], span_mbp[unres[0]]]
            join_y = [n_required[res[-1]], n_required[unres[0]]]
            ax1.loglog(join_x, join_y, "--", color="steelblue", alpha=0.4)
    ax1.axhline(
        n_trials, color="red", linestyle="--", alpha=0.7,
        label=f"current n_trials={n_trials:,}",
    )
    ax1.axvline(full_mbp, color="grey", linestyle=":", alpha=0.6, label="full db span")
    ax1.set_xlabel("Effective genome size (Mbp)")
    ax1.set_ylabel(f"Required trials (CV < {target_cv})")
    ax1.set_title("MC trials required vs genome size")
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3, which="both")

    # p-value panel
    if res:
        ax2.loglog(
            [span_mbp[i] for i in res], [p_hats[i] for i in res],
            "s-", color="darkorange", label="resolved",
        )
    if unres:
        ax2.scatter(
            [span_mbp[i] for i in unres], [p_hats[i] for i in unres],
            marker=7,  # caretdown = upper bound
            color="darkorange", s=80, zorder=5, label="upper bound (unresolved)",
        )
        if res:
            join_x = [span_mbp[res[-1]], span_mbp[unres[0]]]
            join_y = [p_hats[res[-1]], p_hats[unres[0]]]
            ax2.loglog(join_x, join_y, "--", color="darkorange", alpha=0.4)
    ax2.axvline(full_mbp, color="grey", linestyle=":", alpha=0.6, label="full db span")
    ax2.set_xlabel("Effective genome size (Mbp)")
    ax2.set_ylabel("Estimated p-value")
    ax2.set_title("p-value vs genome size")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3, which="both")

    fig.suptitle(
        f"MC stability calibration — "
        f"query: {Path(query_path).name}  db: {Path(db_path).name}\n"
        f"observed distinct db hits: {observed}  "
        f"n_trials per point: {n_trials:,}",
        fontsize=10,
    )
    fig.tight_layout()

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
        print(f"Saved: {output}")
    else:
        plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

_METHODS = {
    "mean": score_mean_overlaps,
    "bm25-merge": score_bm25_merge,
    "bm25-sum": score_bm25_sum,
    "mc-fix-size": score_mc_fix_size,
    "mc-fix-size-distinct": score_mc_fix_size_distinct,
}

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Score a query BED file against one or more database BED files."
    )
    parser.add_argument("query", help="Query BED file (.bed or .bed.gz)")
    parser.add_argument("databases", nargs="+", help="Database BED file(s)")
    parser.add_argument(
        "--method",
        choices=[*_METHODS, "calibrate", "all"],
        default="mean",
        help="Scoring method, 'calibrate' to plot MC stability, or 'all' (default: mean)",
    )
    parser.add_argument(
        "--n-trials",
        type=int,
        default=1000,
        help="MC trials for mc-fix-size methods (default: 1000); trials per span point for calibrate (default overridden to 5000)",
    )
    parser.add_argument(
        "--n-points",
        type=int,
        default=10,
        help="Number of span sizes to sweep for calibrate (default: 10)",
    )
    parser.add_argument(
        "--target-cv",
        type=float,
        default=0.2,
        help="Target CV for stability estimate in calibrate (default: 0.2)",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output path for calibrate plot (PNG/PDF); show interactively if omitted",
    )
    args = parser.parse_args()

    if args.method == "calibrate":
        n_trials_calibrate = args.n_trials if args.n_trials != 1000 else 5000
        for db in args.databases:
            print(f"calibrate: {db}")
            calibrate_mc_stability(
                args.query,
                db,
                n_points=args.n_points,
                n_trials=n_trials_calibrate,
                target_cv=args.target_cv,
                output=args.output,
            )
    else:
        methods = _METHODS if args.method == "all" else {args.method: _METHODS[args.method]}
        for name, score_fn in methods.items():
            print(name)
            for db in args.databases:
                if name in ("mc-fix-size", "mc-fix-size-distinct"):
                    score = score_fn(args.query, db, n_trials=args.n_trials)
                else:
                    score = score_fn(args.query, db)
                print(f"{db}\t{score:.6f}")
