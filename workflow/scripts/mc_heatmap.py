"""
Heatmap of chuckle vs Monte Carlo Spearman correlation.

Rows = MC trial count, columns = chuckle p-value threshold (most permissive
to most stringent).  Each cell = Spearman r between chuckle -log10(p) and
MC -log10(p) restricted to beds where chuckle p < threshold.  Grey cells
have fewer than MIN_N beds.

Usage:
    uv run workflow/scripts/mc_heatmap.py \\
        --chuckle data/chuckle/Rheumatoid_arthritis.json \\
        --mc data/montecarlo/Rheumatoid_arthritis.tsv \\
        --trait Rheumatoid_arthritis \\
        -o data/plots/mc_heatmap_Rheumatoid_arthritis.png --no-show
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

_TINY = float(np.finfo(np.float64).tiny)

# Chuckle p-value thresholds: beds with chuckle p < threshold are included.
# Ordered most-permissive to most-stringent.
P_THRESHOLDS = [1.0, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 1e-4, 1e-5]

MIN_N = 5  # minimum beds needed to report a correlation


# ──────────────────────────────────────────────────────────────
# Loaders
# ──────────────────────────────────────────────────────────────


def _strip_name(filename: str) -> str:
    for ext in (".clean.bed.gz", ".bed.gz", ".clean.bed", ".bed"):
        if filename.endswith(ext):
            return filename[: -len(ext)]
    return filename


def load_chuckle(json_path: str) -> dict[str, float]:
    import json
    with open(json_path) as f:
        data = json.load(f)
    results = data.get("results", data) if isinstance(data, dict) else data
    out: dict[str, float] = {}
    for r in results:
        try:
            name = _strip_name(Path(r["db_name"]).name)
            p = max(float(r["p_value"]), _TINY)
            out[name] = -math.log10(p)
        except (KeyError, ValueError):
            pass
    return out


def load_mc_all_milestones(tsv_path: str) -> dict[int, dict[str, float]]:
    """Parse all milestone sections from a MC TSV. Returns {n_trials: {name: -log10(p)}}."""
    milestones: dict[int, dict[str, float]] = {}
    current_trials: int | None = None
    current_data: dict[str, float] = {}

    with open(tsv_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("# milestone="):
                if current_trials is not None:
                    milestones[current_trials] = current_data
                parts = dict(tok.split("=") for tok in line[2:].split())
                current_trials = int(parts.get("trials", 0))
                current_data = {}
                continue
            if not line or line.startswith("db_name"):
                continue
            if current_trials is None:
                continue
            cols = line.split("\t")
            if len(cols) >= 3:
                name = _strip_name(Path(cols[0]).name)
                try:
                    p = max(float(cols[2]), _TINY)
                    current_data[name] = -math.log10(p)
                except ValueError:
                    pass

    if current_trials is not None:
        milestones[current_trials] = current_data

    return milestones


# ──────────────────────────────────────────────────────────────
# Correlation grid
# ──────────────────────────────────────────────────────────────


def compute_grid(
    chuckle: dict[str, float],
    mc_milestones: dict[int, dict[str, float]],
    thresholds: list[float],
) -> tuple[np.ndarray, np.ndarray, list[int]]:
    """
    Returns:
        corr_grid  shape (n_trials, n_thresholds) — Spearman r, NaN where n < MIN_N
        n_grid     shape (n_trials, n_thresholds) — bed count per cell
        sorted_trials
    """
    sorted_trials = sorted(mc_milestones)
    n_t, n_p = len(sorted_trials), len(thresholds)
    corr_grid = np.full((n_t, n_p), np.nan)
    n_grid = np.zeros((n_t, n_p), dtype=int)

    for i, trials in enumerate(sorted_trials):
        mc = mc_milestones[trials]
        common = sorted(set(chuckle) & set(mc))
        for j, threshold in enumerate(thresholds):
            log_thresh = -math.log10(max(threshold, _TINY))
            subset = [n for n in common if chuckle[n] >= log_thresh]
            n = len(subset)
            n_grid[i, j] = n
            if n < MIN_N:
                continue
            xs = np.array([chuckle[n] for n in subset])
            ys = np.array([mc[n] for n in subset])
            r, _ = stats.spearmanr(xs, ys)
            corr_grid[i, j] = float(r)

    return corr_grid, n_grid, sorted_trials


# ──────────────────────────────────────────────────────────────
# Plot
# ──────────────────────────────────────────────────────────────


def plot_heatmap(
    corr_grid: np.ndarray,
    n_grid: np.ndarray,
    sorted_trials: list[int],
    thresholds: list[float],
    trait: str,
    output_path: str | None = None,
    show: bool = True,
) -> None:
    plt.style.use("dark_background")

    n_t, n_p = corr_grid.shape
    fig, ax = plt.subplots(
        figsize=(n_p * 1.25 + 1.5, n_t * 0.9 + 2.5),
        layout="constrained",
    )

    cmap = plt.cm.RdYlGn.copy()
    cmap.set_bad(color="#2a2a2a")

    im = ax.imshow(
        corr_grid,
        aspect="auto",
        cmap=cmap,
        vmin=0,
        vmax=1,
        origin="lower",
    )

    x_labels = []
    for t in thresholds:
        if t >= 0.1:
            x_labels.append(str(t))
        else:
            x_labels.append(f"{t:.0e}")
    ax.set_xticks(range(n_p))
    ax.set_xticklabels(x_labels, rotation=40, ha="right", fontsize=8)

    ax.set_yticks(range(n_t))
    ax.set_yticklabels([f"{t:,}" for t in sorted_trials], fontsize=8)

    ax.set_xlabel("Chuckle p-value threshold  (include beds with p < threshold)", fontsize=9)
    ax.set_ylabel("Monte Carlo trial count", fontsize=9)

    for i in range(n_t):
        for j in range(n_p):
            n = int(n_grid[i, j])
            r = corr_grid[i, j]
            if np.isnan(r):
                cell_text = f"n={n}\n—"
                fg = "#666666"
            else:
                cell_text = f"n={n}\nr={r:.2f}"
                fg = "black" if r > 0.6 else "white"
            ax.text(j, i, cell_text, ha="center", va="center", fontsize=6, color=fg)

    cbar = fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
    cbar.set_label("Spearman r  (chuckle vs MC)", fontsize=9)

    ax.set_title(
        f"{trait.replace('_', ' ')}  |  chuckle vs MC correlation\n"
        "by chuckle p-value threshold and MC trial count",
        fontsize=11,
        fontweight="bold",
        pad=10,
    )

    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=180, bbox_inches="tight")
        print(f"Saved to {output_path}")

    if show:
        plt.show()

    plt.close(fig)


# ──────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--chuckle", required=True, help="Chuckle JSON for the trait")
    parser.add_argument(
        "--mc", required=True, help="MC multi-milestone TSV for the trait"
    )
    parser.add_argument("--trait", required=True, help="Trait name (for plot title)")
    parser.add_argument("-o", "--output", help="Output image path")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()

    chuckle = load_chuckle(args.chuckle)
    mc_milestones = load_mc_all_milestones(args.mc)

    if not mc_milestones:
        raise ValueError(f"No milestone data found in {args.mc}")

    corr_grid, n_grid, sorted_trials = compute_grid(chuckle, mc_milestones, P_THRESHOLDS)

    plot_heatmap(
        corr_grid,
        n_grid,
        sorted_trials,
        P_THRESHOLDS,
        trait=args.trait,
        output_path=args.output,
        show=not args.no_show,
    )


if __name__ == "__main__":
    main()
