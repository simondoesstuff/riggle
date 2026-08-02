"""
Per-bed -log10(p) scatter for one GWAS trait across all RME beds.

Each point = one RME (cell × state) bed.  Produces one panel per tool pair
(chuckle/giggle, chuckle/BITS, giggle/BITS).

Usage:
    uv run workflow/scripts/tool_scatter.py \\
        --per-bed Rheumatoid_arthritis \\
        --chuckle data/chuckle \\
        --giggle data/giggle \\
        --bits data/bits/1000 \\
        -o data/plots/scatter_bed_pvals_Rheumatoid_arthritis_1000.png --no-show

    uv run workflow/scripts/tool_scatter.py \\
        --per-bed Rheumatoid_arthritis \\
        --chuckle data/chuckle \\
        --giggle data/giggle \\
        --mc data/montecarlo \\
        -o data/plots/scatter_mc_Rheumatoid_arthritis.png --no-show
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

sys.path.insert(0, str(Path(__file__).parent))
from utils import FLOAT_TINY, read_chuckle_records, read_tsv, strip_bed_name

_TOOL_COLORS = {
    "bits": "#e07b54",
    "giggle": "#54a0e0",
    "chuckle": "#7be07b",
    "mc": "#e0c854",
}


# ──────────────────────────────────────────────────────────────
# Per-bed loaders — return dict[bed_name → -log10(p)]
# ──────────────────────────────────────────────────────────────


def load_bits_per_bed(tsv_path: str) -> dict[str, float]:
    _, rows = read_tsv(Path(tsv_path))
    result: dict[str, float] = {}
    for r in rows:
        try:
            p = max(float(r["p_value"]), FLOAT_TINY)
            result[r["name"]] = -math.log10(p)
        except (KeyError, ValueError):
            pass
    return result


def load_giggle_per_bed(tsv_path: str, col: str = "fishers_right_tail") -> dict[str, float]:
    _, rows = read_tsv(Path(tsv_path))
    result: dict[str, float] = {}
    for r in rows:
        try:
            name = strip_bed_name(Path(r["file"]).name)
            p = max(float(r[col]), FLOAT_TINY)
            result[name] = -math.log10(p)
        except (KeyError, ValueError):
            pass
    return result


def load_mc_per_bed(tsv_path: str, trials: int) -> dict[str, float]:
    """Load MC p-values for a specific trial milestone from a multi-milestone TSV."""
    # Floor at 1/trials: empirical p can't be smaller than one exceedance out of N.
    # FLOAT_TINY would produce -log10(p) ≈ 308 for p=0 hits, blowing up the axis.
    p_floor = 1.0 / trials
    result: dict[str, float] = {}
    in_section = False
    with open(tsv_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("# milestone="):
                parts = dict(tok.split("=") for tok in line[2:].split())
                in_section = int(parts.get("trials", 0)) == trials
                continue
            if not in_section or not line or line.startswith("db_name"):
                continue
            cols = line.split("\t")
            if len(cols) >= 3:
                name = strip_bed_name(Path(cols[0]).name)
                try:
                    p = max(float(cols[2]), p_floor)
                    result[name] = -math.log10(p)
                except ValueError:
                    pass
    return result


def load_chuckle_per_bed(json_path: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for r in read_chuckle_records(json_path):
        try:
            name = strip_bed_name(Path(r["db_name"]).name)
            p = max(float(r["p_value"]), FLOAT_TINY)
            out[name] = -math.log10(p)
        except (KeyError, ValueError):
            pass
    return out


# ──────────────────────────────────────────────────────────────
# Plotting
# ──────────────────────────────────────────────────────────────


def _pval_panel(
    ax: Axes,
    x_beds: dict[str, float],
    y_beds: dict[str, float],
    x_label: str,
    y_label: str,
    color: str,
    raw_pval: bool = False,
) -> None:
    common = sorted(set(x_beds) & set(y_beds))
    xs = np.array([x_beds[n] for n in common])
    ys = np.array([y_beds[n] for n in common])

    if raw_pval:
        xs = np.power(10.0, -xs)
        ys = np.power(10.0, -ys)
        x_label = x_label.replace("(-log₁₀ p)", "(p-value)")
        y_label = y_label.replace("(-log₁₀ p)", "(p-value)")
        ax.plot([0, 1], [0, 1], color="white", alpha=0.25, linewidth=0.8, linestyle="--")
    else:
        lo = min(xs.min(), ys.min())
        hi = max(xs.max(), ys.max())
        pad = (hi - lo) * 0.04
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
                color="white", alpha=0.25, linewidth=0.8, linestyle="--")

    ax.scatter(xs, ys, s=12, alpha=0.55, linewidths=0, color=color)

    if len(xs) > 2:
        r = float(np.corrcoef(xs, ys)[0, 1])
        ax.text(0.97, 0.04, f"r = {r:.2f}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=8, color="white", alpha=0.7)

    ax.set_xlabel(x_label, fontsize=9)
    ax.set_ylabel(y_label, fontsize=9)
    ax.set_title(
        f"{x_label.split()[0]} vs {y_label.split()[0]}  (n={len(common)} beds)",
        fontsize=10, pad=6,
    )


def plot_pval_scatter(
    trait: str,
    chuckle_json: str,
    giggle_tsv: str | None = None,
    giggle_col: str = "fishers_right_tail",
    bits_tsv: str | None = None,
    mc_tsv: str | None = None,
    mc_trials: int = 10000,
    output_path: str | None = None,
    show: bool = True,
    raw_pval: bool = False,
) -> None:
    """Per-bed -log10(p) scatter for one trait.  All provided tools are compared pairwise."""
    from itertools import combinations

    all_beds: dict[str, dict[str, float]] = {
        "chuckle": load_chuckle_per_bed(chuckle_json),
    }
    if giggle_tsv:
        all_beds["giggle"] = load_giggle_per_bed(giggle_tsv, col=giggle_col)
    if bits_tsv:
        all_beds["bits"] = load_bits_per_bed(bits_tsv)
    if mc_tsv:
        all_beds["mc"] = load_mc_per_bed(mc_tsv, mc_trials)

    tools = list(all_beds)
    if len(tools) < 2:
        raise ValueError("Provide at least one of giggle_tsv or bits_tsv.")

    pairs = list(combinations(tools, 2))

    plt.style.use("dark_background")
    fig, axes = plt.subplots(
        1, len(pairs), figsize=(6 * len(pairs), 5.5),
        squeeze=False, layout="constrained",
    )

    for col, (x_tool, y_tool) in enumerate(pairs):
        _pval_panel(
            axes[0, col],
            all_beds[x_tool], all_beds[y_tool],
            f"{x_tool}  (-log₁₀ p)", f"{y_tool}  (-log₁₀ p)",
            color=_TOOL_COLORS.get(y_tool, "#aaaaaa"),
            raw_pval=raw_pval,
        )

    fig.suptitle(
        f"{trait.replace('_', ' ')}  |  per-bed p-value comparison across RME",
        fontsize=12, fontweight="bold",
    )

    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=180, bbox_inches="tight")
        print(f"Saved to {output_path}")

    if show:
        plt.show()


# ──────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--per-bed", metavar="TRAIT", required=True,
                        help="Trait name (e.g. Rheumatoid_arthritis).")
    parser.add_argument("--chuckle", metavar="DIR", required=True,
                        help="Directory of chuckle query JSONs (data/chuckle/)")
    parser.add_argument("--giggle", metavar="DIR",
                        help="Directory of giggle search TSVs (data/giggle/)")
    parser.add_argument("--giggle-col", metavar="COL", default="fishers_right_tail",
                        help="giggle column to use as score (default: fishers_right_tail)")
    parser.add_argument("--bits", metavar="DIR",
                        help="Directory of BITS sweep TSVs (data/bits/<trials>/)")
    parser.add_argument("--mc", metavar="DIR",
                        help="Directory of MC multi-milestone TSVs (data/montecarlo/)")
    parser.add_argument("--mc-trials", metavar="N", type=int, default=10000,
                        help="Which trial milestone to use from the MC TSV (default: 10000)")
    parser.add_argument("-o", "--output", help="Output image path")
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument("--raw-pval", action="store_true",
                        help="Plot raw p-values instead of -log10(p)")
    args = parser.parse_args()

    trait = args.per_bed
    chuckle_json = str(Path(args.chuckle) / f"{trait}.json")
    giggle_tsv = str(Path(args.giggle) / f"{trait}.tsv") if args.giggle else None
    bits_tsv = str(Path(args.bits) / f"{trait}.tsv") if args.bits else None
    mc_tsv = str(Path(args.mc) / f"{trait}.tsv") if args.mc else None
    plot_pval_scatter(
        trait, chuckle_json,
        giggle_tsv=giggle_tsv,
        giggle_col=args.giggle_col,
        bits_tsv=bits_tsv,
        mc_tsv=mc_tsv,
        mc_trials=args.mc_trials,
        output_path=args.output,
        show=not args.no_show,
        raw_pval=args.raw_pval,
    )


if __name__ == "__main__":
    main()
