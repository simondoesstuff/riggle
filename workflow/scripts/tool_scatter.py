"""
Per-bed -log10(p) scatter for one GWAS trait across all RME beds.

Each point = one RME (cell × state) bed.  Produces one panel per tool pair
(chuckle/giggle, chuckle/BITS, giggle/BITS).

Usage:
    uv run scripts/tool_scatter.py \\
        --per-bed Rheumatoid_arthritis \\
        --chuckle data/chuckle \\
        --giggle data/giggle \\
        --bits data/bits/1000 \\
        --output imgs/scatter_bed_pvals_ra_1000.png
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

_TINY = float(np.finfo(np.float64).tiny)

_TOOL_COLORS = {
    "bits": "#e07b54",
    "giggle": "#54a0e0",
    "chuckle": "#7be07b",
}


# ──────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────


def _strip_name(filename: str) -> str:
    for ext in (".clean.bed.gz", ".bed.gz", ".clean.bed", ".bed"):
        if filename.endswith(ext):
            return filename[: -len(ext)]
    return filename


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Return (header, rows) from a TSV, stripping leading '#' from any line."""
    with open(path) as f:
        lines = [l.rstrip("\n\t") for l in f if l.strip()]
    if not lines:
        return [], []
    header = lines[0].lstrip("#").split("\t")
    n = len(header)
    rows = []
    for line in lines[1:]:
        parts = line.lstrip("#").split("\t")
        if len(parts) >= n:
            rows.append(dict(zip(header, parts[:n])))
    return header, rows


# ──────────────────────────────────────────────────────────────
# Per-bed loaders — return dict[bed_name → -log10(p)]
# ──────────────────────────────────────────────────────────────


def load_bits_per_bed(tsv_path: str) -> dict[str, float]:
    _, rows = _read_tsv(Path(tsv_path))
    result: dict[str, float] = {}
    for r in rows:
        try:
            p = max(float(r["p_value"]), _TINY)
            result[r["name"]] = -math.log10(p)
        except (KeyError, ValueError):
            pass
    return result


def load_giggle_per_bed(tsv_path: str, col: str = "fishers_right_tail") -> dict[str, float]:
    _, rows = _read_tsv(Path(tsv_path))
    result: dict[str, float] = {}
    for r in rows:
        try:
            name = _strip_name(Path(r["file"]).name)
            p = max(float(r[col]), _TINY)
            result[name] = -math.log10(p)
        except (KeyError, ValueError):
            pass
    return result


def load_chuckle_per_bed(json_path: str) -> dict[str, float]:
    import json as _json

    with open(json_path) as f:
        data = _json.load(f)
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
    parser.add_argument("-o", "--output", help="Output image path")
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument("--raw-pval", action="store_true",
                        help="Plot raw p-values instead of -log10(p)")
    args = parser.parse_args()

    trait = args.per_bed
    chuckle_json = str(Path(args.chuckle) / f"{trait}.json")
    giggle_tsv = str(Path(args.giggle) / f"{trait}.tsv") if args.giggle else None
    bits_tsv = str(Path(args.bits) / f"{trait}.tsv") if args.bits else None
    plot_pval_scatter(
        trait, chuckle_json,
        giggle_tsv=giggle_tsv,
        giggle_col=args.giggle_col,
        bits_tsv=bits_tsv,
        output_path=args.output,
        show=not args.no_show,
        raw_pval=args.raw_pval,
    )


if __name__ == "__main__":
    main()
