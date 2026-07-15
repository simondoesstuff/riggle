"""
Scatter plots comparing chuckle p-values against the three regioners methods
(shuffle, circle, novl).

Produces a 2×3 figure: rows = scale (-log10 / raw p-value), columns = method.

Modes:
    --trait TRAIT      one GWAS trait (one point per RME bed)
    --all-traits       all GWAS traits combined (one point per trait×bed)

Usage:
    uv run workflow/scripts/regioners_scatter.py \\
        --trait Rheumatoid_arthritis \\
        --chuckle data/chuckle \\
        --regioners data/regioners \\
        -o data/plots/scatter_regioners_Rheumatoid_arthritis.png --no-show

    uv run workflow/scripts/regioners_scatter.py \\
        --all-traits \\
        --chuckle data/chuckle \\
        --regioners data/regioners \\
        -o data/plots/scatter_regioners_all.png --no-show
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

_TINY = float(np.finfo(np.float64).tiny)
_METHOD_COLORS = {
    "shuffle": "#e07b54",
    "circle": "#54a0e0",
    "novl": "#c97be0",
}
_METHODS = ["shuffle", "circle", "novl"]


def _strip_name(filename: str) -> str:
    for ext in (".clean.bed.gz", ".bed.gz", ".clean.bed", ".bed"):
        if filename.endswith(ext):
            return filename[: -len(ext)]
    return filename


def load_chuckle(json_path: Path) -> dict[str, float]:
    """Return {bed_name: enrichment_p_value}."""
    with open(json_path) as f:
        data = json.load(f)
    results = data.get("results", data) if isinstance(data, dict) else data
    out: dict[str, float] = {}
    for r in results:
        try:
            name = _strip_name(Path(r["db_name"]).name)
            out[name] = max(float(r["p_value"]), _TINY)
        except (KeyError, ValueError):
            pass
    return out


def load_regioners(tsv_path: Path) -> dict[str, float]:
    """Return {bed_name: enrichment_p_value}, converting alt='l' to enrichment direction."""
    out: dict[str, float] = {}
    with open(tsv_path) as f:
        header = f.readline().lstrip("#").rstrip().split("\t")
        p_idx = header.index("p_val")
        alt_idx = header.index("alt")
        name_idx = header.index("name")
        for line in f:
            cols = line.rstrip().split("\t")
            if len(cols) <= max(p_idx, alt_idx, name_idx):
                continue
            try:
                raw_p = float(cols[p_idx])
                alt = cols[alt_idx]
                # Convert to enrichment p-value: "g" is already enrichment;
                # "l" tests depletion so enrichment p ≈ 1 - p.
                if alt == "g":
                    p = raw_p
                elif alt == "l":
                    p = 1.0 - raw_p
                else:
                    p = raw_p  # "t" or unknown — use as-is
                out[cols[name_idx]] = max(p, _TINY)
            except (ValueError, IndexError):
                pass
    return out


def _pairs(chuckle: dict[str, float], other: dict[str, float]) -> tuple[np.ndarray, np.ndarray]:
    """Return aligned (chuckle_p, other_p) arrays over the common bed names."""
    common = sorted(set(chuckle) & set(other))
    return (
        np.array([chuckle[n] for n in common]),
        np.array([other[n] for n in common]),
    )


def _collect_all_pairs(
    chuckle_dir: Path, regioners_dir: Path, method: str
) -> tuple[np.ndarray, np.ndarray]:
    """Concatenate (chuckle_p, regioners_p) pairs across all available traits."""
    cx_all, ry_all = [], []
    for json_path in sorted(chuckle_dir.glob("*.json")):
        trait = json_path.stem
        tsv_path = regioners_dir / method / f"{trait}.tsv"
        if not tsv_path.exists():
            continue
        cx, ry = _pairs(load_chuckle(json_path), load_regioners(tsv_path))
        cx_all.append(cx)
        ry_all.append(ry)
    if not cx_all:
        raise ValueError(f"No matching trait files found for method={method!r}")
    return np.concatenate(cx_all), np.concatenate(ry_all)


def _scatter_panel(
    ax: Axes,
    cx: np.ndarray,
    ry: np.ndarray,
    method: str,
    log10: bool,
    color: str,
) -> None:
    if log10:
        xs = -np.log10(cx)
        ys = -np.log10(ry)
        ax.set_xlabel("chuckle  (−log₁₀ p)", fontsize=9)
        ax.set_ylabel(f"{method}  (−log₁₀ p)", fontsize=9)
        lo = min(xs.min(), ys.min())
        hi = max(xs.max(), ys.max())
        pad = (hi - lo) * 0.04
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
                color="white", alpha=0.25, linewidth=0.8, linestyle="--")
    else:
        xs, ys = cx, ry
        ax.set_xlabel("chuckle  (p-value)", fontsize=9)
        ax.set_ylabel(f"{method}  (p-value)", fontsize=9)
        ax.plot([0, 1], [0, 1], color="white", alpha=0.25, linewidth=0.8, linestyle="--")

    ax.scatter(xs, ys, s=4, alpha=0.2, linewidths=0, color=color)

    if len(xs) > 2:
        r = float(np.corrcoef(xs, ys)[0, 1])
        ax.text(0.97, 0.04, f"r = {r:.2f}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=8, color="white", alpha=0.7)

    ax.set_title(f"{method}  (n={len(xs):,})", fontsize=10, pad=5)


def _make_figure(
    pairs_per_method: dict[str, tuple[np.ndarray, np.ndarray]],
    title: str,
    output: Path | None,
    show: bool,
) -> None:
    plt.style.use("dark_background")
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), layout="constrained")

    for col, method in enumerate(_METHODS):
        cx, ry = pairs_per_method[method]
        color = _METHOD_COLORS[method]
        _scatter_panel(axes[0, col], cx, ry, method, log10=True,  color=color)
        _scatter_panel(axes[1, col], cx, ry, method, log10=False, color=color)

    fig.suptitle(title, fontsize=13, fontweight="bold")

    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=180, bbox_inches="tight")
        print(f"Saved to {output}")

    if show:
        plt.show()
    else:
        plt.close(fig)


def plot(
    trait: str,
    chuckle_dir: Path,
    regioners_dir: Path,
    output: Path | None,
    show: bool,
) -> None:
    chuckle = load_chuckle(chuckle_dir / f"{trait}.json")
    pairs = {m: _pairs(chuckle, load_regioners(regioners_dir / m / f"{trait}.tsv"))
             for m in _METHODS}
    _make_figure(pairs, f"{trait.replace('_', ' ')}  |  chuckle vs regioners", output, show)


def plot_all_traits(
    chuckle_dir: Path,
    regioners_dir: Path,
    output: Path | None,
    show: bool,
) -> None:
    pairs = {m: _collect_all_pairs(chuckle_dir, regioners_dir, m) for m in _METHODS}
    _make_figure(pairs, "All traits  |  chuckle vs regioners", output, show)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--trait", metavar="TRAIT", help="Single GWAS trait name")
    mode.add_argument("--all-traits", action="store_true",
                      help="Combine all traits in --chuckle / --regioners dirs")
    parser.add_argument("--chuckle", required=True, metavar="DIR",
                        help="Directory of chuckle JSON files")
    parser.add_argument("--regioners", required=True, metavar="DIR",
                        help="Base regioners directory (method TSVs in <DIR>/<method>/<trait>.tsv)")
    parser.add_argument("-o", "--output", metavar="PATH")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()

    chuckle_dir = Path(args.chuckle)
    regioners_dir = Path(args.regioners)
    output = Path(args.output) if args.output else None
    show = not args.no_show

    if args.all_traits:
        plot_all_traits(chuckle_dir, regioners_dir, output, show)
    else:
        plot(args.trait, chuckle_dir, regioners_dir, output, show)


if __name__ == "__main__":
    main()
