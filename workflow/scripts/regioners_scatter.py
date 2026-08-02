"""
Scatter plots comparing chuckle p-values against the three regioners methods
(shuffle, circle, novl).

Produces a 3×3 figure: rows = scale (rank / -log10 / raw p-value), columns = method.

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
from scipy.stats import rankdata

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


def load_chuckle(json_path: Path) -> tuple[dict[str, float], dict[str, float]]:
    """Return ({bed_name: p_value}, {bed_name: llr})."""
    with open(json_path) as f:
        data = json.load(f)
    results = data.get("results", data) if isinstance(data, dict) else data
    p_out: dict[str, float] = {}
    llr_out: dict[str, float] = {}
    for r in results:
        try:
            name = _strip_name(Path(r["db_name"]).name)
            p_out[name] = max(float(r["p_value"]), _TINY)
            llr_out[name] = float(r["llr"])
        except (KeyError, ValueError):
            pass
    return p_out, llr_out


def load_giggle(tsv_path: Path) -> tuple[dict[str, float], dict[str, float]]:
    """Return ({bed_name: fishers_right_tail}, {bed_name: combo_score})."""
    p_out: dict[str, float] = {}
    score_out: dict[str, float] = {}
    with open(tsv_path) as f:
        header = f.readline().lstrip("#").rstrip().split("\t")
        rt_idx = header.index("fishers_right_tail")
        cs_idx = header.index("combo_score")
        file_idx = header.index("file")
        for line in f:
            cols = line.rstrip().split("\t")
            if len(cols) <= max(rt_idx, cs_idx, file_idx):
                continue
            try:
                name = _strip_name(Path(cols[file_idx]).name)
                p_out[name] = max(float(cols[rt_idx]), _TINY)
                score_out[name] = float(cols[cs_idx])
            except (ValueError, IndexError):
                pass
    return p_out, score_out


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


def _pairs(
    tool_p: dict[str, float],
    tool_score: dict[str, float],
    other: dict[str, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return aligned (tool_p, regioners_p, score_rank, p_rank) — ranks within this set."""
    common = sorted(set(tool_p) & set(other))
    cx = np.array([tool_p[n] for n in common])
    ry = np.array([other[n] for n in common])
    sc = np.array([tool_score[n] for n in common])
    score_rank = rankdata(-sc)  # rank 1 = highest score = most enriched
    p_rank     = rankdata(ry)   # rank 1 = lowest p     = most enriched
    return cx, ry, score_rank, p_rank


def _collect_all_pairs(
    tool_dir: Path,
    regioners_dir: Path,
    method: str,
    load_tool,           # callable: path -> (p_dict, score_dict)
    tool_glob: str,
    trait_from_path,     # callable: Path -> str (trait name)
    regioners_path,      # callable: (regioners_dir, method, trait) -> Path
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Concatenate per-trait (tool_p, regioners_p, score_rank, p_rank) across all traits."""
    cx_all, ry_all, sr_all, pr_all = [], [], [], []
    for tool_path in sorted(tool_dir.glob(tool_glob)):
        trait = trait_from_path(tool_path)
        tsv_path = regioners_path(regioners_dir, method, trait)
        if not tsv_path.exists():
            continue
        tool_p, tool_score = load_tool(tool_path)
        cx, ry, sr, pr = _pairs(tool_p, tool_score, load_regioners(tsv_path))
        cx_all.append(cx)
        ry_all.append(ry)
        sr_all.append(sr)
        pr_all.append(pr)
    if not cx_all:
        raise ValueError(f"No matching trait files found for method={method!r}")
    return (
        np.concatenate(cx_all),
        np.concatenate(ry_all),
        np.concatenate(sr_all),
        np.concatenate(pr_all),
    )


def _scatter_panel_rank(
    ax: Axes,
    score_rank: np.ndarray,
    p_rank: np.ndarray,
    method: str,
    tool_label: str,
    color: str,
) -> None:
    ax.scatter(score_rank, p_rank, s=4, alpha=0.2, linewidths=0, color=color)

    if len(score_rank) > 2:
        r = float(np.corrcoef(score_rank, p_rank)[0, 1])
        ax.text(0.97, 0.04, f"ρ = {r:.2f}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=8, color="white", alpha=0.7)

    ax.set_xlabel(f"{tool_label}  (rank by score)", fontsize=9)
    ax.set_ylabel(f"{method}  (rank by p)", fontsize=9)
    ax.set_title(f"{method}  (n={len(score_rank):,})", fontsize=10, pad=5)


def _scatter_panel(
    ax: Axes,
    cx: np.ndarray,
    ry: np.ndarray,
    method: str,
    log10: bool,
    tool_label: str,
    color: str,
) -> None:
    if log10:
        xs = -np.log10(cx)
        ys = -np.log10(ry)
        ax.set_xlabel(f"{tool_label}  (−log₁₀ p)", fontsize=9)
        ax.set_ylabel(f"{method}  (−log₁₀ p)", fontsize=9)
        lo = min(xs.min(), ys.min())
        hi = max(xs.max(), ys.max())
        pad = (hi - lo) * 0.04
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
                color="white", alpha=0.25, linewidth=0.8, linestyle="--")
    else:
        xs, ys = cx, ry
        ax.set_xlabel(f"{tool_label}  (p-value)", fontsize=9)
        ax.set_ylabel(f"{method}  (p-value)", fontsize=9)
        ax.plot([0, 1], [0, 1], color="white", alpha=0.25, linewidth=0.8, linestyle="--")

    ax.scatter(xs, ys, s=4, alpha=0.2, linewidths=0, color=color)

    if len(xs) > 2:
        r = float(np.corrcoef(xs, ys)[0, 1])
        ax.text(0.97, 0.04, f"r = {r:.2f}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=8, color="white", alpha=0.7)

    ax.set_title(f"{method}  (n={len(xs):,})", fontsize=10, pad=5)


def _make_figure(
    pairs_per_method: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    title: str,
    tool_label: str,
    output: Path | None,
    show: bool,
) -> None:
    plt.style.use("dark_background")
    fig, axes = plt.subplots(3, 3, figsize=(15, 13), layout="constrained")

    for col, method in enumerate(_METHODS):
        cx, ry, score_rank, p_rank = pairs_per_method[method]
        color = _METHOD_COLORS[method]
        _scatter_panel_rank(axes[0, col], score_rank, p_rank, method, tool_label, color=color)
        _scatter_panel(axes[1, col], cx, ry, method, log10=True,  tool_label=tool_label, color=color)
        _scatter_panel(axes[2, col], cx, ry, method, log10=False, tool_label=tool_label, color=color)

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
    tool_dir: Path,
    regioners_dir: Path,
    load_tool,
    tool_label: str,
    trait_file: str,
    output: Path | None,
    show: bool,
) -> None:
    tool_p, tool_score = load_tool(tool_dir / trait_file)
    pairs = {
        m: _pairs(tool_p, tool_score, load_regioners(regioners_dir / m / f"{trait}.tsv"))
        for m in _METHODS
    }
    _make_figure(
        pairs,
        f"{trait.replace('_', ' ')}  |  {tool_label} vs regioners",
        tool_label,
        output,
        show,
    )


def plot_all_traits(
    tool_dir: Path,
    regioners_dir: Path,
    load_tool,
    tool_label: str,
    tool_glob: str,
    trait_from_path,
    regioners_path,
    output: Path | None,
    show: bool,
) -> None:
    pairs = {
        m: _collect_all_pairs(
            tool_dir, regioners_dir, m, load_tool, tool_glob, trait_from_path, regioners_path
        )
        for m in _METHODS
    }
    _make_figure(pairs, f"All traits  |  {tool_label} vs regioners", tool_label, output, show)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--trait", metavar="TRAIT", help="Single GWAS trait name")
    mode.add_argument("--all-traits", action="store_true",
                      help="Combine all traits in --chuckle / --regioners dirs")
    tool_grp = parser.add_mutually_exclusive_group(required=True)
    tool_grp.add_argument("--chuckle", metavar="DIR", help="Directory of chuckle JSON files")
    tool_grp.add_argument("--giggle", metavar="DIR", help="Directory of giggle TSV files")
    parser.add_argument("--regioners", required=True, metavar="DIR",
                        help="Base regioners directory (method TSVs in <DIR>/<method>/<trait>.tsv)")
    parser.add_argument("-o", "--output", metavar="PATH")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()

    regioners_dir = Path(args.regioners)
    output = Path(args.output) if args.output else None
    show = not args.no_show
    reg_path = lambda d, m, t: d / m / f"{t}.tsv"

    if args.chuckle:
        tool_dir = Path(args.chuckle)
        load_tool = load_chuckle
        tool_label = "chuckle"
        tool_glob, trait_from_path = "*.json", lambda p: p.stem
        trait_file_fn = lambda t: f"{t}.json"
    else:
        tool_dir = Path(args.giggle)
        load_tool = load_giggle
        tool_label = "giggle"
        tool_glob, trait_from_path = "*.tsv", lambda p: p.stem
        trait_file_fn = lambda t: f"{t}.tsv"

    if args.all_traits:
        plot_all_traits(
            tool_dir, regioners_dir, load_tool, tool_label,
            tool_glob, trait_from_path, reg_path, output, show,
        )
    else:
        plot(
            args.trait, tool_dir, regioners_dir, load_tool, tool_label,
            trait_file_fn(args.trait), output, show,
        )


if __name__ == "__main__":
    main()
