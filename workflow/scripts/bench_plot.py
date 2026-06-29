"""Plot asymptotic complexity from at-scale benchmark trial JSONs.

Each trial JSON (produced by bench_trial.py via Snakemake) contains structured
per-tool timing and size data.  This script aggregates all trials and fits
asymptotic complexity models to each tool's curves.

Usage:
    uv run workflow/scripts/bench_plot.py data/bench/*.json \\
        --index data/plots/bench_index.png \\
        --query data/plots/bench_query.png \\
        --size  data/plots/bench_size.png
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# One free parameter (a): f(x, a) = a * g(x)
_LINEAR = lambda x, a: a * x  # noqa: E731
_NLOGN = lambda x, a: a * x * np.log2(x)  # noqa: E731

_TOOL_FIT = {"Chuckle": _LINEAR, "IGD": _LINEAR, "Giggle": _NLOGN}


def _fit_and_plot(ax: plt.Axes, xs: list, ys: list, label: str, color: str) -> None:
    xs_arr = np.array(xs, dtype=float)
    ys_arr = np.array(ys, dtype=float)
    valid = np.isfinite(xs_arr) & np.isfinite(ys_arr) & (xs_arr > 0) & (ys_arr > 0)
    xs_arr, ys_arr = xs_arr[valid], ys_arr[valid]
    if len(xs_arr) == 0:
        return
    ax.scatter(xs_arr, ys_arr, color=color, label=label, zorder=5, s=40)
    fit_fn = _TOOL_FIT.get(label, _LINEAR)
    if len(xs_arr) >= 2:
        try:
            (a,), _ = curve_fit(fit_fn, xs_arr, ys_arr, p0=[1.0], maxfev=10_000)
            xf = np.logspace(np.log10(xs_arr.min()), np.log10(xs_arr.max()), 300)
            ax.plot(xf, fit_fn(xf, a), color=color, linestyle="--", linewidth=1.2)
        except (RuntimeError, ValueError):
            pass


def _save_plot(
    title: str,
    xlabel: str,
    ylabel: str,
    tool_data: dict[str, tuple[list, list]],
    path: str,
) -> None:
    colors = [p["color"] for p in plt.rcParams["axes.prop_cycle"]]
    fig, ax = plt.subplots(figsize=(9, 5))
    for i, (tool, (xs, ys)) in enumerate(tool_data.items()):
        _fit_and_plot(ax, xs, ys, tool, colors[i % len(colors)])
    ax.set_xscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=8)
    ax.grid(True, which="both", linestyle=":", alpha=0.5)
    fig.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("trials", nargs="+", help="Benchmark trial JSON files")
    parser.add_argument("--index", required=True, help="Output path for index-time plot")
    parser.add_argument("--query", required=True, help="Output path for query-time plot")
    parser.add_argument("--size", required=True, help="Output path for index-size plot")
    args = parser.parse_args()

    trials = []
    for path in args.trials:
        with open(path) as f:
            trials.append(json.load(f))
    trials.sort(key=lambda t: t["total_intervals"])

    index_data: dict[str, tuple[list, list]] = {}
    query_data: dict[str, tuple[list, list]] = {}
    size_data: dict[str, tuple[list, list]] = {}

    for t in trials:
        n = t["total_intervals"]
        for tool, res in t.get("results", {}).items():
            if res.get("index_s") is not None:
                xs, ys = index_data.setdefault(tool, ([], []))
                xs.append(n)
                ys.append(res["index_s"])
            if res.get("query_s") is not None:
                xs, ys = query_data.setdefault(tool, ([], []))
                xs.append(n)
                ys.append(res["query_s"])
            if res.get("index_bytes") is not None:
                xs, ys = size_data.setdefault(tool, ([], []))
                xs.append(n)
                ys.append(res["index_bytes"] / 1024**3)

    _save_plot(
        "Index Time vs Total Intervals",
        "total source intervals (n)",
        "wall-clock time (s)",
        index_data,
        args.index,
    )
    _save_plot(
        "Query Time vs Total Intervals",
        "total source intervals (n)",
        "wall-clock time (s)",
        query_data,
        args.query,
    )
    _save_plot(
        "Index Size vs Total Intervals",
        "total source intervals (n)",
        "index size (GB)",
        size_data,
        args.size,
    )


if __name__ == "__main__":
    main()
