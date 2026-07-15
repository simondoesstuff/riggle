"""Plot ROC-AUC benchmark for GWAS disease-tissue enrichment tools.

Reads the JSON artifact produced by gwas_roc.py and renders one panel per
method.  Each panel shows all per-disease ROC curves (thin, semi-transparent)
overlaid with the rank-normalised combined ROC curve (bold), AUC in legend.

Usage:
    uv run workflow/scripts/gwas_roc_plot.py
        --roc-data data/gwas/roc_data.json
        --output data/plots/gwas_roc.png
        [--no-show]
"""

from __future__ import annotations

import argparse
import json
import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Method display config: key → (label, color)
METHOD_STYLE: dict[str, tuple[str, str]] = {
    "giggle": ("Giggle (combo score)", "#4C72B0"),
    "chuckle": ("Chuckle (LLR)", "#DD8452"),
}


def _plot_method(ax: plt.Axes, method: str, data: dict) -> None:
    label_base, color = METHOD_STYLE.get(method, (method, "#888888"))
    diseases = data.get("diseases", {})
    comb = data.get("combined", {})

    ax.plot([0, 1], [0, 1], color="0.7", lw=0.8, ls="--", zorder=1)

    # Individual disease curves — thin, semi-transparent
    for disease, ddata in diseases.items():
        fpr = ddata["fpr"]
        tpr = ddata["tpr"]
        auc = ddata["auc"]
        ax.plot(fpr, tpr, lw=0.7, alpha=0.35, color=color, zorder=2)

    # Combined curve — bold, with AUC in legend
    if comb:
        n = len(diseases)
        auc = comb["auc"]
        ax.plot(
            comb["fpr"], comb["tpr"],
            lw=2.5, color=color, zorder=3,
            label=f"Combined  AUC = {auc:.3f}  (n={n})",
        )

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xlabel("False Positive Rate", fontsize=9)
    ax.set_ylabel("True Positive Rate", fontsize=9)
    ax.set_title(label_base, fontsize=10, fontweight="bold")
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.25))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.25))
    ax.tick_params(labelsize=8)
    ax.legend(loc="lower right", fontsize=8)


def plot_roc(
    roc_data: dict,
    output: str | None = None,
    show: bool = True,
) -> None:
    methods = roc_data.get("methods", {})
    if not methods:
        print("No method data found.", file=sys.stderr)
        return

    n_methods = len(methods)
    fig, axes = plt.subplots(
        1, n_methods,
        figsize=(4.2 * n_methods, 4.4),
        layout="constrained",
        squeeze=False,
    )

    for ax, (method, data) in zip(axes[0], methods.items()):
        _plot_method(ax, method, data)

    if output:
        fig.savefig(output, dpi=200, bbox_inches="tight")
        print(f"Saved to {output}", file=sys.stderr)
    if show:
        plt.show()
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--roc-data", required=True, help="Path to roc_data.json")
    ap.add_argument("--output", help="Output plot path (e.g. gwas_roc.png)")
    ap.add_argument("--no-show", action="store_true")
    args = ap.parse_args()

    with open(args.roc_data) as f:
        roc_data = json.load(f)

    plot_roc(roc_data, output=args.output, show=not args.no_show)


if __name__ == "__main__":
    main()
