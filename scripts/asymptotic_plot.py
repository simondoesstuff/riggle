"""
Asymptotic complexity plotter.

Usage:
    from asymptotic_plot import plot_asymptotic, Curve

    nlogn = lambda x, a: a * x * np.log2(x)

    plot_asymptotic([
        Curve(label="algo A", xs=[100,1000,10000], ys=[0.1, 1.5, 18.0], fit=nlogn),
        Curve(label="algo B", xs=[100,1000,10000], ys=[0.05, 0.6, 7.2], fit=lambda x, a: a * x),
    ])
"""

from dataclasses import dataclass
from typing import Callable, Optional
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


@dataclass
class Curve:
    label: str
    xs: list[float]
    ys: list[float]
    # f(x, a) -> y, where a is the single scale parameter to fit.
    # NaN/None/inf values in xs or ys are silently dropped; curves with
    # fewer than 2 valid points are plotted without a fit line.
    fit: Optional[Callable] = None
    color: Optional[str] = None


def plot_asymptotic(
    curves: list[Curve],
    title: str = "Asymptotic Complexity",
    xlabel: str = "n",
    ylabel: str = "time",
    show: bool = True,
    save_path: Optional[str] = None,
):
    fig, ax = plt.subplots(figsize=(9, 5))
    prop_cycle = plt.rcParams["axes.prop_cycle"]
    colors = [p["color"] for p in prop_cycle]

    # Pre-process all curves so we can determine the global x range for extrapolation.
    processed = []
    for curve in curves:
        pairs = list(zip(curve.xs, curve.ys))
        if not pairs:
            continue
        xs_raw, ys_raw = zip(*pairs)
        xs = np.array(xs_raw, dtype=float)
        ys = np.array(ys_raw, dtype=float)
        valid = np.isfinite(xs) & np.isfinite(ys) & (xs > 0)
        xs, ys = xs[valid], ys[valid]
        if len(xs) == 0:
            continue
        processed.append((curve, xs, ys))

    if not processed:
        return fig, ax

    global_xmin = min(xs.min() for _, xs, _ in processed)
    global_xmax = max(xs.max() for _, xs, _ in processed)
    x_fit_global = np.logspace(np.log10(global_xmin), np.log10(global_xmax), 300)

    for i, (curve, xs, ys) in enumerate(processed):
        color = curve.color or colors[i % len(colors)]

        if curve.fit is not None and len(xs) >= 2:
            try:
                (a,), _ = curve_fit(curve.fit, xs, ys, p0=[1.0], maxfev=10_000)
                y_fit = curve.fit(x_fit_global, a)
                ax.plot(x_fit_global, y_fit, color=color, linestyle="--", linewidth=1.2)
            except (RuntimeError, ValueError):
                pass

        ax.scatter(xs, ys, color=color, label=curve.label, zorder=5, s=40)

    ax.set_xscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=8)
    ax.grid(True, which="both", linestyle=":", alpha=0.5)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150)
    if show:
        plt.show()
    return fig, ax


if __name__ == "__main__":
    rng = np.random.default_rng(42)
    ns = np.array([100, 300, 1000, 3000, 10_000, 30_000])

    plot_asymptotic([
        Curve(
            label="algo A",
            xs=ns,
            ys=2e-7 * ns * np.log2(ns) + rng.normal(0, 0.002, len(ns)),
            fit=lambda x, a: a * x * np.log2(x),
        ),
        Curve(
            label="algo B",
            xs=ns,
            ys=5e-7 * ns + rng.normal(0, 0.001, len(ns)),
            fit=lambda x, a: a * x,
        ),
    ], title="Smoke test")
