"""
Parse at_scale bench output and plot index or query timing.

Usage:
    python scripts/bench_plot.py scripts/at_scale_many.txt --phase index
    python scripts/bench_plot.py scripts/at_scale_many.txt --phase query --save plot.png
"""

import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from asymptotic_plot import plot_asymptotic, Curve

# ── Parsing ───────────────────────────────────────────────────────────────────

_TIME_RE = re.compile(r"real\s+(\d+)m([\d.]+)s")
_RUN_RE = re.compile(r"RUN: src_files=(\d+), src_records=(\d+)")
_SECTION_RE = re.compile(r"={10,}\s+(\w+)\s+={10,}")


def _parse_time(line: str) -> float | None:
    m = _TIME_RE.search(line)
    return int(m.group(1)) * 60 + float(m.group(2)) if m else None


def parse_bench(text: str) -> dict[str, dict[str, list[tuple[int, float]]]]:
    """
    Parse benchmark output.

    Returns:
        { algo_name: { 'index': [(n, seconds), ...], 'query': [(n, seconds), ...] } }

    where n = total source intervals for that run.  Algos absent from a run or
    that failed/timed out are omitted from that run's data.
    """
    results: dict = {}
    current_n: int | None = None
    current_algo: str | None = None
    current_phase: str | None = None   # 'index' | 'query'
    phase_failed = False
    pending_time: float | None = None

    def commit() -> None:
        nonlocal pending_time, phase_failed
        if (
            pending_time is not None
            and not phase_failed
            and current_algo
            and current_phase
            and current_n is not None
        ):
            results.setdefault(current_algo, {"index": [], "query": []})[
                current_phase
            ].append((current_n, pending_time))
        pending_time = None
        phase_failed = False

    for line in text.splitlines():
        # New run block
        m = _RUN_RE.search(line)
        if m:
            commit()
            current_n = int(m.group(1)) * int(m.group(2))
            current_algo = None
            current_phase = None
            continue

        # Tool section header
        m = _SECTION_RE.search(line)
        if m:
            commit()
            algo = m.group(1).strip()
            if algo == "Native":
                algo = "Riggle"
            current_algo = algo
            current_phase = None
            continue

        # Phase markers
        if ">>>   Building index   <<<" in line:
            commit()
            current_phase = "index"
            continue

        if ">>>   Querying   <<<" in line:
            commit()
            current_phase = "query"
            continue

        # Failure / timeout marker — discard any pending time for this phase
        if "[!]" in line:
            phase_failed = True
            pending_time = None
            continue

        t = _parse_time(line)
        if t is not None:
            pending_time = t  # keep the last 'real' seen for this phase

    commit()
    return results


# ── Plotting ──────────────────────────────────────────────────────────────────

_LINEAR = lambda x, a: a * x              # noqa: E731
_NLOGN  = lambda x, a: a * x * np.log2(x)  # noqa: E731

# Riggle and IGD are linear; Giggle is NlogN.
_ALGO_FIT: dict[str, object] = {
    "Riggle": _LINEAR,
    "IGD":    _LINEAR,
    "Giggle": _NLOGN,
}


def bench_curves(
    parsed: dict[str, dict[str, list[tuple[int, float]]]],
    phase: str = "index",
) -> list[Curve]:
    """
    Convert parsed bench data into Curve objects ready for plot_asymptotic.

    Args:
        parsed:  output of parse_bench()
        phase:   'index' or 'query'
    """
    curves = []
    for algo, data in parsed.items():
        pts = data.get(phase, [])
        if not pts:
            continue
        xs, ys = zip(*pts)
        curves.append(Curve(
            label=algo,
            xs=list(xs),
            ys=list(ys),
            fit=_ALGO_FIT.get(algo, _LINEAR),
        ))
    return curves


# ── CLI ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Benchmark output file (e.g. at_scale_many.txt)")
    parser.add_argument(
        "--phase",
        choices=["index", "query"],
        default="index",
        help="Which timing to plot (default: index)",
    )
    parser.add_argument("--save", metavar="PATH", help="Save plot to file instead of showing")
    args = parser.parse_args()

    text = Path(args.input).read_text()
    parsed = parse_bench(text)
    curves = bench_curves(parsed, phase=args.phase)  # fits per-algo

    phase_label = args.phase.capitalize()
    plot_asymptotic(
        curves,
        title=f"{phase_label} Time vs Total Source Intervals",
        xlabel="total source intervals (n)",
        ylabel="wall-clock time (s)",
        show=args.save is None,
        save_path=args.save,
    )
