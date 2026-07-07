#!/usr/bin/env python3

"""
bed_stats.py

Report statistics about a series of BED files, including histograms
for interval sizes and file sizes (intervals per file).
"""

import argparse
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import gzip
from typing import Generator

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from numpy.typing import NDArray
from tqdm import tqdm


def load_bed(path: Path) -> Generator[tuple[str, int, int], None, None]:
    open_fn = gzip.open if str(path).endswith(".gz") else open
    with open_fn(path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            yield parts[0], int(parts[1]), int(parts[2])


@dataclass
class BedFileStats:
    """Statistics for a single BED file."""

    path: Path
    n_intervals: int
    interval_sizes: list[int]
    total_bp: int
    chromosomes: Counter[str]


@dataclass
class AggregateStats:
    """Aggregate statistics across all BED files."""

    n_files: int
    total_intervals: int
    total_bp: int
    file_sizes: NDArray[np.int64]  # intervals per file
    all_interval_sizes: NDArray[np.int64]
    chromosome_counts: Counter[str]
    per_file: list[BedFileStats]


def compute_file_stats(path: Path) -> BedFileStats:
    """Compute statistics for a single BED file."""
    intervals = list(load_bed(path))
    sizes = [end - start for _, start, end in intervals]
    chroms: Counter[str] = Counter(chrom for chrom, _, _ in intervals)

    return BedFileStats(
        path=path,
        n_intervals=len(intervals),
        interval_sizes=sizes,
        total_bp=sum(sizes),
        chromosomes=chroms,
    )


def compute_aggregate_stats(files: list[Path]) -> AggregateStats:
    """Compute aggregate statistics across all BED files."""
    per_file = [compute_file_stats(f) for f in tqdm(files, desc="Processing BED files")]

    all_sizes: list[int] = []
    all_chroms: Counter[str] = Counter()
    file_sizes = []

    for stats in per_file:
        all_sizes.extend(stats.interval_sizes)
        all_chroms.update(stats.chromosomes)
        file_sizes.append(stats.n_intervals)

    return AggregateStats(
        n_files=len(files),
        total_intervals=sum(file_sizes),
        total_bp=sum(s.total_bp for s in per_file),
        file_sizes=np.array(file_sizes, dtype=np.int64),
        all_interval_sizes=np.array(all_sizes, dtype=np.int64),
        chromosome_counts=all_chroms,
        per_file=per_file,
    )


def print_summary(stats: AggregateStats) -> None:
    """Print summary statistics to stdout."""
    print("=" * 60)
    print("BED FILE STATISTICS SUMMARY")
    print("=" * 60)
    print()

    print(f"Files analyzed:      {stats.n_files:,}")
    print(f"Total intervals:     {stats.total_intervals:,}")
    print(f"Total base pairs:    {stats.total_bp:,}")
    print()

    print("INTERVALS PER FILE:")
    print(f"  Mean:              {np.mean(stats.file_sizes):,.1f}")
    print(f"  Median:            {np.median(stats.file_sizes):,.1f}")
    print(f"  Std:               {np.std(stats.file_sizes):,.1f}")
    print(f"  Min:               {np.min(stats.file_sizes):,}")
    print(f"  Max:               {np.max(stats.file_sizes):,}")
    print()

    print("INTERVAL SIZES (bp):")
    print(f"  Mean:              {np.mean(stats.all_interval_sizes):,.1f}")
    print(f"  Median:            {np.median(stats.all_interval_sizes):,.1f}")
    print(f"  Std:               {np.std(stats.all_interval_sizes):,.1f}")
    print(f"  Min:               {np.min(stats.all_interval_sizes):,}")
    print(f"  Max:               {np.max(stats.all_interval_sizes):,}")
    print()

    # Percentiles
    percentiles = [25, 50, 75, 90, 95, 99]
    print("INTERVAL SIZE PERCENTILES:")
    for p in percentiles:
        val = np.percentile(stats.all_interval_sizes, p)
        print(f"  {p}th:             {val:,.0f} bp")
    print()

    # Chromosome distribution
    print("CHROMOSOME DISTRIBUTION (top 10):")
    for chrom, count in stats.chromosome_counts.most_common(10):
        pct = 100 * count / stats.total_intervals
        print(f"  {chrom:8s}         {count:>10,} ({pct:5.1f}%)")
    print()


def plot_histograms(
    stats: AggregateStats,
    output_path: Path | None = None,
    log_scale: bool = False,
) -> None:
    """Generate histogram plots for interval sizes and file sizes."""
    _fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    sns.set_style("whitegrid")

    # 1. Interval size histogram
    ax1 = axes[0, 0]
    sns.histplot(stats.all_interval_sizes, bins=50, ax=ax1, color="steelblue")
    ax1.set_xlabel("Interval Size (bp)")
    ax1.set_ylabel("Count")
    ax1.set_title("Interval Size Distribution")
    if log_scale:
        ax1.set_yscale("log")

    # 2. Interval size histogram (log x-scale for wide distributions)
    ax2 = axes[0, 1]
    log_sizes = np.log10(stats.all_interval_sizes + 1)
    sns.histplot(log_sizes, bins=50, ax=ax2, color="darkorange")
    ax2.set_xlabel("log10(Interval Size + 1)")
    ax2.set_ylabel("Count")
    ax2.set_title("Interval Size Distribution (log scale)")
    if log_scale:
        ax2.set_yscale("log")

    # 3. File size (intervals per file) histogram
    ax3 = axes[1, 0]
    sns.histplot(
        stats.file_sizes, bins=min(50, stats.n_files), ax=ax3, color="seagreen"
    )
    ax3.set_xlabel("Intervals per File")
    ax3.set_ylabel("Count")
    ax3.set_title("BED File Size Distribution")
    if log_scale:
        ax3.set_yscale("log")

    # 4. Chromosome distribution bar chart
    ax4 = axes[1, 1]
    top_chroms = stats.chromosome_counts.most_common(25)
    chroms, counts = zip(*top_chroms) if top_chroms else ([], [])
    ax4.bar(chroms, counts, color="mediumpurple")
    ax4.set_xlabel("Chromosome")
    ax4.set_ylabel("Interval Count")
    ax4.set_title("Chromosome Distribution")
    ax4.tick_params(axis="x", rotation=45)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"Saved plot to: {output_path}")
    else:
        plt.show()


def collect_bed_files(paths: list[Path]) -> list[Path]:
    """Collect BED files from paths (files or directories)."""
    bed_files: list[Path] = []

    for path in paths:
        if path.is_file():
            bed_files.append(path)
        elif path.is_dir():
            bed_files.extend(path.glob("*.bed"))
            bed_files.extend(path.glob("*.bed.gz"))
        else:
            print(f"Warning: Path does not exist: {path}", file=sys.stderr)

    return sorted(set(bed_files))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Report statistics about BED files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  %(prog)s data/beds/
  %(prog)s *.bed.gz --output stats.png
  %(prog)s data/beds/ --log-scale --output bed_stats.pdf
""",
    )

    parser.add_argument(
        "paths",
        type=Path,
        nargs="+",
        help="BED files or directories containing BED files",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output path for histogram plot (e.g., stats.png). Shows interactive plot if not specified.",
    )
    parser.add_argument(
        "--log-scale",
        action="store_true",
        help="Use log scale for y-axis in histograms",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip histogram generation, print stats only",
    )
    parser.add_argument(
        "--per-file",
        action="store_true",
        help="Print per-file statistics",
    )

    args = parser.parse_args()

    bed_files = collect_bed_files(args.paths)

    if not bed_files:
        print("Error: No BED files found", file=sys.stderr)
        sys.exit(1)

    print(f"Analyzing {len(bed_files)} BED files...")
    stats = compute_aggregate_stats(bed_files)

    print_summary(stats)

    if args.per_file:
        print("=" * 60)
        print("PER-FILE STATISTICS")
        print("=" * 60)
        print(f"{'File':<50} {'Intervals':>12} {'Total BP':>15}")
        print("-" * 80)
        for fs in sorted(stats.per_file, key=lambda x: x.n_intervals, reverse=True):
            name = fs.path.name[:48]
            print(f"{name:<50} {fs.n_intervals:>12,} {fs.total_bp:>15,}")
        print()

    if not args.no_plot:
        print("Preparing plots...")
        plot_histograms(stats, args.output, args.log_scale)


if __name__ == "__main__":
    main()
