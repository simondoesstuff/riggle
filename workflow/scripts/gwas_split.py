"""Download and split GWAS XLS into per-trait hg19 BED files.

Source: Supplementary Table 1 from Maurano et al. (2012) / Kundaje et al. (2015),
Nature 518, doi:10.1038/nature13835.

Columns (0-indexed): 0=Disease, 3=Chr, 4=Position_hg19.
Output: one BED file per unique disease name in the output directory.

Usage:
    uv run workflow/scripts/gwas_split.py data/gwas/raw/nature13835-s1.xls data/gwas/hg19
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import xlrd


def _safe_name(s: str) -> str:
    """Convert a disease name to a filesystem-safe filename stem."""
    s = s.strip()
    s = s.replace("'", "").replace("/", "_").replace(" ", "_")
    s = re.sub(r"[^\w-]", "", s)
    return s


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("xls", help="Downloaded GWAS XLS file")
    parser.add_argument("output_dir", help="Directory to write per-disease hg19 BED files")
    args = parser.parse_args()

    wb = xlrd.open_workbook(args.xls)
    ws = wb.sheet_by_index(0)

    # Columns are 0-indexed: Disease=0, Chr=3, Pos=4
    by_disease: dict[str, list[tuple[str, int]]] = defaultdict(list)

    for row in range(1, ws.nrows):  # skip header
        disease = str(ws.cell_value(row, 0)).strip()
        chrom = str(ws.cell_value(row, 3)).strip()
        pos_raw = ws.cell_value(row, 4)

        if not disease or not chrom or pos_raw == "":
            continue

        # Chromosome: normalise to chrN format
        chrom_str = str(chrom).strip()
        # xlrd may read integers as floats (e.g. "1.0"), strip .0
        if chrom_str.endswith(".0"):
            chrom_str = chrom_str[:-2]
        if not chrom_str.startswith("chr"):
            chrom_str = f"chr{chrom_str}"

        # Position: xlrd reads numbers as float
        try:
            pos = int(float(str(pos_raw).rstrip(".0") or pos_raw))
        except (ValueError, TypeError):
            continue

        by_disease[_safe_name(disease)].append((chrom_str, pos))

    if not by_disease:
        sys.exit("No data rows found — check column indices.")

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    for name, intervals in sorted(by_disease.items()):
        with open(out / f"{name}.bed", "w") as f:
            for chrom, pos in sorted(intervals):
                f.write(f"{chrom}\t{pos}\t{pos + 1}\n")

    print(f"Wrote {len(by_disease)} BED files to {out}")


if __name__ == "__main__":
    main()
