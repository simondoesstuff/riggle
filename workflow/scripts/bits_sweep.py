"""
Run bits_test for one GWAS query against all RME beds, producing a TSV.

Usage:
    uv run scripts/bits_sweep.py --query QUERY.bed --beds-dir DIR
                                  --universe UNIVERSE.bed --trials N
                                  --output OUT.tsv [--workers N]
"""

from __future__ import annotations

import argparse
import gzip
import os
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


def _strip_name(filename: str) -> str:
    for ext in (".clean.bed.gz", ".bed.gz", ".clean.bed", ".bed"):
        if filename.endswith(ext):
            return filename[: -len(ext)]
    return filename


def _run_one(args_tuple: tuple[str, str, str, int]) -> tuple[str, tuple | None]:
    query, bed_gz, universe, trials = args_tuple
    name = _strip_name(Path(bed_gz).name)

    with tempfile.NamedTemporaryFile(suffix=".bed", delete=False) as tmp:
        tmp_path = tmp.name
        with gzip.open(bed_gz, "rb") as f:
            tmp.write(f.read())

    try:
        proc = subprocess.run(
            ["bits_test", "-a", query, "-b", tmp_path, "-g", universe, "-n", str(trials)],
            capture_output=True,
            text=True,
        )
        line = proc.stdout.strip()
        if not line:
            return name, None
        # Parse: "O:10  E:4.950000  sd:3.007748  p:0.076923"
        parts: dict[str, float] = {}
        for token in line.split():
            k, v = token.split(":", 1)
            parts[k] = float(v)
        return name, (parts["O"], parts["E"], parts["sd"], parts["p"])
    except Exception:
        return name, None
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", required=True, help="GWAS query BED file")
    parser.add_argument("--beds-dir", required=True, help="Directory of RME *.bed.gz files")
    parser.add_argument("--universe", required=True, help="Universe/whitelist BED for BITS")
    parser.add_argument("--trials", type=int, default=1000, help="Monte Carlo iterations")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("--workers", type=int, default=os.cpu_count(), help="Parallel workers")
    args = parser.parse_args()

    beds = sorted(Path(args.beds_dir).glob("*.bed.gz"))
    if not beds:
        raise SystemExit(f"No *.bed.gz files found in {args.beds_dir}")

    tasks = [(args.query, str(bed), args.universe, args.trials) for bed in beds]

    results: dict[str, tuple] = {}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(_run_one, t): t for t in tasks}
        for future in as_completed(futures):
            name, result = future.result()
            if result is not None:
                results[name] = result

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as out:
        out.write("name\tobserved\texpected\tsd\tp_value\n")
        for name in sorted(results):
            O, E, sd, p = results[name]
            out.write(f"{name}\t{O}\t{E}\t{sd}\t{p}\n")


if __name__ == "__main__":
    main()
