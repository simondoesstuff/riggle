"""
Run regioners for one GWAS query against all RME beds with a given randomization method.

Usage:
    uv run workflow/scripts/regioners_sweep.py --query QUERY.bed --beds-dir DIR
        --genome CHROM_SIZES --method METHOD --trials N --output OUT.tsv [--workers N]
"""

from __future__ import annotations

import argparse
import gzip
import json
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


def _run_one(args_tuple: tuple) -> tuple[str, dict | None]:
    query, bed_gz, genome, method, trials, bin_path = args_tuple
    name = _strip_name(Path(bed_gz).name)

    with tempfile.NamedTemporaryFile(suffix=".bed", delete=False) as tmp_b:
        tmp_b_path = tmp_b.name
        with gzip.open(bed_gz, "rb") as f:
            tmp_b.write(f.read())

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tmp_out:
        tmp_out_path = tmp_out.name

    try:
        proc = subprocess.run(
            [
                bin_path,
                "--genome", genome,
                "-A", query,
                "-B", tmp_b_path,
                "--num-times", str(trials),
                "--random", method,
                "--output", tmp_out_path,
            ],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            return name, None
        with open(tmp_out_path) as f:
            data = json.load(f)
        test = data["test"]
        return name, {
            "observed": test["observed"],
            "mean": test["mean"],
            "std_dev": test["std_dev"],
            "z_score": test["z_score"],
            "p_val": test["p_val"],
            "alt": test["alt"],
        }
    except json.JSONDecodeError:
        return name, None
    finally:
        Path(tmp_b_path).unlink(missing_ok=True)
        Path(tmp_out_path).unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", required=True)
    parser.add_argument("--beds-dir", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--method", required=True, choices=["shuffle", "circle", "novl"])
    parser.add_argument("--trials", type=int, default=100)
    parser.add_argument("--output", required=True)
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--bin", default="regioners", dest="bin_path")
    args = parser.parse_args()

    # Fail fast if the binary isn't callable.
    try:
        subprocess.run([args.bin_path, "--version"], capture_output=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        raise SystemExit(f"regioners binary not found or broken: {args.bin_path!r} ({e})")

    beds = sorted(Path(args.beds_dir).glob("*.bed.gz"))
    if not beds:
        raise SystemExit(f"No *.bed.gz files found in {args.beds_dir}")

    tasks = [(args.query, str(bed), args.genome, args.method, args.trials, args.bin_path) for bed in beds]

    results: dict[str, dict] = {}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(_run_one, t): t for t in tasks}
        for future in as_completed(futures):
            name, result = future.result()
            if result is not None:
                results[name] = result

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as out:
        out.write("name\tobserved\tmean\tstd_dev\tz_score\tp_val\talt\n")
        for name in sorted(results):
            r = results[name]
            out.write(f"{name}\t{r['observed']}\t{r['mean']}\t{r['std_dev']}\t{r['z_score']}\t{r['p_val']}\t{r['alt']}\n")


if __name__ == "__main__":
    main()
