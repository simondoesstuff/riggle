"""At-scale benchmark trial — Snakemake script: protocol.

Times index build, query, and measures index size for Riggle, IGD, and Giggle
on a single (src_files, src_records, qry_records) configuration.  Writes a
structured JSON to snakemake.output[0] so results can be aggregated and plotted
without parsing terminal text.
"""

import json
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

src_files = int(snakemake.wildcards.src_files)
src_records = int(snakemake.wildcards.src_records)
qry_records = int(snakemake.wildcards.qry_records)
timeout = snakemake.params.timeout
batch_size = snakemake.params.batch_size


def run_timed(cmd: list[str], timeout_s: int) -> tuple[float | None, bool, bool]:
    """Run cmd and return (elapsed_seconds, timed_out, failed)."""
    try:
        t0 = time.perf_counter()
        result = subprocess.run(cmd, timeout=timeout_s, capture_output=True)
        elapsed = time.perf_counter() - t0
        if result.returncode != 0:
            print(f"  FAILED ({result.returncode}): {' '.join(cmd[:4])}")
            print(result.stderr.decode(errors="replace")[-500:])
            return None, False, True
        return elapsed, False, False
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT: {' '.join(cmd[:4])}")
        return None, True, False
    except FileNotFoundError:
        print(f"  NOT FOUND: {cmd[0]}")
        return None, False, True


def dir_bytes(path: str) -> int:
    total = 0
    for dirpath, _, files in os.walk(path):
        for name in files:
            fp = os.path.join(dirpath, name)
            if not os.path.islink(fp):
                total += os.path.getsize(fp)
    return total


work = tempfile.mkdtemp(prefix="bench_")
print(f"Working in {work}")
try:
    src_dir = os.path.join(work, "src")
    qry_file = os.path.join(work, "query.bed.gz")
    os.makedirs(src_dir)

    print(f"Generating {src_files} source beds × {src_records} records each...")
    subprocess.run(
        [
            "cargo", "run", "--release", "--bin", "gen", "--",
            "-f", str(src_files), "-n", str(src_records),
            "--min-len", "1", "--max-len", "50000", "-c", "--sort", "-o", src_dir,
        ],
        check=True, capture_output=True,
    )

    print(f"Generating query bed with {qry_records} records...")
    subprocess.run(
        [
            "cargo", "run", "--release", "--bin", "gen", "--",
            "-n", str(qry_records), "-s", "99",
            "--min-len", "1", "--max-len", "50000", "-c", "--sort", "-o", qry_file,
        ],
        check=True, capture_output=True,
    )

    results: dict[str, dict] = {}

    # --- Chuckle ---
    print("Chuckle: building index...")
    idx = os.path.join(work, "riggle_idx")
    idx_s, idx_to, idx_fail = run_timed(
        [
            "cargo", "run", "--release", "--bin", "riggle", "--",
            "add", "-i", src_dir, "-d", idx, "--batch-size", str(batch_size),
        ],
        timeout,
    )
    idx_bytes = dir_bytes(idx) if os.path.exists(idx) else None
    if not idx_fail and not idx_to and os.path.exists(idx):
        print("Chuckle: querying...")
        qry_s, qry_to, _ = run_timed(
            [
                "cargo", "run", "--release", "--bin", "riggle", "--",
                "query", "-d", idx, "-q", qry_file, "-o", "/dev/null",
                "--batch-size", str(batch_size),
            ],
            timeout,
        )
    else:
        qry_s, qry_to = None, False
    results["Chuckle"] = {
        "index_s": idx_s, "query_s": qry_s, "index_bytes": idx_bytes,
        "index_timeout": idx_to, "query_timeout": qry_to,
    }

    # --- IGD ---
    print("IGD: building index...")
    idx = os.path.join(work, "igd_idx")
    os.makedirs(idx, exist_ok=True)
    idx_s, idx_to, idx_fail = run_timed(
        ["igd", "create", src_dir, idx, "-s", "0"],
        timeout,
    )
    idx_bytes = dir_bytes(idx) if os.path.exists(idx) else None
    igd_files = list(Path(idx).rglob("*.igd")) if os.path.exists(idx) else []
    if not idx_fail and not idx_to and igd_files:
        print("IGD: querying...")
        qry_s, qry_to, _ = run_timed(
            ["igd", "search", str(igd_files[0]), "-q", qry_file],
            timeout,
        )
    else:
        qry_s, qry_to = None, False
    results["IGD"] = {
        "index_s": idx_s, "query_s": qry_s, "index_bytes": idx_bytes,
        "index_timeout": idx_to, "query_timeout": qry_to,
    }

    # --- Giggle ---
    print("Giggle: building index...")
    idx = os.path.join(work, "giggle_idx")
    idx_s, idx_to, idx_fail = run_timed(
        ["giggle", "index", "-s", "-f", "-i", f"{src_dir}/*", "-o", idx],
        timeout,
    )
    idx_bytes = dir_bytes(idx) if os.path.exists(idx) else None
    giggle_ok = os.path.exists(os.path.join(idx, "root_ids.dat"))
    if not idx_fail and not idx_to and giggle_ok:
        print("Giggle: querying...")
        qry_s, qry_to, _ = run_timed(
            ["giggle", "search", "-i", idx, "-q", qry_file],
            timeout,
        )
    else:
        qry_s, qry_to = None, False
    results["Giggle"] = {
        "index_s": idx_s, "query_s": qry_s, "index_bytes": idx_bytes,
        "index_timeout": idx_to, "query_timeout": qry_to,
    }

    output = {
        "src_files": src_files,
        "src_records": src_records,
        "qry_records": qry_records,
        "total_intervals": src_files * src_records,
        "results": results,
    }
    Path(snakemake.output[0]).parent.mkdir(parents=True, exist_ok=True)
    with open(snakemake.output[0], "w") as f:
        json.dump(output, f, indent=2)
    print(f"Written: {snakemake.output[0]}")

finally:
    shutil.rmtree(work, ignore_errors=True)
