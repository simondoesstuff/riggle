"""
eigenbed_score.py — eigenbed-based cosine similarity scoring

Usage (build index from corpus):
    eigenbed_score.py --corpus A.bed B.bed ... --reference R1.bed R2.bed ... --query Q.bed [options]

Usage (reuse a pre-built index):
    eigenbed_score.py --index DIR --reference R1.bed R2.bed ... --query Q.bed [options]

Steps:
  1. Build an eigenbed index from the corpus BED files  (skipped if --index points to
     an existing index and --corpus is omitted).
  2. Embed the reference BED files through that index.
  3. Embed the query BED file through that index.
  4. Report cosine similarity between the query and each reference embedding,
     sorted descending.

Options:
  --corpus BED [...]      BED files to build the index from (omit to reuse --index)
  --reference BED [...]   BED files to score against the query (required)
  --query BED             Single BED file to use as the query (required)
  --index DIR             Directory for the eigenbed index; if --corpus is omitted this
                          must point to an already-built index [default: temp dir]
  --components K          Max SVD components [default: 50]
  --eigenbed PATH         Path to the eigenbed binary [default: auto-detect]
  --top N                 Show only top N results [default: all]
  --keep-index            Do not delete the index directory after scoring
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

# ── helpers ──────────────────────────────────────────────────────────────────


def find_eigenbed(override: str | None) -> str:
    if override:
        return override
    # Prefer a pre-built binary next to this script or on PATH.
    candidates = [
        os.path.join(os.path.dirname(__file__), "..", "target", "release", "eigenbed"),
        os.path.join(os.path.dirname(__file__), "..", "target", "debug", "eigenbed"),
        shutil.which("eigenbed") or "",
    ]
    for c in candidates:
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    # Fall back to cargo run (slow but always works from repo root).
    return None  # signals to use cargo run


def run(cmd: list[str], desc: str) -> str:
    """Run a command, stream stderr to our stderr, return stdout."""
    print(f"[eigenbed] {desc}", file=sys.stderr)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr, end="")
        sys.exit(f"error: command failed: {' '.join(cmd)}")
    # Echo eigenbed's progress lines.
    if result.stderr:
        print(result.stderr, file=sys.stderr, end="")
    return result.stdout


def parse_embeddings(tsv: str) -> dict[str, list[float]]:
    """Parse eigenbed embed TSV output → {path: [f32, ...]}."""
    embeddings: dict[str, list[float]] = {}
    for line in tsv.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        path = parts[0]
        values = [float(v) for v in parts[1:]]
        embeddings[path] = values
    return embeddings


def cosine_similarity(a: list[float], b: list[float]) -> float:
    """Cosine similarity between two vectors.
    Embeddings are already L2-normalised by eigenbed, so this is just the dot product."""
    dot = sum(x * y for x, y in zip(a, b))
    # Clamp to [-1, 1] to guard against tiny floating-point overflows.
    return max(-1.0, min(1.0, dot))


# ── main ─────────────────────────────────────────────────────────────────────


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--corpus", nargs="+", default=None, metavar="BED")
    p.add_argument("--reference", nargs="+", required=True, metavar="BED")
    p.add_argument("--query", required=True, metavar="BED")
    p.add_argument("--index", default=None, metavar="DIR")
    p.add_argument("--components", type=int, default=50, metavar="K")
    p.add_argument("--eigenbed", default=None, metavar="PATH")
    p.add_argument("--top", type=int, default=None, metavar="N")
    p.add_argument("--keep-index", action="store_true")
    args = p.parse_args()

    # Validate: need either --corpus (to build) or --index (pre-built).
    if not args.corpus and not args.index:
        p.error("one of --corpus or --index (pre-built) is required")

    # Resolve eigenbed binary / cargo-run fallback.
    eb = find_eigenbed(args.eigenbed)

    def eigenbed_cmd(subcmd: list[str]) -> list[str]:
        if eb:
            return [eb] + subcmd
        # Use cargo run from the repo root.
        return ["cargo", "run", "--release", "--bin", "eigenbed", "--"] + subcmd

    # Index directory: temp dir when building from corpus without --index;
    # or the user-specified dir (pre-built or newly built).
    own_tmpdir = None
    if args.index:
        index_dir = args.index
    else:
        own_tmpdir = tempfile.mkdtemp(prefix="eigenbed_pca_")
        index_dir = own_tmpdir

    try:
        # 1. Build index from corpus (skipped when reusing a pre-built index).
        if args.corpus:
            run(
                eigenbed_cmd(
                    [
                        "index",
                        "--output",
                        index_dir,
                        "--components",
                        str(args.components),
                    ]
                    + args.corpus
                ),
                desc=f"indexing {len(args.corpus)} corpus BED file(s) → {index_dir}",
            )
        else:
            print(f"[eigenbed] reusing pre-built index at {index_dir}", file=sys.stderr)

        # 2. Embed reference BED files.
        ref_tsv = run(
            eigenbed_cmd(["embed", "--index", index_dir] + args.reference),
            desc=f"embedding {len(args.reference)} reference BED file(s)",
        )
        ref_embeddings = parse_embeddings(ref_tsv)

        # 3. Embed query BED file.
        qry_tsv = run(
            eigenbed_cmd(["embed", "--index", index_dir, args.query]),
            desc=f"embedding query: {args.query}",
        )
        qry_embeddings = parse_embeddings(qry_tsv)

        if not qry_embeddings:
            sys.exit("error: no embedding produced for query")

        query_path, query_vec = next(iter(qry_embeddings.items()))

        # 4. Compute and report cosine similarities.
        scores: list[tuple[float, str]] = []
        for ref_path, ref_vec in ref_embeddings.items():
            sim = cosine_similarity(query_vec, ref_vec)
            scores.append((sim, ref_path))

        scores.sort(key=lambda x: -x[0])
        if args.top is not None:
            scores = scores[: args.top]

        print(f"# query: {query_path}")
        print("# reference_file\tcosine_similarity")
        for sim, ref_path in scores:
            print(f"{ref_path}\t{sim:.6f}")

    finally:
        if own_tmpdir and not args.keep_index:
            shutil.rmtree(own_tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
