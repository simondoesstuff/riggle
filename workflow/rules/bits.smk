# BITS Monte Carlo interval intersection significance testing
#
# Outputs data/bits/{trials}/{trait}.tsv  — one row per RME (cell×state) bed.
# `trials` is the number of Monte Carlo permutations.  Use 100 for exploration,
# 1000 for publication-quality p-values.  Results at different trial counts are
# independent artifacts so both can coexist under different subdirectories.
#
# Each bits_test call uses ~5 MB RAM; the sweep script parallelises internally
# across beds for a given trait.  Multiple traits can safely run in parallel.

import os

BITS_DIR = "data/bits"


def _bits_all_beds(wildcards):
    split_dir = checkpoints.rme_split.get(**wildcards).output[0]
    names = glob_wildcards(os.path.join(split_dir, "{name}.bed")).name
    return expand(f"{RME_BEDS}/{{name}}.bed.gz", name=names)


rule bits_sweep:
    """BITS Monte Carlo sweep: one GWAS trait × all RME beds → TSV."""
    input:
        query="data/gwas/{trait}.bed",
        beds=_bits_all_beds,
        universe=f"{RME_DIR}/whitelist.bed",
    output:
        f"{BITS_DIR}/{{trials}}/{{trait}}.tsv",
    threads: 8
    params:
        beds_dir=RME_BEDS,
    shell:
        "uv run workflow/scripts/bits_sweep.py"
        " --query {input.query}"
        " --beds-dir {params.beds_dir}"
        " --universe {input.universe}"
        " --trials {wildcards.trials}"
        " --output {output}"
        " --workers {threads}"
