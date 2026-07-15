# regioners permutation test (shuffle / circle / novl) against all RME beds.
#
# Outputs: data/regioners/{method}/{trait}.tsv  — one row per RME (cell×state) bed.
# shuffle and circle run at 1000 trials; novl at 100 (it is ~10x slower).
#
# Usage:
#   snakemake -j8 regioners_all
#   snakemake -j8 data/regioners/shuffle/Rheumatoid_arthritis.tsv

import os
import shutil
from glob import glob
from pathlib import Path

REGIONERS_DIR = "data/regioners"
REGIONERS_METHODS = ["shuffle", "circle", "novl"]
REGIONERS_TRIALS = {"shuffle": 1000, "circle": 1000, "novl": 100}
REGIONERS_TRAITS = [Path(f).stem for f in glob("data/gwas/*.bed")]


def _regioners_bin(_):
    path = shutil.which("regioners")
    if path is None:
        raise ValueError("regioners not found in PATH — run from within the nix dev shell")
    return path


def _regioners_all_beds(wildcards):
    split_dir = checkpoints.rme_split.get(**wildcards).output[0]
    names = glob_wildcards(os.path.join(split_dir, "{name}.bed")).name
    return expand(f"{RME_BEDS}/{{name}}.bed.gz", name=names)


def _regioners_trials(wildcards):
    return REGIONERS_TRIALS[wildcards.method]


rule regioners_all:
    """All GWAS traits × all three randomization methods."""
    input:
        expand(
            f"{REGIONERS_DIR}/{{method}}/{{trait}}.tsv",
            method=REGIONERS_METHODS,
            trait=REGIONERS_TRAITS,
        ),


rule regioners_sweep:
    """regioners permutation test: one GWAS trait × all RME beds × one method → TSV."""
    input:
        query="data/queries/{trait}.bed",
        beds=_regioners_all_beds,
        genome=f"{RME_DIR}/hg38.chrom.sizes",
    output:
        f"{REGIONERS_DIR}/{{method}}/{{trait}}.tsv",
    threads: 8
    params:
        beds_dir=RME_BEDS,
        trials=_regioners_trials,
        bin=_regioners_bin,
    wildcard_constraints:
        method="shuffle|circle|novl",
    shell:
        "uv run workflow/scripts/regioners_sweep.py"
        " --query {input.query}"
        " --beds-dir {params.beds_dir}"
        " --genome {input.genome}"
        " --method {wildcards.method}"
        " --trials {params.trials}"
        " --output {output}"
        " --workers {threads}"
        " --bin {params.bin}"
