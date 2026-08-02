# Chuckle (analytic NB interval intersection significance testing)
#
# Generate all GWAS trait queries (runs in parallel):
#   snakemake -j39 chuckle_gwas_all
#
# The RME index (data/roadmap_epigenomics/chuckle_index) is managed by the chuckle_index rule below.
# Building it is a one-time cost; query jobs depend on it automatically.

import os

CHUCKLE_DIR = "data/chuckle/hg38"
CHUCKLE_INDEX = f"{RME_DIR}/chuckle_index"
CHUCKLE_BIN = "target/release/chuckle"


rule chuckle_build:
    """Ensure the chuckle release binary is up to date."""
    input:
        ancient("Cargo.toml"),
        ancient("src/main.rs"),
    output:
        CHUCKLE_BIN,
    shell:
        "cargo build --release --bin chuckle"


rule chuckle_index:
    """Build chuckle depth-map index from all processed RME beds.

    This is an expensive one-time build kept as a managed Snakemake artifact
    so query rules can depend on it without rebuilding unnecessarily.
    """
    input:
        beds=rme_all_beds,
        bin=ancient(CHUCKLE_BIN),
    output:
        directory(CHUCKLE_INDEX),
    params:
        beds_dir=RME_BEDS,
    shell:
        "{input.bin} add -i {params.beds_dir} -d {output} --stats"


rule chuckle_query:
    """Query chuckle RME index for one GWAS trait → JSON with per-bed p-values."""
    input:
        query="data/queries/{trait}.bed",
        bin=CHUCKLE_BIN,
        index=CHUCKLE_INDEX,
    output:
        f"{CHUCKLE_DIR}/{{trait}}.json",
    shell:
        "{input.bin} query"
        " --db {CHUCKLE_INDEX}"
        " --query {input.query}"
        " --output {output}"
        " --stats"


def _all_chuckle_gwas(wildcards):
    """Chuckle hg38 JSONs for all GWAS traits (resolved after gwas_split checkpoint)."""
    split_dir = checkpoints.gwas_split.get(**wildcards).output[0]
    names = glob_wildcards(f"{split_dir}/{{name}}.bed").name
    return expand(f"{CHUCKLE_DIR}/{{name}}.json", name=names)


rule chuckle_gwas_all:
    """Generate chuckle hg38 queries for all GWAS traits.

    Example:
        snakemake -j39 chuckle_gwas_all
    """
    input:
        _all_chuckle_gwas,
