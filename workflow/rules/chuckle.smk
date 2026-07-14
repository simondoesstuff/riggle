# Chuckle (FFT-based interval intersection significance testing)
#
# WARNING: querying the full RME index requires ~250 GB RAM — allocate
# --resources mem_mb=260000 and run at most one query job at a time.
# Run as:
#   snakemake --resources mem_mb=260000 -j1 data/chuckle/<trait>.json
#
# The RME index (data/rme_chuckle) is managed by the chuckle_index rule below.
# Building it is a one-time cost; query jobs depend on it automatically.

import os

CHUCKLE_DIR = "data/chuckle"
CHUCKLE_INDEX = "data/rme_chuckle"
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
    """Build chuckle Fourier-domain index from all processed RME beds.

    This is an expensive one-time build kept as a managed Snakemake artifact
    so query rules can depend on it without rebuilding unnecessarily.
    """
    input:
        beds=rme_all_beds,
        bin=CHUCKLE_BIN,
    output:
        directory(CHUCKLE_INDEX),
    params:
        beds_dir=RME_BEDS,
    shell:
        "{input.bin} add -i {params.beds_dir} -d {output}"


rule chuckle_query:
    """Query chuckle RME index for one GWAS trait → JSON with per-bed p-values."""
    input:
        query="data/queries/{trait}.bed",
        bin=CHUCKLE_BIN,
        index=CHUCKLE_INDEX,
        whitelist=ancient(f"{RME_DIR}/whitelist.bed"),
    output:
        f"{CHUCKLE_DIR}/{{trait}}.json",
    resources:
        mem_mb=260000,
    shell:
        "{input.bin} query"
        " --db {CHUCKLE_INDEX}"
        " --query {input.query}"
        " --output {output}"
        " --stats"
        " --whitelist {input.whitelist}"
