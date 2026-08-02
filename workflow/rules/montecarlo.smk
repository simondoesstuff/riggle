# Monte Carlo shift permutation test against the RME chuckle index.
#
# Intermediate data: data/montecarlo/{trait}.tsv
#   Multi-milestone TSV; each milestone section is delimited by a
#   "# milestone=N trials=N" header line, then a tab-separated data block.
#
# Usage:
#   snakemake data/montecarlo/Rheumatoid_arthritis.tsv
#   snakemake data/plots/scatter_mc_Rheumatoid_arthritis.png
#   snakemake data/plots/mc_heatmap_Rheumatoid_arthritis.png

MC_DIR = "data/montecarlo"
MC_BIN = "target/release/montecarlo"
MC_MILESTONES = "100,500,1000,5000,10000"
MC_BATCH_SIZE = 100
MC_SEED = 42


rule montecarlo_build:
    """Ensure the montecarlo release binary is up to date."""
    input:
        ancient("Cargo.toml"),
        ancient("src/bin/montecarlo.rs"),
    output:
        MC_BIN,
    shell:
        "cargo build --release --bin montecarlo"


rule montecarlo_query:
    """Run Monte Carlo permutation test for one trait across the RME index.

    Outputs a multi-milestone TSV with p-values at each trial milestone.
    No statistics (FFT) are computed during trials — overlap counts only.
    The whitelist restricts shifted intervals to annotated regions, matching
    chuckle's null model.
    """
    input:
        query="data/queries/{trait}.bed",
        bin=ancient(MC_BIN),
        index=CHUCKLE_INDEX,
        whitelist=ancient(f"{RME_DIR}/whitelist.bed"),
    output:
        f"{MC_DIR}/{{trait}}.tsv",
    params:
        milestones=MC_MILESTONES,
        batch_size=MC_BATCH_SIZE,
        seed=MC_SEED,
    shell:
        "{input.bin}"
        " --db {CHUCKLE_INDEX}"
        " --query {input.query}"
        " --whitelist {input.whitelist}"
        " --trials {params.milestones}"
        " --batch-size {params.batch_size}"
        " --seed {params.seed}"
        " > {output}"
