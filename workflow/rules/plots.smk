# Scatter plot rules for tool comparison.
# Outputs go to data/plots/.
#
# Usage example:
#   snakemake data/plots/scatter_bed_pvals_Rheumatoid_arthritis_1000.png

GWAS_TRAITS = glob_wildcards("data/gwas/{trait}.bed").trait


rule scatter_bed_pvals:
    """Per-bed -log10(p) scatter for one trait: all pairwise tool comparisons."""
    input:
        chuckle=f"{CHUCKLE_DIR}/{{trait}}.json",
        giggle="data/giggle/{trait}.tsv",
        bits="data/bits/{trials}/{trait}.tsv",
    output:
        "data/plots/scatter_bed_pvals_{trait}_{trials}.png",
    shell:
        "uv run workflow/scripts/tool_scatter.py"
        " --per-bed {wildcards.trait}"
        " --chuckle {CHUCKLE_DIR}"
        " --giggle data/giggle"
        " --bits data/bits/{wildcards.trials}"
        " -o {output} --no-show"


rule scatter_mc:
    """Per-bed -log10(p) scatter replacing BITS with Monte Carlo (10k trials).

    Three panels: chuckle vs giggle, chuckle vs MC, giggle vs MC.
    """
    input:
        chuckle=f"{CHUCKLE_DIR}/{{trait}}.json",
        giggle="data/giggle/{trait}.tsv",
        mc=f"{MC_DIR}/{{trait}}.tsv",
    output:
        "data/plots/scatter_mc_{trait}.png",
    params:
        mc_trials=10000,
    shell:
        "uv run workflow/scripts/tool_scatter.py"
        " --per-bed {wildcards.trait}"
        " --chuckle {CHUCKLE_DIR}"
        " --giggle data/giggle"
        " --mc {MC_DIR}"
        " --mc-trials {params.mc_trials}"
        " -o {output} --no-show"


rule mc_heatmap:
    """Heatmap of chuckle vs MC Spearman correlation across p-value thresholds and trial counts."""
    input:
        chuckle=f"{CHUCKLE_DIR}/{{trait}}.json",
        mc=f"{MC_DIR}/{{trait}}.tsv",
    output:
        "data/plots/mc_heatmap_{trait}.png",
    shell:
        "uv run workflow/scripts/mc_heatmap.py"
        " --chuckle {input.chuckle}"
        " --mc {input.mc}"
        " --trait {wildcards.trait}"
        " -o {output} --no-show"
