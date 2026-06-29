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
