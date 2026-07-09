# Scatter plot rules for tool comparison.
# Outputs go to data/plots/.
#
# Usage example:
#   snakemake data/plots/scatter_bed_pvals_Rheumatoid_arthritis_1000.png

GWAS_TRAITS = glob_wildcards("data/gwas/{trait}.bed").trait

RME_HEAT_TARGETS = [
    "data/plots/rme_heat_Rheumatoid_arthritis.png",
    "data/plots/rme_heat_Alzheimers_combined.png",
    "data/plots/rme_heat_myod_myotube.png",
]


rule rme_heat_all:
    """Produce RME heatmaps for the default trait set."""
    input:
        RME_HEAT_TARGETS,


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


# INFO: ------------------------
# 			  RME Heatmaps
# ------------------------------


RME_HEAT_STATES = [
    "Enhancers",
    "Genic_enhancers",
    "Flanking_Active_TSS",
    "Bivalent_Enhancer",
    "Repressed_PolyComb",
    "Weak_Repressed_PolyComb",
    "Flanking_Bivalent_TSS_Enh",
    "Bivalent_Poised_TSS",
    "Strong_transcription",
    "ZNF_genes_and_repeats",
]


rule rme_heat:
    """Side-by-side RME heatmap: chuckle (−log p-value) and giggle (combo score).

    Tissue types on the vertical axis (grouped by cell category), chromatin
    states on the horizontal axis.  Any GWAS trait with both chuckle and giggle
    outputs is supported via the wildcard.

    Example:
        snakemake data/plots/rme_heat_Rheumatoid_arthritis.png
        snakemake data/plots/rme_heat_Alzheimers_combined.png
    """
    input:
        chuckle=f"{CHUCKLE_DIR}/{{trait}}.json",
        giggle="data/giggle/{trait}.tsv",
    output:
        "data/plots/rme_heat_{trait}.png",
    params:
        states=" ".join(RME_HEAT_STATES),
    shell:
        "uv run workflow/scripts/rme_heat.py"
        " --scores {input.chuckle} {input.giggle}"
        " --names 'Chuckle (LLR)' 'Giggle (combo score)'"
        " --score-fields llr combo_score"
        " --states {params.states}"
        " -o {output} --no-show"


rule rme_heat_full:
    """Four-panel RME heatmap: giggle p-value, giggle combo score,
    chuckle p-value, chuckle LLR (NB saddlepoint).

    Columns restricted to the enhancer-centric state subset in RME_HEAT_STATES.
    Requires a chuckle JSON built with the current binary (llr field present).

    Example:
        snakemake data/plots/rme_heat_full_Rheumatoid_arthritis.png
    """
    input:
        chuckle=f"{CHUCKLE_DIR}/{{trait}}.json",
        giggle="data/giggle/{trait}.tsv",
    output:
        "data/plots/rme_heat_full_{trait}.png",
    params:
        states=" ".join(RME_HEAT_STATES),
    shell:
        "uv run workflow/scripts/rme_heat.py"
        " --scores {input.giggle} {input.giggle} {input.chuckle} {input.chuckle}"
        " --names 'Giggle (p-value)' 'Giggle (combo score)' 'Chuckle (p-value)' 'Chuckle (LLR)'"
        " --score-fields fishers_right_tail combo_score p_value llr"
        " --states {params.states}"
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
