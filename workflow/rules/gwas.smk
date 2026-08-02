# GWAS loci data (Maurano et al. / Kundaje et al., Nature 2015, doi:10.1038/nature13835)
#
# Downloads Supplementary Table 1 as XLS, splits into per-trait hg19 BED files,
# then liftovers each to hg38.  Final artifacts: data/gwas/{trait}.bed
#
# Disease→tissue mapping (Disease Ontology):
#   snakemake data/gwas/disease_tissue.json
#
# Regenerate all GWAS beds:
#   snakemake gwas_all
#
# ROC-AUC benchmark (giggle + chuckle vs disease-tissue ground truth):
#   snakemake data/gwas/roc_data.json
#   snakemake data/plots/gwas_roc.png

GWAS_DIR = "data/gwas"
GWAS_RAW = f"{GWAS_DIR}/raw"
DOID_DIR = "data/doid"
DOID_OBO_URL = "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/doid.obo"
GWAS_XLS_URL = (
    "https://www.nature.com/nature/journal/v518/n7539/extref/nature13835-s1.xls"
)
GWAS_CHAIN_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
)


rule gwas_download_xls:
    """Download GWAS supplementary XLS from Nature."""
    output:
        f"{GWAS_RAW}/nature13835-s1.xls",
    params:
        url=GWAS_XLS_URL,
    shell:
        "mkdir -p {GWAS_RAW} && wget -q --show-progress '{params.url}' -O {output}"


rule gwas_download_chain:
    """Download hg19→hg38 liftOver chain file."""
    output:
        f"{GWAS_RAW}/hg19ToHg38.over.chain.gz",
    params:
        url=GWAS_CHAIN_URL,
    shell:
        "wget -q '{params.url}' -O {output}"


checkpoint gwas_split:
    """Parse GWAS XLS and write per-trait hg19 BED files."""
    input:
        f"{GWAS_RAW}/nature13835-s1.xls",
    output:
        directory(f"{GWAS_DIR}/hg19"),
    shell:
        "uv run workflow/scripts/gwas_split.py {input} {output}"


rule gwas_liftover:
    """Liftover one GWAS trait BED from hg19 to hg38."""
    input:
        bed=f"{GWAS_DIR}/hg19/{{trait}}.bed",
        chain=f"{GWAS_RAW}/hg19ToHg38.over.chain.gz",
    output:
        f"{GWAS_DIR}/{{trait}}.bed",
    shell:
        "liftOver {input.bed} {input.chain} {output} /dev/null"


def _all_gwas_beds(wildcards):
    """Resolve all lifted GWAS beds after the split checkpoint runs."""
    split_dir = checkpoints.gwas_split.get(**wildcards).output[0]
    names = glob_wildcards(f"{split_dir}/{{name}}.bed").name
    return expand(f"{GWAS_DIR}/{{name}}.bed", name=names)


rule gwas_all:
    """Build all hg38 GWAS trait BED files."""
    input:
        _all_gwas_beds,


rule gwas_download_doid:
    """Download Disease Ontology OBO (stored outside data/gwas)."""
    output:
        f"{DOID_DIR}/doid.obo",
    params:
        url=DOID_OBO_URL,
    shell:
        "mkdir -p {DOID_DIR} && wget -q --show-progress '{params.url}' -O {output}"


rule gwas_disease_tissue:
    """Map each GWAS trait to macro tissue categories via DOID is_a hierarchy."""
    input:
        doid=f"{DOID_DIR}/doid.obo",
        hg19_dir=f"{GWAS_DIR}/hg19",
    output:
        f"{GWAS_DIR}/disease_tissue.json",
    shell:
        "uv run workflow/scripts/gwas_disease_tissue.py {input.doid} {input.hg19_dir} {output}"


# ---------------------------------------------------------------------------
# ROC-AUC benchmark
#
# Giggle results are available for all GWAS diseases; chuckle results are
# picked up opportunistically from data/chuckle/hg38/ (expensive to generate —
# see chuckle.smk / chuckle_gwas_all).  The data artifact is rebuilt whenever
# giggle results change; re-run manually after adding new chuckle results.
# ---------------------------------------------------------------------------

import json as _json
from pathlib import Path as _Path

# Traits with ground-truth tissue labels and existing giggle results
_dt_path = f"{GWAS_DIR}/disease_tissue.json"
if _Path(_dt_path).exists():
    with open(_dt_path) as _f:
        _DT = _json.load(_f)
    GWAS_ROC_TRAITS = [
        t for t in _DT
        if _Path(f"data/giggle/{t}.tsv").exists()
    ]
else:
    GWAS_ROC_TRAITS = []


rule gwas_roc_data:
    """Compute per-disease and combined ROC curves for giggle and chuckle."""
    input:
        disease_tissue=f"{GWAS_DIR}/disease_tissue.json",
        giggle=expand("data/giggle/{trait}.tsv", trait=GWAS_ROC_TRAITS),
        chuckle=expand("data/chuckle/hg38/{trait}.json", trait=GWAS_ROC_TRAITS),
    output:
        f"{GWAS_DIR}/roc_data.json",
    params:
        giggle_dir="data/giggle",
        chuckle_dir="data/chuckle/hg38",
    shell:
        "uv run workflow/scripts/gwas_roc.py"
        " --disease-tissue {input.disease_tissue}"
        " --giggle-dir {params.giggle_dir}"
        " --chuckle-dir {params.chuckle_dir}"
        " --output {output}"


rule gwas_roc_plot:
    """Plot ROC-AUC benchmark: combined curve + per-disease AUC bar chart."""
    input:
        f"{GWAS_DIR}/roc_data.json",
    output:
        "data/plots/gwas_roc.png",
    shell:
        "uv run workflow/scripts/gwas_roc_plot.py"
        " --roc-data {input}"
        " --output {output}"
        " --no-show"
