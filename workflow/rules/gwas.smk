# GWAS loci data (Maurano et al. / Kundaje et al., Nature 2015, doi:10.1038/nature13835)
#
# Downloads Supplementary Table 1 as XLS, splits into per-trait hg19 BED files,
# then liftovers each to hg38.  Final artifacts: data/gwas/{trait}.bed
#
# Regenerate all GWAS beds:
#   snakemake gwas_all

GWAS_DIR = "data/gwas"
GWAS_RAW = f"{GWAS_DIR}/raw"
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
