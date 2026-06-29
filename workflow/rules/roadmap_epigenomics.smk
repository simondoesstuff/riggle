# Roadmap Epigenomics (RME) data preparation
#
# Source: 127 epigenomes × 15 chromatin states from the WUSTL coreMarks model.
# Each source file is an EDACC-named bed.gz with all 15 states interleaved.
# Coordinates are hg38 (liftover already applied at source by WUSTL).
#
# This workflow produces:
#   data/roadmap_epigenomics/beds/{cellType}_{chromatinState}.bed.gz   (1905 files)
#   data/roadmap_epigenomics/whitelist.bed             — all-state union, ±10 kbp slop, merged
#   data/roadmap_epigenomics/enhancers_whitelist.bed   — enhancer states only, ±10 kbp slop, merged
#   data/roadmap_epigenomics/quiescent_low_blacklist.bed — Quiescent_Low only, no slop, merged
#
# Processing per bed:
#   1. Remove intervals overlapping hg38 blacklist v2 (Boyle Lab)
#   2. Drop intervals shorter than MIN_SIZE bp
#   3. Sort and bgzip
#
# Expected totals after processing: ~1905 files, ~55M intervals
#
# hg19 alternative: swap DATA_URL for the non-lift tarball:
#   http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/
#     ChmmModels/coreMarks/jointModel/final/all.mnemonics.bedFiles.tgz

import os

RME_DIR = "data/roadmap_epigenomics"
RME_BEDS = f"{RME_DIR}/beds"

RME_DATA_URL = (
    "http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations"
    "/ChmmModels/coreMarks/jointModel/final/all_hg38lift.mnemonics.bedFiles.tgz"
)
RME_BLACKLIST_URL = (
    "https://github.com/Boyle-Lab/Blacklist/raw/refs/heads/master/lists"
    "/hg38-blacklist.v2.bed.gz"
)
RME_CHROM_SIZES_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
)

RME_MIN_SIZE = 190
RME_WHITELIST_SLOP = 10000

# Enhancer state names (subset of the 15 chromatin states)
RME_ENHANCER_STATES = frozenset(["Enhancers", "Genic_enhancers", "Bivalent_Enhancer"])


# --- helper functions for checkpoint-dependent input ---


def rme_all_beds(wildcards):
    """All processed bed.gz files, resolved after the split checkpoint."""
    split_dir = checkpoints.rme_split.get(**wildcards).output[0]
    names = glob_wildcards(os.path.join(split_dir, "{name}.bed")).name
    return expand(f"{RME_BEDS}/{{name}}.bed.gz", name=names)


def rme_enhancer_beds(wildcards):
    """Processed bed.gz files for enhancer chromatin states only."""
    split_dir = checkpoints.rme_split.get(**wildcards).output[0]
    names = glob_wildcards(os.path.join(split_dir, "{name}.bed")).name
    enhancer_names = [
        n for n in names if any(s in n for s in RME_ENHANCER_STATES)
    ]
    return expand(f"{RME_BEDS}/{{name}}.bed.gz", name=enhancer_names)


def rme_quiescent_beds(wildcards):
    """Processed bed.gz files for the Quiescent_Low state only."""
    split_dir = checkpoints.rme_split.get(**wildcards).output[0]
    names = glob_wildcards(os.path.join(split_dir, "{name}.bed")).name
    return expand(
        f"{RME_BEDS}/{{name}}.bed.gz",
        name=[n for n in names if "Quiescent_Low" in n],
    )


# --- rules ---


rule rme_download_raw:
    output:
        temp(f"{RME_BEDS}/raw.tgz"),
    params:
        url=RME_DATA_URL,
    shell:
        "wget -q --show-progress '{params.url}' -O {output}"


rule rme_download_blacklist:
    output:
        temp(f"{RME_BEDS}/hg38-blacklist.v2.bed.gz"),
    params:
        url=RME_BLACKLIST_URL,
    shell:
        "wget -q --show-progress '{params.url}' -O {output}"


rule rme_blacklist_decompress:
    input:
        f"{RME_BEDS}/hg38-blacklist.v2.bed.gz",
    output:
        temp(f"{RME_BEDS}/hg38-blacklist.v2.bed"),
    shell:
        "gunzip -c {input} > {output}"


rule rme_extract_raw:
    input:
        f"{RME_BEDS}/raw.tgz",
    output:
        temp(directory(f"{RME_BEDS}/raw")),
    shell:
        "mkdir -p {output} && tar -xzf {input} -C {output}"


checkpoint rme_split:
    """Split each EDACC source file into per-(cellType, chromatinState) beds."""
    input:
        f"{RME_BEDS}/raw",
    output:
        directory(f"{RME_BEDS}/split"),
    script:
        "../scripts/rme_split.py"


rule rme_process_bed:
    """Filter blacklisted regions, enforce minimum interval size, sort, bgzip."""
    input:
        bed=f"{RME_BEDS}/split/{{name}}.bed",
        blacklist=f"{RME_BEDS}/hg38-blacklist.v2.bed",
    output:
        f"{RME_BEDS}/{{name}}.bed.gz",
    params:
        min_size=RME_MIN_SIZE,
    shell:
        "bedtools intersect -v -a {input.bed} -b {input.blacklist}"
        " | awk '($3-$2) >= {params.min_size}'"
        " | bedtools sort -i -"
        " | bgzip"
        " > {output}"


rule rme_chrom_sizes:
    output:
        f"{RME_DIR}/hg38.chrom.sizes",
    params:
        url=RME_CHROM_SIZES_URL,
    shell:
        "wget -q '{params.url}' -O {output}"


rule rme_whitelist:
    """Non-overlapping union of all RME intervals with ±10 kbp slop.
    Restricts analyses to genomic regions with any epigenomic annotation
    across the full RME compendium. Expected: ~353 merged intervals, ~2.83 Gb covered.
    """
    input:
        beds=rme_all_beds,
        chrom_sizes=f"{RME_DIR}/hg38.chrom.sizes",
    output:
        f"{RME_DIR}/whitelist.bed",
    params:
        slop=RME_WHITELIST_SLOP,
    shell:
        "zcat {input.beds}"
        " | cut -f1-3"
        " | LC_ALL=C sort -k1,1 -k2,2n -S 2G -T {resources.tmpdir}"
        " | bedtools slop -i stdin -g {input.chrom_sizes} -b {params.slop}"
        " | bedtools merge -i stdin"
        " > {output}"


rule rme_enhancers_whitelist:
    """Non-overlapping union of Enhancers, Genic_enhancers, and Bivalent_Enhancer
    intervals with ±10 kbp slop. Expected: ~6858 merged intervals, ~2.72 Gb covered.
    """
    input:
        beds=rme_enhancer_beds,
        chrom_sizes=f"{RME_DIR}/hg38.chrom.sizes",
    output:
        f"{RME_DIR}/enhancers_whitelist.bed",
    params:
        slop=RME_WHITELIST_SLOP,
    shell:
        "zcat {input.beds}"
        " | cut -f1-3"
        " | LC_ALL=C sort -k1,1 -k2,2n -S 2G -T {resources.tmpdir}"
        " | bedtools slop -i stdin -g {input.chrom_sizes} -b {params.slop}"
        " | bedtools merge -i stdin"
        " > {output}"


rule rme_quiescent_blacklist:
    """Non-overlapping union of all Quiescent_Low state intervals (no slop).
    Excludes constitutively quiescent regions from downstream scoring.
    Expected: ~17641 merged intervals, ~2.77 Gb covered.
    """
    input:
        rme_quiescent_beds,
    output:
        f"{RME_DIR}/quiescent_low_blacklist.bed",
    shell:
        "zcat {input}"
        " | cut -f1-3"
        " | LC_ALL=C sort -k1,1 -k2,2n -S 2G -T {resources.tmpdir}"
        " | bedtools merge -i stdin"
        " > {output}"
