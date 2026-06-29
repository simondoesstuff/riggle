# Giggle interval search index and per-trait significance scoring
#
# Giggle uses Fisher's exact test on interval overlaps.  Output columns:
#   file  file_size  overlaps  odds_ratio  fishers_two_tail
#   fishers_left_tail  fishers_right_tail  combo_score
#
# `fishers_right_tail` is the enrichment p-value (observed >= expected);
# `combo_score` is giggle's composite log-odds enrichment metric.
# Both are retained in the TSV so downstream scripts can choose.
#
# Note: giggle requires bgzipped, sorted input files and bgzipped query files.

import os

GIGGLE_DIR = "data/giggle"
# Giggle requires index and input to share a parent directory.
# Beds are under {RME_DIR}/beds/, so index must also live under {RME_DIR}/.
GIGGLE_INDEX_DIR = f"{RME_DIR}/giggle_index"


def _giggle_all_beds(wildcards):
    split_dir = checkpoints.rme_split.get(**wildcards).output[0]
    names = glob_wildcards(os.path.join(split_dir, "{name}.bed")).name
    return expand(f"{RME_BEDS}/{{name}}.bed.gz", name=names)


rule giggle_index:
    """Build giggle search index from all processed RME beds."""
    input:
        _giggle_all_beds,
    output:
        directory(GIGGLE_INDEX_DIR),
    params:
        beds_glob=f"{RME_BEDS}/*.bed.gz",
    shell:
        "giggle index -i '{params.beds_glob}' -o {output} -s -f"


rule giggle_bgzip_query:
    """Bgzip a GWAS query BED so giggle can read it."""
    input:
        "data/gwas/{trait}.bed",
    output:
        temp(f"{GIGGLE_DIR}/queries/{{trait}}.bed.gz"),
    shell:
        "mkdir -p $(dirname {output}) && bgzip -c {input} > {output}"


rule giggle_search:
    """Search giggle index for one GWAS trait → significance TSV."""
    input:
        query=f"{GIGGLE_DIR}/queries/{{trait}}.bed.gz",
        index=GIGGLE_INDEX_DIR,
    output:
        f"{GIGGLE_DIR}/{{trait}}.tsv",
    params:
        index=GIGGLE_INDEX_DIR,
    shell:
        "giggle search -i {params.index} -q {input.query} -s > {output}"
