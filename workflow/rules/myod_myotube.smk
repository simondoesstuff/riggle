# MyoD ChIP-Seq peaks in human myotubes — GEO accession GSM1218850
#
# Source: Sarthy et al., Seattle Children's Research Institute.
#   Part of series GSE50413 / GSE50415.
#   Cell line: MB135 (primary human myoblasts) differentiated into myotubes
#   by serum withdrawal (1% horse serum + insulin + transferring, 5 days).
#   Antibody: MyoD (ChIP-Seq).
#
# Reads were aligned to hg19 with BWA 0.5.9; peaks called with an in-house
# algorithm and deposited in ENCODE broadPeak format.  This rule downloads
# the supplementary peak file, strips to BED3, liftovers to hg38, and
# produces a sorted bgzipped BED.
#
# FTP:  ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1218nnn/GSM1218850/suppl/
#       GSM1218850_MB135DMMD.peak.txt.gz
# GEO:  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1218850
#
# Output:
#   data/myod_myotube/peaks.bed.gz   — hg38 BED3, sorted, bgzipped

MYOD_DIR = "data/myod_myotube"
MYOD_FTP_URL = (
    "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1218nnn/GSM1218850/suppl"
    "/GSM1218850_MB135DMMD.peak.txt.gz"
)

QUERY_BEDS["myod_myotube"] = f"{MYOD_DIR}/peaks.bed.gz"


rule myod_peaks:
    """Download, liftover hg19→hg38, sort, and bgzip MyoD ChIP-Seq peaks."""
    input:
        chain=f"{GWAS_RAW}/hg19ToHg38.over.chain.gz",
    output:
        f"{MYOD_DIR}/peaks.bed.gz",
    params:
        url=MYOD_FTP_URL,
    shell:
        "wget -q -O - '{params.url}'"
        " | zcat"
        " | awk '!/^#/ && $1 ~ /^chr/'"
        " | cut -f1-3"
        " | liftOver stdin {input.chain} /dev/stdout /dev/null"
        " | bedtools sort -i -"
        " | bgzip"
        " > {output}"
