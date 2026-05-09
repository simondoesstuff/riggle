Roadmap epigenomics chmm marks are provided as a series of 127 bed files named by edacc IDs where each region is annotated with one of fifteen chromatin states. To make the data easier to work with, we: split, sort and rename files into the form `{cellType}_{chromatinState}.bed`.

```bash
WORK_DIR="data/roadmap_epigenomics/beds"
DATA_URL="http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all_hg38lift.mnemonics.bedFiles.tgz"
BLACKLIST_URL="https://github.com/Boyle-Lab/Blacklist/raw/refs/heads/master/lists/hg38-blacklist.v2.bed.gz"
MIN_SIZE=190

mkdir -p "${WORK_DIR}"

echo "Fetching data..."
wget -nc -q --show-progress "${DATA_URL}" -O "${WORK_DIR}/raw.tgz"
tar -xzf "${WORK_DIR}/raw.tgz" -C "${WORK_DIR}"

wget -nc -q --show-progress "${BLACKLIST_URL}" -O "${WORK_DIR}/hg38-blacklist.v2.bed.gz"
gunzip -f "${WORK_DIR}/hg38-blacklist.v2.bed.gz"

echo "Splitting..."
uv run scripts/rme_split.py "${WORK_DIR}/*_mnemonics.bed.gz" "${WORK_DIR}/"

echo "Processing pipeline..."
find "${WORK_DIR}" -maxdepth 1 -name "*.bed" ! -name "hg38-blacklist.v2.bed" | parallel --bar \
    "bedtools intersect -v -a {} -b ${WORK_DIR}/hg38-blacklist.v2.bed | \
    awk '(\$3-\$2) >= ${MIN_SIZE}' | \
    bedtools sort -i - | \
    bgzip > {.}.bed.gz && rm {}"

echo "Cleaning up..."
rm -f "${WORK_DIR}/raw.tgz" "${WORK_DIR}/"*_mnemonics.bed.gz "${WORK_DIR}/hg38-blacklist.v2.bed"

echo "Total files: $(find "${WORK_DIR}" -maxdepth 1 -name "*.bed.gz" | wc -l)"
echo "Total intervals: $(zcat "${WORK_DIR}"/*.bed.gz | wc -l)"
# Total files: 1905
# Total intervals: 55563321
```

**For hg19** compatible, avoid the liftover:

```bash
# all_hg38lift.mnemonics.bedFiles.tgz  -->  all.mnemonics.bedFiles.tgz
wget http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all.mnemonics.bedFiles.tgz -O raw.tgz
```
