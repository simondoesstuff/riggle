# Data

## Default Set

2026/04/16

Targeting human interval files.

Taking specifically:

- hg38 (homo-sapien)
- non-control experiments
- "experiments" as opposed to "datasets"
- not perturbed
- status: released
- conservative & optimal IDR thresholded peaks
- the following file formats:

  - bed narrowPeak
  - bed idr_ranked_peak
  - bed bed3+
  - bed broadPeak
  - bed bed3
  - bed bed9+
  - bed bedRnaElements
  - bed tss_peak
  - bed idr_peak
  - bed bedMethyl
  - bed bed9
  - bed bed6+

https://www.encodeproject.org/search/?type=Experiment&control_type!=*&status=released&perturbed=false&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assembly=GRCh38&files.file_type=bed+narrowPeak&files.file_type=bed+idr_ranked_peak&files.file_type=bed+bed3+&files.file_type=bed+broadPeak&files.file_type=bed+bed3&files.file_type=bed+bed9+&files.file_type=bed+bedRnaElements&files.file_type=bed+tss_peak&files.file_type=bed+idr_peak&files.file_type=bed+bedMethyl&files.file_type=bed+bed6+&files.file_type=bed+bed9&files.output_type=conservative+IDR+thresholded+peaks&files.output_type=optimal+IDR+thresholded+peaks

Then, "Processed" files were downloaded.

### Clean

- filter for normal chromosomes (`chr[0-9XY]+`)
- apply the Boyle Lab hg38 [blacklist](https://github.com/Boyle-Lab/Blacklist)

```sh
data=data/encode

wget https://github.com/Boyle-Lab/Blacklist/raw/refs/heads/master/lists/hg38-blacklist.v2.bed.gz -P $data

if [ -f $data/"ENCODE Search Files.txt" ]; then
    mkdir $data/beds
    tail -n +2 "$data/ENCODE Search Files.txt" | parallel --bar -j 4 wget {} -q -c --retry-connrefused -P "$data/beds"
    parallel --bar "zcat {} | \
        grep -E '^chr[0-9XY]+[[:blank:]]' | \
        bedtools intersect -v -a - -b $data/hg38-blacklist.v2.bed.gz | \
        gzip > {}.tmp && \
        mv {}.tmp {}" ::: "$data/beds/"*
    find "$data/beds/" -maxdepth 1 -name "*.gz" | parallel "zcat {} | wc -l" | awk '{total += $1} END {print total}'
    # 169967213
    find "$data/beds/" -maxdepth 1 -name "*.gz" | wc -l
    # 9456
else
    echo "Missing 'ENCODE Search Files.txt'"
fi
```
