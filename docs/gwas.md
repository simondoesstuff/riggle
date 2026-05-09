https://www.nature.com/articles/nature13835

some dependencies are required

```bash
nix shell nixpkgs#csvkit nixpkgs#qsv
```

and you need liftover, see the [flake.nix](flake.nix)

```bash
data=data/gwas

mkdir "$data"
wget https://www.nature.com/nature/journal/v518/n7539/extref/nature13835-s1.xls -O "$data/gwas.xls"
# convert to csv, subset cols, remove integer .0 decimals
in2csv "$data/gwas.xls" | csvcut -c 1,4,5 | sed 's/\.0$//g' > "$data/gwas.csv"
# split, grouped by col 1 "Disease"
qsv partition Disease "$data" "$data/gwas.csv"
# remove header, remove col 1, add new col ($col3 + 1), convert to tsv (bed)
parallel 'awk -F, "NR>1 {OFS=\"\t\"; print \$2, \$3, \$3+1}" {} > {.}.bed' ::: "$data"/*.csv
# liftover
wget -P $data https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
parallel --bar "liftover {} \"$data/hg19ToHg38.over.chain.gz\" {}_ /dev/null 2> /dev/null && mv {}_ {}" ::: "$data"/*.bed

# cleanup
rm "$data"/*.csv "$data"/gwas.bed "$data"/gwas.xls "$data"/hg19ToHg38.over.chain.gz



wc -l "$data"/*.bed

   160 data/gwas/Allergy.bed
    66 data/gwas/Alopecia_areata.bed
   175 data/gwas/Alzheimers_combined.bed
   207 data/gwas/Ankylosing_spondylitis.bed
   121 data/gwas/Asthma.bed
   117 data/gwas/Atopic_dermatitis.bed
    48 data/gwas/Autoimmune_thyroiditis.bed
    66 data/gwas/Behcets_disease.bed
   430 data/gwas/Bone_mineral_density.bed
   141 data/gwas/C_reactive_protein.bed
   370 data/gwas/Celiac_disease.bed
   197 data/gwas/Chronic_kidney_disease.bed
    47 data/gwas/Creatinine_levels.bed
   787 data/gwas/Crohns_disease.bed
   110 data/gwas/Fasting_glucose_related_traits.bed
   329 data/gwas/HDL_cholesterol.bed
   236 data/gwas/Juvenile_idiopathic_arthritis.bed
    56 data/gwas/Kawasaki_disease.bed
   261 data/gwas/LDL_cholesterol.bed
   173 data/gwas/Liver_enzyme_levels_gamma_glutamyl_transferase.bed
    63 data/gwas/Migraine.bed
   696 data/gwas/Multiple_sclerosis.bed
   368 data/gwas/Platelet_counts.bed
   275 data/gwas/Primary_biliary_cirrhosis.bed
    80 data/gwas/Primary_sclerosing_cholangitis.bed
    74 data/gwas/Progressive_supranuclear_palsy.bed
   264 data/gwas/Psoriasis.bed
   298 data/gwas/Red_blood_cell_traits.bed
   128 data/gwas/Renal_function_related_traits_BUN.bed
    46 data/gwas/Restless_legs_syndrome.bed
   195 data/gwas/Rheumatoid_arthritis.bed
   221 data/gwas/Systemic_lupus_erythematosus.bed
    43 data/gwas/Systemic_sclerosis.bed
   234 data/gwas/Triglycerides.bed
   339 data/gwas/Type_1_diabetes.bed
   522 data/gwas/Type_2_diabetes.bed
   447 data/gwas/Ulcerative_colitis.bed
   192 data/gwas/Urate_levels.bed
   155 data/gwas/Vitiligo.bed
  8737 total
```
