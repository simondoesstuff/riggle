"""Split each EDACC ChromHMM source file into per-(cellType, chromatinState) beds.

Snakemake script: protocol — accesses snakemake.input[0] (raw dir),
snakemake.output[0] (split dir), snakemake.threads (workers).
"""

import glob
import gzip
import multiprocessing
import os
from contextlib import ExitStack
from functools import partial
from pathlib import Path
from typing import TextIO, cast

states_map = {
    "1_TssA": "Active_TSS",
    "2_TssAFlnk": "Flanking_Active_TSS",
    "3_TxFlnk": "Transcr_at_gene_5_and_3",
    "4_Tx": "Strong_transcription",
    "5_TxWk": "Weak_transcription",
    "6_EnhG": "Genic_enhancers",
    "7_Enh": "Enhancers",
    "8_ZNF/Rpts": "ZNF_genes_and_repeats",
    "9_Het": "Heterochromatin",
    "10_TssBiv": "Bivalent_Poised_TSS",
    "11_BivFlnk": "Flanking_Bivalent_TSS_Enh",
    "12_EnhBiv": "Bivalent_Enhancer",
    "13_ReprPC": "Repressed_PolyComb",
    "14_ReprPCWk": "Weak_Repressed_PolyComb",
    "15_Quies": "Quiescent_Low",
}

edacc_map = {
    "E001": "ES_I3_Cell_Line",
    "E002": "ES_WA7_Cell_Line",
    "E003": "H1_Cell_Line",
    "E004": "H1_BMP4_Derived_Mesendoderm_Cultured_Cells",
    "E005": "H1_BMP4_Derived_Trophoblast_Cultured_Cells",
    "E006": "H1_Derived_Mesenchymal_Stem_Cells",
    "E007": "H1_Derived_Neuronal_Progenitor_Cultured_Cells",
    "E008": "H9_Cell_Line",
    "E009": "H9_Derived_Neuronal_Progenitor_Cultured_Cells",
    "E010": "H9_Derived_Neuron_Cultured_Cells",
    "E011": "hESC_Derived_CD184pp_Endoderm_Cultured_Cells",
    "E012": "hESC_Derived_CD56pp_Ectoderm_Cultured_Cells",
    "E013": "hESC_Derived_CD56pp_Mesoderm_Cultured_Cells",
    "E014": "HUES48_Cell_Line",
    "E015": "HUES6_Cell_Line",
    "E016": "HUES64_Cell_Line",
    "E017": "IMR90_Cell_Line",
    "E018": "iPS_15b_Cell_Line",
    "E019": "iPS_18_Cell_Line",
    "E020": "iPS_20b_Cell_Line",
    "E021": "iPS_DF_6.9_Cell_Line",
    "E022": "iPS_DF_19.11_Cell_Line",
    "E023": "Mesenchymal_Stem_Cell_Derived_Adipocyte_Cultured_Cells",
    "E024": "4star",
    "E025": "Adipose_Derived_Mesenchymal_Stem_Cell_Cultured_Cells",
    "E026": "Bone_Marrow_Derived_Mesenchymal_Stem_Cell_Cultured_Cells",
    "E027": "Breast_Myoepithelial_Cells",
    "E028": "Breast_vHMEC",
    "E029": "CD14_Primary_Cells",
    "E030": "CD15_Primary_Cells",
    "E031": "CD19_Primary_Cells_Cord_BI",
    "E032": "CD19_Primary_Cells_Peripheral_UW",
    "E033": "CD3_Primary_Cells_Cord_BI",
    "E034": "CD3_Primary_Cells_Peripheral_UW",
    "E035": "CD34_Primary_Cells",
    "E036": "CD34_Cultured_Cells",
    "E037": "CD4_Memory_Primary_Cells",
    "E038": "CD4_Naive_Primary_Cells",
    "E039": "CD4pp_CD25__CD45RApp_Naive_Primary_Cells",
    "E040": "CD4pp_CD25__CD45ROpp_Memory_Primary_Cells",
    "E041": "CD4pp_CD25__IL17__PMA_Ionomycin_stimulated_MACS_purified_Th_Primary_Cells",
    "E042": "CD4pp_CD25__IL17pp_PMA_Ionomcyin_stimulated_Th17_Primary_Cells",
    "E043": "CD4pp_CD25__Th_Primary_Cells",
    "E044": "CD4pp_CD25pp_CD127__Treg_Primary_Cells",
    "E045": "CD4pp_CD25int_CD127pp_Tmem_Primary_Cells",
    "E046": "CD56_Primary_Cells",
    "E047": "CD8_Naive_Primary_Cells",
    "E048": "CD8_Memory_Primary_Cells",
    "E050": "Mobilized_CD34_Primary_Cells_Female",
    "E051": "Mobilized_CD34_Primary_Cells_Male",
    "E049": "Chondrocytes_from_Bone_Marrow_Derived_Mesenchymal_Stem_Cell_Cultured_Cells",
    "E052": "Muscle_Satellite_Cultured_Cells",
    "E053": "Neurosphere_Cultured_Cells_Cortex_Derived",
    "E054": "Neurosphere_Cultured_Cells_Ganglionic_Eminence_Derived",
    "E055": "Penis_Foreskin_Fibroblast_Primary_Cells_skin01",
    "E056": "Penis_Foreskin_Fibroblast_Primary_Cells_skin02",
    "E057": "Penis_Foreskin_Keratinocyte_Primary_Cells_skin02",
    "E058": "Penis_Foreskin_Keratinocyte_Primary_Cells_skin03",
    "E059": "Penis_Foreskin_Melanocyte_Primary_Cells_skin01",
    "E061": "Penis_Foreskin_Melanocyte_Primary_Cells_skin03",
    "E062": "Peripheral_Blood_Mononuclear_Primary_Cells",
    "E063": "Adipose_Nuclei",
    "E065": "Aorta",
    "E066": "Adult_Liver",
    "E067": "Brain_Angular_Gyrus",
    "E068": "Brain_Anterior_Caudate",
    "E069": "Brain_Cingulate_Gyrus",
    "E070": "Brain_Germinal_Matrix",
    "E071": "Brain_Hippocampus_Middle",
    "E072": "Brain_Inferior_Temporal_Lobe",
    "E073": "Brain_Mid_Frontal_Lobe",
    "E074": "Brain_Substantia_Nigra",
    "E075": "Colonic_Mucosa",
    "E076": "Colon_Smooth_Muscle",
    "E077": "Duodenum_Mucosa",
    "E078": "Duodenum_Smooth_Muscle",
    "E079": "Esophagus",
    "E080": "Fetal_Adrenal_Gland",
    "E081": "Fetal_Brain_Male",
    "E082": "Fetal_Brain_Female",
    "E083": "Fetal_Heart",
    "E084": "Fetal_Intestine_Large",
    "E085": "Fetal_Intestine_Small",
    "E086": "Fetal_Kidney",
    "E087": "Pancreatic_Islets",
    "E088": "Fetal_Lung",
    "E089": "Fetal_Muscle_Trunk",
    "E090": "Fetal_Muscle_Leg",
    "E091": "Fetal_Placenta",
    "E092": "Fetal_Stomach",
    "E093": "Fetal_Thymus",
    "E094": "Gastric",
    "E095": "Left_Ventricle",
    "E096": "Lung",
    "E097": "Ovary",
    "E098": "Pancreas",
    "E099": "Placenta_Amnion",
    "E100": "Psoas_Muscle",
    "E101": "Rectal_Mucosa.Donor_29",
    "E102": "Rectal_Mucosa.Donor_31",
    "E103": "Rectal_Smooth_Muscle",
    "E104": "Right_Atrium",
    "E105": "Right_Ventricle",
    "E106": "Sigmoid_Colon",
    "E107": "Skeletal_Muscle_Male",
    "E108": "Skeletal_Muscle_Female",
    "E109": "Small_Intestine",
    "E110": "Stomach_Mucosa",
    "E111": "Stomach_Smooth_Muscle",
    "E112": "Thymus",
    "E113": "Spleen",
    "E114": "A549_EtOH_0.02pct_Lung_Carcinoma",
    "E115": "Dnd41_TCell_Leukemia",
    "E116": "GM12878_Lymphoblastoid",
    "E117": "HeLa_S3_Cervical_Carcinoma",
    "E118": "HepG2_Hepatocellular_Carcinoma",
    "E119": "HMEC_Mammary_Epithelial",
    "E120": "HSMM_Skeletal_Muscle_Myoblasts",
    "E121": "HSMMtube_Skeletal_Muscle_Myotubes_Derived_from_HSMM",
    "E122": "HUVEC_Umbilical_Vein_Endothelial_Cells",
    "E123": "K562_Leukemia",
    "E124": "Monocytes_CD14pp_RO01746",
    "E125": "NH_A_Astrocytes",
    "E126": "NHDF_Ad_Adult_Dermal_Fibroblasts",
    "E127": "NHEK_Epidermal_Keratinocytes",
    "E128": "NHLF_Lung_Fibroblasts",
    "E129": "Osteoblasts",
}


def _smart_open(path: Path, mode: str = "rt") -> TextIO:
    if path.suffix == ".gz":
        return cast(TextIO, gzip.open(path, mode))
    return cast(TextIO, open(path, mode, encoding="utf-8"))


def _process_file(
    input_file: Path,
    output_dir: Path,
    states_map: dict[str, str],
    edacc_map: dict[str, str],
) -> tuple[int, set[Path]]:
    print(f"Processing: {input_file.name}")
    open_files: dict[str, TextIO] = {}
    output_filepaths: set[Path] = set()
    processed_lines = 0

    eid = input_file.name.split("_")[0]
    sample_name = edacc_map.get(eid)
    if not sample_name:
        raise ValueError(f"No EDACC name for ID '{eid}' in {input_file.name}")

    with ExitStack() as stack, _smart_open(input_file, "rt") as f_in:
        for line_number, line in enumerate(f_in, 1):
            toks = line.rstrip().split("\t")
            if len(toks) < 4:
                raise ValueError(
                    f"Malformed line in {input_file.name} (line {line_number}): "
                    f"expected 4+ columns, got {len(toks)}"
                )
            state_id = toks[3]
            state_name = states_map.get(state_id)
            if not state_name:
                raise ValueError(
                    f"Unknown state '{state_id}' in {input_file.name} (line {line_number})"
                )
            output_filename = f"{sample_name}_{state_name}.bed"
            if output_filename not in open_files:
                out_path = output_dir / output_filename
                f_out = stack.enter_context(open(out_path, "wt", encoding="utf-8"))
                open_files[output_filename] = f_out
                output_filepaths.add(out_path)
            open_files[output_filename].write(line)
            processed_lines += 1

    return processed_lines, output_filepaths


# --- Snakemake script: entry point ---

raw_dir = snakemake.input[0]
output_dir = Path(snakemake.output[0])
num_processes = os.cpu_count()

output_dir.mkdir(parents=True, exist_ok=True)

input_files = [Path(f) for f in glob.glob(f"{raw_dir}/*_mnemonics.bed.gz")]
if not input_files:
    raise SystemExit(f"No *_mnemonics.bed.gz files found in {raw_dir}")

print(f"Found {len(input_files)} files to process with {num_processes} worker(s).")

_process_partial = partial(
    _process_file,
    output_dir=output_dir,
    states_map=states_map,
    edacc_map=edacc_map,
)

ctx = multiprocessing.get_context("fork")
with ctx.Pool(processes=num_processes) as pool:
    results = pool.map(_process_partial, input_files)

total_lines = sum(r for r, _ in results)
total_files = len({p for _, paths in results for p in paths})
print(f"Done: {total_lines} lines → {total_files} output files in {output_dir}")
