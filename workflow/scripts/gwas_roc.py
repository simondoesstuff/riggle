"""Compute ROC-AUC benchmark for GWAS disease-tissue enrichment tools.

Scores are aggregated from bed-file level to cell-type level before ROC
computation:
  1. Parse each bed file into (cell_type, chromatin_state, score).
  2. Z-score normalise within each chromatin state column — uninformative
     states collapse to zero, informative ones stay contrasted.
  3. Average z-scores across states → one score per cell type.
  4. Label each cell type as TP/FP against the disease's ground-truth tissues.
  5. Compute ROC on the cell-type level.

Usage:
    uv run workflow/scripts/gwas_roc.py
        --disease-tissue data/gwas/disease_tissue.json
        --giggle-dir data/giggle
        --chuckle-dir data/chuckle
        --output data/gwas/roc_data.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from utils import parse_bed_filename, read_chuckle_records, read_tsv

# ---------------------------------------------------------------------------
# Bed-file → tissue classification
# ---------------------------------------------------------------------------

# Cell type prefix (uppercased) keyword → tissue categories.
# All matching rules are applied (multi-label).
# 4star_*, ESC, iPSC, Breast/Mammary, HeLa: no matching category → excluded.
_CELL_TISSUE_RULES: list[tuple[list[str], list[str]]] = [
    # Nervous
    (["BRAIN", "ASTROCYTE", "NEUROSPHERE", "NEURONAL_PROGENITOR",
      "H9_DERIVED_NEURON", "FETAL_BRAIN"], ["nervous"]),
    # Immune (T-cells, B-cells, NK-cells, monocytes, lymphoid lines)
    (["CD3", "CD4", "CD8", "CD14", "CD15", "CD56", "GM12878", "LYMPHOBLASTOID",
      "PERIPHERAL_BLOOD", "MONOCYTE", "DND41", "T_CELL"], ["immune"]),
    # Thymus and Spleen: immune + hematopoietic
    (["THYMUS", "SPLEEN"], ["immune", "hematopoietic"]),
    # Hematopoietic stem / progenitor
    (["CD19", "CD34", "K562"], ["immune", "hematopoietic"]),
    # Cardiovascular
    (["VENTRICLE", "ATRIUM", "AORTA", "HUVEC", "FETAL_HEART"], ["cardiovascular"]),
    # Respiratory
    (["LUNG", "NHLF"], ["respiratory"]),
    # Gastrointestinal
    (["COLON", "INTESTINE", "RECTAL", "DUODENUM", "ESOPHAGUS", "GASTRIC",
      "STOMACH", "MUCOSA", "SIGMOID"], ["gastrointestinal"]),
    # Liver
    (["LIVER", "HEPG2", "HEPATOCELLULAR"], ["liver"]),
    # Renal
    (["KIDNEY"], ["renal"]),
    # Endocrine + metabolic (pancreas, adrenal)
    (["PANCREAS", "PANCREATIC_ISLET", "ADRENAL"], ["endocrine", "metabolic"]),
    # Musculoskeletal (muscle, bone, cartilage)
    (["SKELETAL_MUSCLE", "HSMM", "PSOAS", "FETAL_MUSCLE",
      "OSTEOBLAST", "CHONDROCYTE", "MUSCLE_SATELLITE"], ["musculoskeletal"]),
    # Connective tissue / mesenchymal
    (["MESENCHYMAL", "IMR90"], ["connective_tissue"]),
    # Adipose → metabolic + connective_tissue
    (["ADIPOSE", "ADIPOCYTE"], ["metabolic", "connective_tissue"]),
    # Skin (keratinocytes, melanocytes, foreskin-derived)
    (["KERATINOCYTE", "MELANOCYTE", "NHEK", "EPIDERMAL",
      "FORESKIN", "NHDF"], ["skin", "connective_tissue"]),
]

def _prefix_to_tissues(prefix: str) -> list[str]:
    upper = prefix.upper()
    tissues: set[str] = set()
    for keywords, cats in _CELL_TISSUE_RULES:
        if any(kw in upper for kw in keywords):
            tissues.update(cats)
    return sorted(tissues)


def bed_to_tissues(bed_name: str) -> list[str]:
    parsed = parse_bed_filename(bed_name)
    return _prefix_to_tissues(parsed[0]) if parsed else []


# ---------------------------------------------------------------------------
# Score loaders
# ---------------------------------------------------------------------------


def load_giggle(path: Path) -> dict[str, float]:
    """Parse giggle TSV → {bed_filename: combo_score}."""
    header, rows = read_tsv(path)
    for col in ("file", "combo_score"):
        if col not in header:
            raise ValueError(f"Missing column '{col}' in {path}")
    scores: dict[str, float] = {}
    for r in rows:
        try:
            scores[Path(r["file"]).name] = float(r["combo_score"])
        except (KeyError, ValueError):
            continue
    return scores


def load_chuckle(path: Path) -> dict[str, float]:
    """Parse chuckle JSON → {bed_filename: llr}."""
    scores: dict[str, float] = {}
    for r in read_chuckle_records(path):
        if r.get("llr") is None:
            continue
        scores[Path(r["db_name"]).name] = float(r["llr"])
    return scores


# ---------------------------------------------------------------------------
# Aggregation: bed-file scores → per-cell-type scores
# ---------------------------------------------------------------------------


def aggregate_to_cell_types(scores: dict[str, float]) -> dict[str, float]:
    """Z-score within each chromatin state, then average across states.

    Returns {cell_type_prefix: aggregated_score}.

    Rationale: uninformative chromatin states have near-zero variance across
    cell types → their z-scores collapse to ≈0 and contribute nothing to the
    per-cell-type average.  States with genuine tissue-specific signal keep
    their contrast.
    """
    # Parse all bed names into (prefix, state, score)
    triples: list[tuple[str, str, float]] = []
    for bed, score in scores.items():
        parsed = parse_bed_filename(bed)
        if parsed is not None:
            triples.append((*parsed, score))

    if not triples:
        return {}

    cell_types = sorted({t[0] for t in triples})
    states = sorted({t[1] for t in triples})
    ct_idx = {ct: i for i, ct in enumerate(cell_types)}
    st_idx = {st: i for i, st in enumerate(states)}

    mat = np.full((len(cell_types), len(states)), np.nan)
    for prefix, state, score in triples:
        mat[ct_idx[prefix], st_idx[state]] = score

    # Z-score each column (chromatin state) independently
    col_mean = np.nanmean(mat, axis=0)
    col_std = np.nanstd(mat, axis=0)
    col_std[col_std == 0] = 1.0   # flat columns → no signal → divide by 1
    normed = (mat - col_mean) / col_std
    normed = np.nan_to_num(normed, nan=0.0)

    # Average z-scores across all states for each cell type
    agg = normed.mean(axis=1)

    return {ct: float(agg[ct_idx[ct]]) for ct in cell_types}


# ---------------------------------------------------------------------------
# ROC computation
# ---------------------------------------------------------------------------


def _roc_curve(
    labeled: list[tuple[float, int]],
) -> tuple[list[float], list[float], float] | None:
    n_pos = sum(l for _, l in labeled)
    n_neg = len(labeled) - n_pos
    if n_pos == 0 or n_neg == 0:
        return None

    labeled_sorted = sorted(labeled, key=lambda x: -x[0])
    fpr_pts = [0.0]
    tpr_pts = [0.0]
    tp = fp = 0
    for _, label in labeled_sorted:
        if label:
            tp += 1
        else:
            fp += 1
        fpr_pts.append(fp / n_neg)
        tpr_pts.append(tp / n_pos)

    auc = float(np.trapezoid(tpr_pts, fpr_pts))
    return fpr_pts, tpr_pts, auc


def _roc_auc(
    scores: dict[str, float],
    true_tissues: list[str],
) -> tuple[list[float], list[float], float, int, int, list[tuple[float, int]]] | None:
    """Aggregate to cell-type level, then compute ROC."""
    cell_scores = aggregate_to_cell_types(scores)
    labeled: list[tuple[float, int]] = []
    for prefix, score in cell_scores.items():
        tissues = _prefix_to_tissues(prefix)
        if not tissues:
            continue
        labeled.append((score, 1 if set(tissues) & set(true_tissues) else 0))

    if not labeled:
        return None

    result = _roc_curve(labeled)
    if result is None:
        return None

    fpr, tpr, auc = result
    n_pos = sum(l for _, l in labeled)
    n_neg = len(labeled) - n_pos
    return fpr, tpr, auc, n_pos, n_neg, labeled


def _combined_roc(
    per_disease: dict[str, dict],
) -> tuple[list[float], list[float], float]:
    """Pool rank-normalised per-cell-type scores across diseases."""
    pooled: list[tuple[float, int]] = []
    for ddata in per_disease.values():
        pairs: list[tuple[float, int]] = ddata.get("ranked_pairs", [])
        n = len(pairs)
        if n == 0:
            continue
        for rank, (_, label) in enumerate(sorted(pairs, key=lambda x: -x[0])):
            pooled.append((1.0 - rank / n, label))

    if not pooled:
        return [0.0, 1.0], [0.0, 1.0], 0.5

    result = _roc_curve(pooled)
    if result is None:
        return [0.0, 1.0], [0.0, 1.0], 0.5
    return result


# ---------------------------------------------------------------------------
# Per-method pipeline
# ---------------------------------------------------------------------------


def compute_method(
    disease_tissue: dict[str, dict],
    result_files: dict[str, Path],
    loader,
) -> dict:
    per_disease: dict[str, dict] = {}

    for disease, info in disease_tissue.items():
        if disease not in result_files:
            continue
        true_tissues: list[str] = info.get("tissues", [])
        if not true_tissues:
            continue
        try:
            scores = loader(result_files[disease])
        except Exception as e:
            print(f"  WARNING: failed to load {result_files[disease]}: {e}", file=sys.stderr)
            continue

        result = _roc_auc(scores, true_tissues)
        if result is None:
            print(f"  WARNING: no classifiable cell types for {disease}", file=sys.stderr)
            continue

        fpr, tpr, auc, n_pos, n_neg, labeled = result
        per_disease[disease] = {
            "fpr": fpr,
            "tpr": tpr,
            "auc": auc,
            "n_pos": n_pos,
            "n_neg": n_neg,
            "true_tissues": true_tissues,
            "ranked_pairs": labeled,
        }
        print(f"  {disease}: AUC={auc:.3f}  n_pos={n_pos}  n_neg={n_neg}", file=sys.stderr)

    comb_fpr, comb_tpr, comb_auc = _combined_roc(per_disease)
    print(f"  combined AUC = {comb_auc:.3f}", file=sys.stderr)

    for d in per_disease.values():
        del d["ranked_pairs"]

    return {
        "diseases": per_disease,
        "combined": {"fpr": comb_fpr, "tpr": comb_tpr, "auc": comb_auc},
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--disease-tissue", required=True)
    ap.add_argument("--giggle-dir")
    ap.add_argument("--chuckle-dir")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    with open(args.disease_tissue) as f:
        disease_tissue: dict[str, dict] = json.load(f)

    output: dict[str, dict] = {"methods": {}}

    if args.giggle_dir:
        giggle_dir = Path(args.giggle_dir)
        giggle_files = {p.stem: p for p in giggle_dir.glob("*.tsv") if p.stem in disease_tissue}
        print(f"Giggle: {len(giggle_files)} diseases", file=sys.stderr)
        output["methods"]["giggle"] = compute_method(disease_tissue, giggle_files, load_giggle)

    if args.chuckle_dir:
        chuckle_dir = Path(args.chuckle_dir)
        chuckle_files = {p.stem: p for p in chuckle_dir.glob("*.json") if p.stem in disease_tissue}
        print(f"Chuckle: {len(chuckle_files)} diseases", file=sys.stderr)
        if chuckle_files:
            output["methods"]["chuckle"] = compute_method(disease_tissue, chuckle_files, load_chuckle)
        else:
            print("  (no chuckle results found)", file=sys.stderr)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(output, indent=2))
    print(f"Wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
