"""Map GWAS trait names to macro-level tissue/system categories via DOID.

For each .bed file in the GWAS hg19 directory, finds the best-matching DOID
term via is_a hierarchy walk.  Non-disease quantitative traits (cholesterol
levels, platelet counts, etc.) fall back to keyword-based rules.

Usage:
    uv run workflow/scripts/gwas_disease_tissue.py doid.obo gwas_hg19_dir output.json
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import unicodedata
from pathlib import Path

import pronto


# Maps tissue label → the DOID preferred term name that anchors it.
# A disease inherits every tissue category whose root is one of its ancestors.
TISSUE_ROOTS: dict[str, str] = {
    "immune": "immune system disease",
    "nervous": "nervous system disease",
    "musculoskeletal": "musculoskeletal system disease",
    "cardiovascular": "cardiovascular system disease",
    "hematopoietic": "hematopoietic system disease",
    "skin": "skin disease",
    "gastrointestinal": "gastrointestinal system disease",
    "liver": "liver disease",
    "renal": "urinary system disease",
    "endocrine": "endocrine system disease",
    "metabolic": "disease of metabolism",
    "respiratory": "respiratory system disease",
    "connective_tissue": "connective tissue disease",
}

# Definition-text keywords: scans the matched DOID term's own definition for
# tissue associations that the is_a hierarchy doesn't capture.  DOID's tree
# places IBD under GI (not immune) and T2DM under metabolic (not endocrine),
# but the definitions themselves encode the missing biology.
# Note: Crohn's disease and PSC lack ANY immune signal in DOID (hierarchy,
# def text, subsets, or relationships) so those cannot be recovered here.
DEF_TISSUE_KEYWORDS: list[tuple[list[str], list[str]]] = [
    (["autoimmun", "autoantibodie"], ["immune"]),
    (["insulin", "beta cell", "islet of langerhans"], ["endocrine"]),
    (["vasculit", "blood vessel", "arterit"], ["cardiovascular"]),
]

# Keyword rules for quantitative traits that have no DOID disease match.
# Each entry: (list[keywords_any_of], list[tissues]).
KEYWORD_MAP: list[tuple[list[str], list[str]]] = [
    (["cholesterol", "triglyceride", "hdl", "ldl", "lipid"], ["cardiovascular", "metabolic"]),
    (["glucose", "fasting glucose", "insulin"], ["endocrine", "metabolic"]),
    (["creatinine", "renal function", "blood urea nitrogen", "bun"], ["renal", "metabolic"]),
    (["urate", "uric acid"], ["metabolic", "renal"]),
    (["platelet"], ["hematopoietic"]),
    (["red blood cell", "erythrocyte", "hemoglobin"], ["hematopoietic"]),
    (["bone mineral", "bone density"], ["musculoskeletal", "metabolic"]),
    (["liver enzyme", "gamma glutamyl", "ggt", "alanine aminotransferase", "alt"], ["liver", "metabolic"]),
    (["c reactive protein", "crp"], ["immune"]),
]


def _norm(name: str) -> str:
    """Lowercase, remove accents/apostrophes, collapse whitespace."""
    name = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode()
    name = name.lower().replace("_", " ").replace("'", "")
    return re.sub(r"\s+", " ", name).strip()


def _word_set(s: str) -> frozenset[str]:
    """Words of a normalised string, with trailing 's' stripped for fuzzy matching."""
    return frozenset(w.rstrip("s") for w in s.split())


def _candidates(raw: str) -> list[str]:
    """Normalized forms to try for a trait name, from most to least specific."""
    base = _norm(raw)
    result = [base]
    for suffix_re in (
        r"\s+combined$",
        r"\s+related\s+traits?$",
        r"\s+levels?$",
        r"\s+counts?$",
        r"\s+traits?$",
        r"\s+related$",
        r"\s+measurement$",
    ):
        stripped = re.sub(suffix_re, "", base).strip()
        if stripped and stripped != base and stripped not in result:
            result.append(stripped)
    return result


def build_index(
    ont: pronto.Ontology,
) -> tuple[dict[str, pronto.Term], dict[str, pronto.Term]]:
    name_idx: dict[str, pronto.Term] = {}
    syn_idx: dict[str, pronto.Term] = {}
    for term in ont.terms():
        if not term.name:
            continue
        name_idx[_norm(term.name)] = term
        for syn in term.synonyms:
            syn_idx[_norm(syn.description)] = term
    return name_idx, syn_idx


def resolve_roots(
    name_idx: dict[str, pronto.Term],
) -> dict[str, pronto.Term]:
    roots: dict[str, pronto.Term] = {}
    for label, doid_name in TISSUE_ROOTS.items():
        key = _norm(doid_name)
        if key in name_idx:
            roots[label] = name_idx[key]
        else:
            print(f"WARNING: tissue root '{doid_name}' not found in DOID", file=sys.stderr)
    return roots


def _def_tissues(term: pronto.Term) -> list[str]:
    """Additional tissues inferred from the term's definition text."""
    if not term.definition:
        return []
    text = term.definition.lower()
    result: list[str] = []
    for keywords, labels in DEF_TISSUE_KEYWORDS:
        if any(kw in text for kw in keywords):
            result.extend(labels)
    return result


def tissues_for_term(
    term: pronto.Term,
    roots: dict[str, pronto.Term],
) -> list[str]:
    ancestor_ids = {t.id for t in term.superclasses()}
    from_hierarchy = [label for label, root in roots.items() if root.id in ancestor_ids]
    from_def = _def_tissues(term)
    return sorted(set(from_hierarchy + from_def))


def keyword_tissues(trait: str) -> list[str]:
    norm = _norm(trait)
    for keywords, labels in KEYWORD_MAP:
        if any(kw in norm for kw in keywords):
            return labels
    return []


def match_trait(
    trait: str,
    name_idx: dict[str, pronto.Term],
    syn_idx: dict[str, pronto.Term],
) -> tuple[pronto.Term | None, str]:
    """Try progressively fuzzier lookups; return (term, match_type) or (None, '')."""
    candidates = _candidates(trait)

    # 1. Exact name or synonym match
    for cand in candidates:
        if cand in name_idx:
            return name_idx[cand], "exact_name"
        if cand in syn_idx:
            return syn_idx[cand], "exact_synonym"

    # 2. If a keyword rule matches, this is a quantitative trait — skip fuzzy DOID matching
    #    to avoid spurious hits (e.g. "triglyceride" matching a Wolman disease synonym).
    if keyword_tissues(trait):
        return None, ""

    # 3. Word-subset match: all de-pluralised words in candidate appear in a DOID name/synonym.
    #    Prefer the shortest matching DOID name (most general term).
    #    Require candidate is at least 5 chars to avoid trivial single-letter hits.
    for cand in candidates:
        if len(cand) < 5:
            continue
        cand_words = _word_set(cand)
        hits = [
            (len(nm), term)
            for nm, term in name_idx.items()
            if cand_words <= _word_set(nm)
        ]
        if hits:
            return min(hits, key=lambda x: x[0])[1], f"word_subset:{cand}"
        syn_hits = [
            (len(nm), term)
            for nm, term in syn_idx.items()
            if cand_words <= _word_set(nm)
        ]
        if syn_hits:
            return min(syn_hits, key=lambda x: x[0])[1], f"word_subset_syn:{cand}"

    return None, ""


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("doid_obo", help="Path to doid.obo")
    ap.add_argument("gwas_hg19_dir", help="Directory of per-trait hg19 BED files")
    ap.add_argument("output", help="Output JSON path")
    args = ap.parse_args()

    print("Loading DOID...", file=sys.stderr)
    ont = pronto.Ontology(args.doid_obo)
    name_idx, syn_idx = build_index(ont)
    roots = resolve_roots(name_idx)
    print(
        f"  {len(name_idx)} terms, {len(syn_idx)} synonyms, {len(roots)} tissue roots",
        file=sys.stderr,
    )

    beds = sorted(Path(args.gwas_hg19_dir).glob("*.bed"))
    if not beds:
        sys.exit(f"No .bed files found in {args.gwas_hg19_dir}")

    mapping: dict[str, dict] = {}
    unmatched: list[str] = []

    for bed in beds:
        trait = bed.stem
        term, match_type = match_trait(trait, name_idx, syn_idx)
        if term is not None:
            tis = tissues_for_term(term, roots)
            mapping[trait] = {
                "doid_id": term.id,
                "doid_name": term.name,
                "match": match_type,
                "tissues": tis,
            }
        else:
            tis = keyword_tissues(trait)
            mapping[trait] = {
                "doid_id": None,
                "doid_name": None,
                "match": "keyword" if tis else "none",
                "tissues": tis,
            }
            if not tis:
                unmatched.append(trait)

        flag = "" if mapping[trait]["tissues"] else " [UNMATCHED]"
        doid_label = f" ({term.name})" if term else ""
        print(f"  {trait}: {mapping[trait]['tissues']}{doid_label}{flag}", file=sys.stderr)

    Path(args.output).write_text(json.dumps(mapping, indent=2))

    n_doid = sum(1 for v in mapping.values() if v["doid_id"])
    n_kw = sum(1 for v in mapping.values() if v["match"] == "keyword")
    print(
        f"\nDOID: {n_doid}  keyword: {n_kw}  unmatched: {len(unmatched)}/{len(mapping)}",
        file=sys.stderr,
    )
    if unmatched:
        print("Unmatched traits:", ", ".join(unmatched), file=sys.stderr)
    print(f"Wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
