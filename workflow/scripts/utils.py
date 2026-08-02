"""Shared utilities for workflow/scripts."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

FLOAT_TINY: float = float(np.finfo(np.float64).tiny)

CHROMATIN_STATES: list[str] = [
    "Active_TSS",
    "Flanking_Active_TSS",
    "Strong_transcription",
    "Weak_transcription",
    "Enhancers",
    "Genic_enhancers",
    "ZNF_genes_and_repeats",
    "Heterochromatin",
    "Bivalent_Poised_TSS",
    "Flanking_Bivalent_TSS_Enh",
    "Bivalent_Enhancer",
    "Repressed_PolyComb",
    "Weak_Repressed_PolyComb",
    "Transcr_at_gene_5_and_3",
    "Quiescent_Low",
]


def strip_bed_name(filename: str) -> str:
    """Strip common BED file extensions from a filename."""
    for ext in (".clean.bed.gz", ".bed.gz", ".clean.bed", ".bed"):
        if filename.endswith(ext):
            return filename[: -len(ext)]
    return filename


def parse_bed_filename(name: str) -> tuple[str, str] | None:
    """Parse a RME BED filename into (cell_type_prefix, chromatin_state).

    Returns None if the name doesn't end with a known chromatin state.
    """
    stem = strip_bed_name(Path(name).name)
    for state in sorted(CHROMATIN_STATES, key=len, reverse=True):
        if stem.endswith(state):
            return stem[: -len(state)].rstrip("_"), state
    return None


def read_chuckle_records(path) -> list[dict]:
    """Load a chuckle JSON result and normalize to a list of per-bed records."""
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, list):
        return data
    return next(v for v in data.values() if isinstance(v, list))


def read_tsv(path) -> tuple[list[str], list[dict[str, str]]]:
    """Parse a TSV file (with optional leading '#' on the header) into (header, rows).

    Each row is a dict keyed by column name.  Rows with too few fields are skipped.
    """
    with open(path) as f:
        lines = [line.rstrip("\n\t") for line in f if line.strip()]
    if not lines:
        return [], []
    header = lines[0].lstrip("#").split("\t")
    n = len(header)
    rows = []
    for line in lines[1:]:
        parts = line.lstrip("#").split("\t")
        if len(parts) >= n:
            rows.append(dict(zip(header, parts[:n])))
    return header, rows
