"""Standalone RME heatmap visualization — condensed from giggle-ml."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.colors import Colormap, Normalize
from matplotlib.figure import Figure
from numpy.typing import NDArray

type Pathish = os.PathLike[str] | str

# ==========================================
# Genomics Domain Specifics
# ==========================================

chromatin_states = [
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

cell_categories = [
    "iPSC",
    "Thymus",
    "Sm muscle",
    "Other",
    "Neurosph",
    "Muscle",
    "Mesench",
    "Lung",
    "Heart",
    "HSC and B cell",
    "Epithelial",
    "ESC",
    "ES deriv",
    "Digestive",
    "Cancer cell line",
    "Brain",
    "Blood and T-cell",
]

category_keywords_map = {
    "iPSC": ["IPS"],
    "Thymus": ["THYMUS", "SPLEEN"],
    "Sm muscle": ["SMOOTH_MUSCLE"],
    "Neurosph": ["NEUROSPHERE"],
    "Muscle": ["SKELETAL", "HSMM", "PSOAS", "FETAL_MUSCLE"],
    "Mesench": [
        "MESENCHYMAL",
        "FIBROBLAST",
        "ADIPOSE",
        "OSTEOBLAST",
        "CHONDROCYTE",
        "IMR90",
    ],
    "Lung": ["LUNG", "NHLF"],
    "Heart": ["HEART", "AORTA", "VENTRICLE", "ATRIUM", "HUVEC"],
    "HSC and B cell": ["CD19", "CD34", "GM12878"],
    "Epithelial": ["EPITHELIAL", "HMEC", "BREAST"],
    "ESC": ["H1_CELL_LINE", "H9_CELL_LINE", "HUES", "ES_WA7", "ES_I3"],
    "ES deriv": ["H1_BMP4_DERIVED", "HESC_DERIVED", "H1_DERIVED", "H9_DERIVED"],
    "Digestive": [
        "LIVER",
        "STOMACH",
        "INTESTINE",
        "COLON",
        "RECTAL",
        "DUODENUM",
        "ESOPHAGUS",
        "GASTRIC",
        "MUCOSA",
    ],
    "Cancer cell line": [
        "CARCINOMA",
        "LEUKEMIA",
        "K562",
        "HEPG2",
        "HELA",
        "A549",
        "DND41",
    ],
    "Brain": ["BRAIN", "ASTROCYTE"],
    "Blood and T-cell": [
        "CD3",
        "CD4",
        "CD8",
        "CD14",
        "CD15",
        "CD56",
        "PERIPHERAL_BLOOD",
        "MONOCYTE",
    ],
    "Other": [],
}


def classify_cell_type(cell_type: str) -> str:
    cell_type_upper = cell_type.upper()
    for category, keywords in category_keywords_map.items():
        if any(kw in cell_type_upper for kw in keywords):
            return category
    return "Other"


def classify_bed_file(bed: Pathish) -> tuple[str, str]:
    filename = os.path.basename(str(bed))
    if filename.endswith(".bed.gz"):
        filename = filename[:-7]
    elif filename.endswith(".bed"):
        filename = filename[:-4]
    for state in sorted(chromatin_states, key=len, reverse=True):
        if filename.endswith(state):
            prefix = filename[: -(len(state))].rstrip("_")
            return prefix, state
    raise ValueError(f'Unknown format, expecting "cellType_chrmState", got: {filename}')


def parse_score_file(score_path: Pathish) -> tuple[list[str], list[float]]:
    def parse() -> Iterator[tuple[str, float]]:
        with open(score_path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                split = line.strip().split()
                yield Path(split[0]).name.split(".")[0], float(split[-1])

    parsed_tuples = sorted(parse(), key=lambda x: -x[1])
    if not parsed_tuples:
        return [], []
    keys, values = zip(*parsed_tuples)
    return list(keys), list(values)


def parse_json_file(
    json_path: Pathish, score_field: str = "p_value"
) -> tuple[list[str], list[float]]:
    with open(json_path) as f:
        data = json.load(f)

    # Accept a bare list or find the first list value in an object
    records: list[dict] = (
        data
        if isinstance(data, list)
        else next(v for v in data.values() if isinstance(v, list))
    )

    # Find the file-path-like field with the most unique values (skip score_field)
    first = records[0]
    path_fields = [
        k
        for k, v in first.items()
        if isinstance(v, str) and ("." in v or "/" in v) and k != score_field
    ]
    name_field = max(path_fields, key=lambda k: len({r[k] for r in records}))

    def _strip_extensions(filename: str) -> str:
        for ext in (".clean.bed.gz", ".bed.gz", ".clean.bed", ".bed"):
            if filename.endswith(ext):
                return filename[: -len(ext)]
        return filename.split(".")[0]

    pairs = [
        (_strip_extensions(Path(r[name_field]).name), float(r[score_field]))
        for r in records
    ]

    # p_value-like fields: log
    if "p_val" in score_field.lower() or "pval" in score_field.lower():
        pairs = [(a, -math.log(b) if b > 0 else float("inf")) for a, b in pairs]

    pairs.sort(key=lambda x: -x[1])
    if not pairs:
        return [], []
    keys, values = zip(*pairs)
    return list(keys), list(values)


def parse_input_file(
    path: Pathish, score_field: str = "p_value"
) -> tuple[list[str], list[float]]:
    if str(path).endswith(".json"):
        return parse_json_file(path, score_field)
    return parse_score_file(path)


# ==========================================
# Dense Heatmap
# ==========================================


@dataclass
class _CategoryMapping:
    label_to_category: Callable[[str], str]
    category_order: Sequence[str] | None = None

    def get_category(self, label: str) -> str:
        return self.label_to_category(label)


@dataclass
class _AxisConfig:
    labels: Sequence[str]
    category_mapping: _CategoryMapping | None = None


@dataclass
class _DenseHeatMatrix:
    data: NDArray[np.floating]
    row_config: _AxisConfig
    col_config: _AxisConfig
    title: str | None = None
    cmap: str | Colormap = "RdBu_r"
    cbar_label: str = "Similarity Score"
    line_color: str = "black"

    _row_order: list[int] = field(default_factory=list, repr=False)
    _col_order: list[int] = field(default_factory=list, repr=False)
    _row_breaks: list[int] = field(default_factory=list, repr=False)
    _col_breaks: list[int] = field(default_factory=list, repr=False)
    _row_categories: list[str] = field(default_factory=list, repr=False)
    _col_categories: list[str] = field(default_factory=list, repr=False)

    def _compute_order(
        self, config: _AxisConfig
    ) -> tuple[list[int], list[int], list[str]]:
        labels = list(config.labels)
        if config.category_mapping is None:
            return list(range(len(labels))), [], []
        mapping = config.category_mapping
        cat_to_indices: dict[str, list[int]] = {}
        for i, label in enumerate(labels):
            cat = mapping.get_category(label)
            cat_to_indices.setdefault(cat, []).append(i)
        if mapping.category_order is not None:
            categories = [c for c in mapping.category_order if c in cat_to_indices]
            for cat in cat_to_indices:
                if cat not in categories:
                    categories.append(cat)
        else:
            categories = sorted(cat_to_indices.keys())
        ordered: list[int] = []
        breaks: list[int] = []
        for cat in categories:
            if ordered:
                breaks.append(len(ordered))
            ordered.extend(cat_to_indices[cat])
        return ordered, breaks, categories

    def _tick_positions(
        self, n: int, breaks: list[int], categories: list[str]
    ) -> tuple[list[float], list[str]]:
        if not categories:
            return list(range(n)), []
        boundaries = [0] + breaks + [n]
        positions = [
            (boundaries[i] + boundaries[i + 1] - 1) / 2 for i in range(len(categories))
        ]
        return positions, categories

    def plot(self, fig: Figure, ax: Axes) -> None:
        self._row_order, self._row_breaks, self._row_categories = self._compute_order(
            self.row_config
        )
        self._col_order, self._col_breaks, self._col_categories = self._compute_order(
            self.col_config
        )

        plot_data = self.data[np.ix_(self._row_order, self._col_order)]
        vmin, vmax = float(np.nanmin(self.data)), float(np.nanmax(self.data))

        cmap = plt.get_cmap(self.cmap) if isinstance(self.cmap, str) else self.cmap
        cmap = cmap.copy()
        cmap.set_bad(color="white", alpha=0)

        im = ax.imshow(
            plot_data,
            aspect="auto",
            cmap=cmap,
            norm=Normalize(vmin=vmin, vmax=vmax),
            interpolation="nearest",
        )

        line_kwargs = {"color": self.line_color, "linestyle": ":", "linewidth": 1.5}
        for brk in self._row_breaks:
            ax.axhline(y=brk - 0.5, **line_kwargs)
        for brk in self._col_breaks:
            ax.axvline(x=brk - 0.5, **line_kwargs)

        if self.row_config.category_mapping:
            pos, labels = self._tick_positions(
                len(self._row_order), self._row_breaks, self._row_categories
            )
            ax.set_yticks(pos)
            ax.set_yticklabels(labels)
        else:
            ax.set_yticks(range(len(self.row_config.labels)))
            ax.set_yticklabels([self.row_config.labels[i] for i in self._row_order])

        if self.col_config.category_mapping:
            pos, labels = self._tick_positions(
                len(self._col_order), self._col_breaks, self._col_categories
            )
            ax.set_xticks(pos)
            ax.set_xticklabels(labels, rotation=45, ha="right")
        else:
            ax.set_xticks(range(len(self.col_config.labels)))
            ax.set_xticklabels(
                [self.col_config.labels[i] for i in self._col_order],
                rotation=45,
                ha="right",
            )

        if self.title:
            ax.set_title(self.title, fontsize=14, fontweight="bold")

        fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02).set_label(self.cbar_label)


# ==========================================
# RME Plotting
# ==========================================


@dataclass
class _PlotData:
    name: str
    data: np.ndarray
    ordered_cells: list[str]
    bed_to_category: dict[str, str]


def _prepare_plot_data(
    score_paths: tuple[Pathish, ...],
    names: tuple[str, ...],
    states: list[str] | None = None,
    score_field: str = "p_value",
) -> list[_PlotData]:
    col_labels = states if states is not None else chromatin_states
    result: list[_PlotData] = []

    for path, name in zip(score_paths, names):
        beds, scores = parse_input_file(path, score_field)
        cell_types: set[str] = set()
        parsed_data: dict[tuple[str, str], float] = {}

        for bed, score in zip(beds, scores):
            try:
                cell_id, state = classify_bed_file(bed)
                cell_types.add(cell_id)
                parsed_data[(cell_id, state)] = score
            except ValueError:
                continue

        if not cell_types:
            print(
                f"Warning: No valid parsed data for {name}. Skipping.", file=sys.stderr
            )
            continue

        ordered_cells = sorted(cell_types)
        data = np.full((len(ordered_cells), len(col_labels)), np.nan, dtype=np.float64)
        bed_to_category: dict[str, str] = {}

        for i, cell_id in enumerate(ordered_cells):
            bed_to_category[cell_id] = classify_cell_type(cell_id)
            for j, state in enumerate(col_labels):
                data[i, j] = parsed_data.get((cell_id, state), np.nan)

        result.append(
            _PlotData(
                name=name,
                data=data,
                ordered_cells=ordered_cells,
                bed_to_category=bed_to_category,
            )
        )

    return result


def plot_rme_similarity(
    score_paths: tuple[Pathish, ...],
    names: tuple[str, ...],
    output_path: Pathish | None = None,
    show: bool = True,
    states: list[str] | None = None,
    score_field: str = "p_value",
) -> None:
    if len(score_paths) != len(names):
        raise ValueError("Must provide the same number of score paths and names.")

    plot_data_list = _prepare_plot_data(
        score_paths, names, states=states, score_field=score_field
    )
    if not plot_data_list:
        print("No valid data to plot.", file=sys.stderr)
        return

    col_labels = states if states is not None else chromatin_states
    max_rows = max(len(pd.ordered_cells) for pd in plot_data_list)
    n_plots = len(plot_data_list)

    plt.style.use("dark_background")
    fig, axes = plt.subplots(
        1,
        n_plots,
        figsize=(8 * n_plots, max(3.0, max_rows * 0.1)),
        squeeze=False,
        layout="constrained",
    )

    for idx, pd in enumerate(plot_data_list):
        _DenseHeatMatrix(
            data=pd.data,
            row_config=_AxisConfig(
                labels=pd.ordered_cells,
                category_mapping=_CategoryMapping(
                    lambda label, btc=pd.bed_to_category: btc.get(label, "Other"),
                    cell_categories,
                ),
            ),
            col_config=_AxisConfig(labels=col_labels),
            title=pd.name,
        ).plot(fig, axes[0, idx])

    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"Saved to {output_path}", file=sys.stderr)

    if show:
        plt.show()


# ==========================================
# CLI
# ==========================================


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot comparative dense heatmaps from RME score files."
    )
    parser.add_argument(
        "-s", "--scores", nargs="+", required=True, help="Paths to score files."
    )
    parser.add_argument(
        "-n", "--names", nargs="+", required=True, help="Names for the score files."
    )
    parser.add_argument("-o", "--output", help="Output file path (e.g. plot.png).")
    parser.add_argument(
        "--no-show", action="store_true", help="Do not display interactively."
    )
    parser.add_argument(
        "--states",
        nargs="+",
        choices=chromatin_states,
        metavar="STATE",
        help=f"Chromatin states to include (default: all). Choices: {', '.join(chromatin_states)}",
    )
    parser.add_argument(
        "--score-field",
        default="p_value",
        help=(
            "Field name to use as the score when parsing JSON inputs (default: p_value). "
            "Fields containing 'p_val' are automatically negated so lower p-value → higher score."
        ),
    )
    args = parser.parse_args()

    if len(args.scores) != len(args.names):
        print("Error: --scores and --names must have the same count.", file=sys.stderr)
        sys.exit(1)

    plot_rme_similarity(
        tuple(args.scores),
        tuple(args.names),
        output_path=args.output,
        show=not args.no_show,
        states=args.states,
        score_field=args.score_field,
    )


if __name__ == "__main__":
    main()
