#!/usr/bin/env python3
"""Build final visual summaries from the completed organelle pipeline outputs."""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("submission_figures") / ".matplotlib"))

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from Bio import Phylo
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "submission_figures"
FULL = ROOT / "full_pipeline_run" / "results"

ORGANELLES = {
    "cpDNA": {
        "fst": FULL / "17_population_genetics" / "cpDNA.primary.population_genetics.pairwise_fst.tsv",
        "tree": FULL / "19_bootstrap_phylogenetic_tree" / "cpDNA.primary.iqtree_ml.treefile",
        "pca": FULL / "15_pca" / "cpDNA.primary.pca.png",
        "admixture": FULL / "18_admixture_replicates" / "cpDNA.primary.bestK8.structure.png",
        "map": ROOT / "genome_maps" / "cpDNA.circular_genome_map.png",
    },
    "mtDNA": {
        "fst": FULL / "17_population_genetics" / "mtDNA.primary.population_genetics.pairwise_fst.tsv",
        "tree": FULL / "19_bootstrap_phylogenetic_tree" / "mtDNA.primary.iqtree_ml.treefile",
        "pca": FULL / "15_pca" / "mtDNA.primary.pca.png",
        "admixture": FULL / "18_admixture_replicates" / "mtDNA.primary.bestK8.structure.png",
        "map": ROOT / "genome_maps" / "mtDNA.circular_genome_map.png",
    },
}


def read_sample_groups() -> dict[str, str]:
    path = FULL / "07_downstream_sample_set" / "included_samples.tsv"
    groups: dict[str, str] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample_id = row["sample_id"]
            popcode = row.get("popcode") or ""
            profile = row.get("naming_profile") or "unknown"
            groups[sample_id] = popcode if popcode else profile
    return groups


def group_palette(groups: dict[str, str]) -> dict[str, str]:
    unique = sorted(set(groups.values()))
    palette = sns.color_palette("tab20", n_colors=20) + sns.color_palette("Set3", n_colors=12)
    return {group: palette[i % len(palette)] for i, group in enumerate(unique)}


def build_fst_heatmap(organelle: str, fst_path: Path) -> Path:
    df = pd.read_csv(fst_path, sep="\t")
    populations = sorted(set(df["population_1"]).union(df["population_2"]))
    matrix = pd.DataFrame(0.0, index=populations, columns=populations)
    for row in df.itertuples(index=False):
        fst = float(row.fst)
        matrix.loc[row.population_1, row.population_2] = fst
        matrix.loc[row.population_2, row.population_1] = fst

    fig_size = max(10, len(populations) * 0.34)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=180)
    sns.heatmap(
        matrix,
        ax=ax,
        cmap="viridis",
        vmin=0,
        vmax=1,
        square=True,
        cbar_kws={"label": "Pairwise Fst"},
        xticklabels=True,
        yticklabels=True,
    )
    ax.set_title(f"{organelle} pairwise Fst among populations", fontsize=14, pad=12)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=90, labelsize=6)
    ax.tick_params(axis="y", labelsize=6)
    fig.tight_layout()

    out = OUT_DIR / f"{organelle}.pairwise_fst_heatmap.png"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)
    return out


def build_annotated_tree(organelle: str, tree_path: Path, groups: dict[str, str], colors: dict[str, str]) -> Path:
    tree = Phylo.read(tree_path, "newick")
    terminal_count = len(tree.get_terminals())
    fig_height = max(18, terminal_count * 0.075)
    fig, ax = plt.subplots(figsize=(12, fig_height), dpi=180)

    def label_color(label: str) -> object:
        return colors.get(groups.get(label, "unknown"), "#444444")

    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        show_confidence=False,
        label_colors=label_color,
        label_func=lambda clade: clade.name if clade.is_terminal() else None,
    )
    ax.set_title(f"{organelle} maximum-likelihood tree colored by population group", fontsize=14)
    ax.set_xlabel("Substitutions per site")
    ax.tick_params(axis="y", labelsize=4)

    legend_groups = sorted(set(groups.values()))
    handles = [
        Line2D([0], [0], color=colors[group], lw=3, label=group)
        for group in legend_groups[:28]
    ]
    if handles:
        ax.legend(
            handles=handles,
            title="Group",
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            fontsize=6,
            title_fontsize=7,
            frameon=False,
        )
    fig.tight_layout()

    out = OUT_DIR / f"{organelle}.population_colored_tree.png"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)
    return out


def add_image(ax, path: Path, title: str) -> None:
    image = mpimg.imread(path)
    ax.imshow(image)
    ax.set_title(title, fontsize=12, pad=8)
    ax.axis("off")


def build_submission_panel() -> Path:
    fig, axes = plt.subplots(3, 2, figsize=(16, 22), dpi=180)
    add_image(axes[0, 0], ORGANELLES["cpDNA"]["map"], "A. cpDNA circular genome map")
    add_image(axes[0, 1], ORGANELLES["mtDNA"]["map"], "B. mtDNA circular genome map")
    add_image(axes[1, 0], ORGANELLES["cpDNA"]["pca"], "C. cpDNA PCA")
    add_image(axes[1, 1], ORGANELLES["mtDNA"]["pca"], "D. mtDNA PCA")
    add_image(axes[2, 0], ORGANELLES["cpDNA"]["admixture"], "E. cpDNA ADMIXTURE, best K=8")
    add_image(axes[2, 1], ORGANELLES["mtDNA"]["admixture"], "F. mtDNA ADMIXTURE, best K=8")
    fig.suptitle(
        "Dudleya setchellii cpDNA and mtDNA population-genomics summary",
        fontsize=18,
        weight="bold",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.985))

    out = OUT_DIR / "dudleya_organelle_submission_panel.png"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return out


def write_summary(outputs: dict[str, Path]) -> None:
    with (OUT_DIR / "submission_figures_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["figure", "path"])
        for name, path in outputs.items():
            writer.writerow([name, path.relative_to(ROOT)])


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    groups = read_sample_groups()
    colors = group_palette(groups)
    outputs: dict[str, Path] = {}
    for organelle, paths in ORGANELLES.items():
        outputs[f"{organelle}_fst_heatmap"] = build_fst_heatmap(organelle, paths["fst"])
        outputs[f"{organelle}_population_colored_tree"] = build_annotated_tree(
            organelle,
            paths["tree"],
            groups,
            colors,
        )
    outputs["submission_panel"] = build_submission_panel()
    write_summary(outputs)


if __name__ == "__main__":
    main()
