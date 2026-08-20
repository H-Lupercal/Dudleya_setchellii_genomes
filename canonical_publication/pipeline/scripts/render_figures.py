#!/usr/bin/env python3
"""Render fingerprinted publication figures from canonical analysis outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import tomllib
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path

import matplotlib
import numpy as np
from Bio import Phylo
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from organelle_pipeline.figures import (
    TAXON_COLORS,
    composition_counts,
    distance_aware_layout,
    major_haplotype_labels,
    side_label_layout,
    signed_fst_limit,
    support_is_strong,
    unrooted_tree_layout,
)
from organelle_pipeline.paths import repository_relative, validate_run_id
from organelle_pipeline.provenance import (
    build_stage_fingerprint_from_hashes,
    runtime_provenance,
    sha256_file,
    validate_resume,
    validate_saved_outputs,
)
from organelle_pipeline.references import read_single_fasta

FORMATS = ("png", "pdf", "svg")
ORGANELLES = ("chloroplast", "mitochondria")
FST_CMAP = LinearSegmentedColormap.from_list("signed_fst", ("#2166AC", "#F7F7F7", "#B2182B"))
CLUSTER_COLORS = (
    "#332288",
    "#88CCEE",
    "#44AA99",
    "#117733",
    "#999933",
    "#DDCC77",
    "#CC6677",
    "#882255",
    "#AA4499",
    "#661100",
    "#6699CC",
    "#AA4466",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_state(root: Path, path: Path) -> dict[str, object]:
    if not path.is_file():
        raise RuntimeError(f"Required upstream state is missing: {path.relative_to(root)}")
    state = json.loads(path.read_text())
    if state.get("status") not in {"complete", "PASS"}:
        raise RuntimeError(f"Required upstream state is not complete: {path.relative_to(root)}")
    fingerprint = state.get("fingerprint")
    if not isinstance(fingerprint, dict) or not isinstance(fingerprint.get("digest"), str):
        raise RuntimeError(f"Required upstream state has no fingerprint: {path.relative_to(root)}")
    validate_saved_outputs(root, state)
    return state


def require_declared(root: Path, path: Path, state: dict[str, object]) -> None:
    outputs = state.get("outputs")
    relative = path.relative_to(root).as_posix()
    if not isinstance(outputs, dict) or outputs.get(relative) != sha256_file(path):
        raise RuntimeError(f"Figure input is not checksum-declared by its upstream state: {relative}")


def metadata_for(path: Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(path)
    metadata = {row["sample_id"]: row for row in rows}
    if len(metadata) != len(rows):
        raise RuntimeError(f"Duplicate sample ID in metadata: {path}")
    unknown = sorted({row["species"] for row in rows} - set(TAXON_COLORS))
    if unknown:
        raise RuntimeError(f"Unknown taxa in metadata: {unknown}")
    return metadata


def bed_intervals(path: Path) -> list[tuple[int, int, str]]:
    intervals = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        intervals.append((int(fields[1]), int(fields[2]), fields[3] if len(fields) > 3 else "interval"))
    return intervals


def save_figure(figure: plt.Figure, directory: Path, figure_id: str) -> list[Path]:
    outputs = []
    for extension in FORMATS:
        path = directory / f"{figure_id}.{extension}"
        metadata_by_format: dict[str, object]
        if extension == "pdf":
            metadata_by_format = {"Creator": "canonical organelle pipeline", "CreationDate": None, "ModDate": None}
        elif extension == "svg":
            metadata_by_format = {"Creator": "canonical organelle pipeline", "Date": None}
        else:
            metadata_by_format = {"Software": "canonical organelle pipeline"}
        figure.savefig(path, dpi=300, bbox_inches="tight", metadata=metadata_by_format)
        outputs.append(path)
    plt.close(figure)
    return outputs


def add_taxon_legend(axis: plt.Axes, present: set[str], *, horizontal_bottom: bool = False) -> None:
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=TAXON_COLORS[taxon], markeredgecolor="black", label=taxon)
        for taxon in TAXON_COLORS
        if taxon in present
    ]
    if horizontal_bottom:
        axis.legend(
            handles=handles,
            title="Taxon",
            fontsize=7,
            title_fontsize=8,
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.035),
            ncol=3,
        )
    else:
        axis.legend(handles=handles, title="Taxon", fontsize=7, title_fontsize=8, frameon=False, loc="best")


def render_reference_maps(root: Path, run_id: str, output: Path) -> list[Path]:
    cp_length = len(read_single_fasta(root / "canonical_publication/references/selected/chloroplast.fa")[0][1])
    mt_length = len(read_single_fasta(root / "canonical_publication/references/selected/mitochondria.fa")[0][1])
    cp_annotations = read_tsv(root / "canonical_publication/references/annotations/chloroplast.projected.tsv")
    mt_annotations = read_tsv(root / "canonical_publication/references/annotations/mitochondria.projected.tsv")
    cp_ir = bed_intervals(root / "canonical_publication/references/masks/chloroplast_ir_copies.bed")
    cp_population = bed_intervals(root / "canonical_publication/references/masks/chloroplast_population_sites.bed")
    mt_repeat = bed_intervals(root / "canonical_publication/references/masks/mitochondria_repeat_mask.bed")
    mt_high = bed_intervals(root / f"canonical_publication/references/masks/{run_id}/mitochondria_high_confidence_sites.bed")

    figure = plt.figure(figsize=(12, 5.4), constrained_layout=True)
    cp_axis = figure.add_subplot(1, 2, 1, projection="polar")
    theta = np.linspace(0, 2 * np.pi, 720)
    cp_axis.plot(theta, np.full_like(theta, 1.0), color="#303030", linewidth=1.2)
    for row in cp_annotations:
        if row["feature_type"] != "gene":
            continue
        start = (int(row["start_1based"]) - 1) / cp_length * 2 * np.pi
        end = int(row["end_1based"]) / cp_length * 2 * np.pi
        cp_axis.plot(np.linspace(start, end, 30), np.full(30, 1.08), color="#56B4E9", linewidth=3)
    for start, end, _ in cp_ir:
        cp_axis.plot(
            np.linspace(start / cp_length * 2 * np.pi, end / cp_length * 2 * np.pi, 30),
            np.full(30, 0.88),
            color="#D55E00",
            linewidth=7,
        )
    for start, end, _ in cp_population:
        cp_axis.plot(
            np.linspace(start / cp_length * 2 * np.pi, end / cp_length * 2 * np.pi, 30),
            np.full(30, 0.75),
            color="#009E73",
            linewidth=5,
        )
    cp_axis.set_ylim(0.62, 1.18)
    cp_axis.set_yticks((0.75, 0.88, 1.08), labels=("population mask", "IR copies", "genes"), fontsize=7)
    cp_axis.set_xticks(
        np.linspace(0, 2 * np.pi, 4, endpoint=False), labels=["0", f"{cp_length // 4:,}", f"{cp_length // 2:,}", f"{3 * cp_length // 4:,}"]
    )
    cp_axis.grid(alpha=0.2)
    cp_axis.set_title(f"A  Chloroplast reference (circular; {cp_length:,} bp)", loc="left", fontweight="bold")

    mt_axis = figure.add_subplot(1, 2, 2)
    mt_axis.hlines(0, 0, mt_length, color="#303030", linewidth=1.2)
    for row in mt_annotations:
        if row["feature_type"] == "gene":
            mt_axis.hlines(0.35, int(row["start_1based"]) - 1, int(row["end_1based"]), color="#56B4E9", linewidth=4)
    for start, end, _ in mt_repeat:
        mt_axis.hlines(-0.32, start, end, color="#D55E00", linewidth=8)
    for start, end, _ in mt_high:
        mt_axis.hlines(-0.62, start, end, color="#009E73", linewidth=6)
    mt_axis.text(mt_length * 0.01, 0.4, "projected genes", fontsize=8, va="bottom")
    mt_axis.text(mt_length * 0.01, -0.27, "self-repeat exclusion", fontsize=8, va="bottom")
    mt_axis.text(mt_length * 0.01, -0.57, "read-backed high confidence", fontsize=8, va="bottom")
    mt_axis.set_xlim(0, mt_length)
    mt_axis.set_ylim(-0.85, 0.75)
    mt_axis.set_yticks([])
    mt_axis.set_xlabel("Reference coordinate (bp)")
    mt_axis.spines[["left", "right", "top"]].set_visible(False)
    mt_axis.set_title(f"B  Mitochondrial candidate (linear; {mt_length:,} bp)", loc="left", fontweight="bold")
    figure.suptitle("Reference architecture and inference masks", fontweight="bold")
    return save_figure(figure, output, "reference_callability")


def render_qc(root: Path, run_id: str, config: dict[str, object], output: Path) -> list[Path]:
    breadth = read_tsv(root / f"canonical_publication/results/qc/{run_id}/sample_breadth.tsv")
    preprocessing = read_tsv(root / f"canonical_publication/results/qc/{run_id}/read_preprocessing_summary.tsv")
    if {row["sample_id"] for row in breadth} != {row["sample_id"] for row in preprocessing}:
        raise RuntimeError("QC breadth/preprocessing sample mismatch")
    threshold = float(config["qc"]["minimum_breadth"])  # type: ignore[index]
    depth = int(config["qc"]["eligibility_depth"])  # type: ignore[index]
    figure, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    cp = np.asarray([float(row[f"cp_unique_sites_breadth_dp{depth}"]) for row in breadth])
    mt = np.asarray([float(row[f"mt_unique_sites_breadth_dp{depth}"]) for row in breadth])
    axes[0].scatter(cp, mt, s=25, color="#0072B2", edgecolor="white", linewidth=0.4, alpha=0.8)
    axes[0].axvline(threshold, color="#D55E00", linestyle="--", linewidth=1)
    axes[0].axhline(threshold, color="#D55E00", linestyle="--", linewidth=1)
    axes[0].set(xlabel=f"Chloroplast unique-site breadth at DP{depth}", ylabel=f"Mitochondrial unique-site breadth at DP{depth}")
    axes[0].set_title("A  Organelle-specific eligibility", loc="left", fontweight="bold")
    axes[0].text(0.02, 0.02, f"Eligibility threshold = {threshold:.0%}", transform=axes[0].transAxes, fontsize=8)
    retention = np.asarray([float(row["read_retention"]) for row in preprocessing]) * 100
    q20 = np.asarray([float(row["passing_q20_rate"]) for row in preprocessing]) * 100
    positions = np.arange(len(preprocessing))
    axes[1].scatter(positions, retention, s=18, label="reads retained", color="#009E73", alpha=0.8)
    axes[1].scatter(positions, q20, s=18, label="passing bases Q20", color="#CC79A7", alpha=0.8)
    axes[1].set(xlabel="Samples ordered by read retention", ylabel="Percent", ylim=(0, 102))
    order = np.argsort(retention)
    axes[1].collections[0].set_offsets(np.column_stack((positions, retention[order])))
    axes[1].collections[1].set_offsets(np.column_stack((positions, q20[order])))
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].set_title("B  Read preprocessing", loc="left", fontweight="bold")
    figure.suptitle("Preprocessing and mapping quality control", fontweight="bold")
    return save_figure(figure, output, "preprocessing_qc")


def render_pca(root: Path, run_id: str, organelle: str, output: Path) -> list[Path]:
    coordinates = read_tsv(root / f"canonical_publication/results/pca/{run_id}/{organelle}.coordinates.tsv")
    variance = read_tsv(root / f"canonical_publication/results/pca/{run_id}/{organelle}.variance.tsv")
    metadata = metadata_for(root / f"canonical_publication/metadata/qc/{run_id}/{organelle}_samples.tsv")
    if {row["sample_id"] for row in coordinates} != set(metadata):
        raise RuntimeError(f"PCA/metadata sample mismatch for {organelle}")
    values = np.asarray([[float(row["PC1"]), float(row["PC2"])] for row in coordinates])
    if not np.isfinite(values).all():
        raise RuntimeError(f"PCA contains non-finite coordinates for {organelle}")
    variance_values = {row["component"]: float(row["explained_variance_ratio"]) for row in variance}
    figure, axis = plt.subplots(figsize=(7.5, 6.5), constrained_layout=True)
    present = set()
    for row, position in zip(coordinates, values, strict=True):
        sample = row["sample_id"]
        taxon = metadata[sample]["species"]
        present.add(taxon)
        axis.scatter(*position, color=TAXON_COLORS[taxon], s=35, edgecolor="black", linewidth=0.25, alpha=0.8)
    axis.axhline(0, color="#BBBBBB", linewidth=0.5)
    axis.axvline(0, color="#BBBBBB", linewidth=0.5)
    axis.set_xlabel(f"PC1 ({variance_values['PC1'] * 100:.1f}%)")
    axis.set_ylabel(f"PC2 ({variance_values['PC2'] * 100:.1f}%)")
    axis.set_title(f"{organelle.capitalize()} PCA — haploid SNPs, MAC≥2", fontweight="bold")
    add_taxon_legend(axis, present)
    axis.text(
        0.01,
        0.01,
        "Points are samples; population identifiers and exact coordinates are in the canonical PCA table",
        transform=axis.transAxes,
        fontsize=7,
    )
    return save_figure(figure, output, f"{organelle}.pca")


def render_haplotypes(root: Path, run_id: str, organelle: str, output: Path) -> list[Path]:
    base = root / f"canonical_publication/results/haplotypes/{run_id}"
    assignments = read_tsv(base / f"{organelle}.sample_haplotypes.tsv")
    haplotypes = read_tsv(base / f"{organelle}.haplotypes.tsv")
    edges = read_tsv(base / f"{organelle}.network_edges.tsv")
    metadata = metadata_for(root / f"canonical_publication/metadata/qc/{run_id}/{organelle}_samples.tsv")
    if {row["sample_id"] for row in assignments} != set(metadata):
        raise RuntimeError(f"Haplotype/metadata sample mismatch for {organelle}")
    sample_ids: dict[str, list[str]] = defaultdict(list)
    for row in assignments:
        if row["haplotype"] != "AMBIGUOUS":
            sample_ids[row["haplotype"]].append(row["sample_id"])
    nodes = tuple(row["haplotype"] for row in haplotypes)
    edge_values = tuple((row["haplotype_1"], row["haplotype_2"], float(row["mutational_distance"])) for row in edges)
    positions = distance_aware_layout(nodes, edge_values)
    counts_by_haplotype = {row["haplotype"]: int(row["sample_count"]) for row in haplotypes}
    if set(counts_by_haplotype) != set(nodes):
        raise RuntimeError(f"Haplotype count/node mismatch for {organelle}")
    labels = set(major_haplotype_labels(counts_by_haplotype, minimum_count=5))
    label_positions = side_label_layout(positions, labels)
    figure, axis = plt.subplots(figsize=(8, 7), constrained_layout=True)
    for left, right, _mutational_distance in edge_values:
        x = (positions[left][0], positions[right][0])
        y = (positions[left][1], positions[right][1])
        axis.plot(x, y, color="#666666", linewidth=0.8, zorder=1)
    taxon_by_sample = {sample: row["species"] for sample, row in metadata.items()}
    scale = max((abs(value) for position in positions.values() for value in position), default=1.0)
    for node in nodes:
        counts = composition_counts(sample_ids[node], taxon_by_sample)
        radius = 0.012 * scale * math.sqrt(max(1, sum(counts)))
        axis.pie(
            counts,
            colors=list(TAXON_COLORS.values()),
            radius=radius,
            center=positions[node],
            wedgeprops={"edgecolor": "white", "linewidth": 0.5},
            normalize=True,
        )
    for node, (label_x, label_y, alignment) in label_positions.items():
        axis.annotate(
            f"{node} (n={counts_by_haplotype[node]})",
            xy=positions[node],
            xytext=(label_x, label_y),
            textcoords="data",
            ha=alignment,
            va="center",
            fontsize=6.5,
            fontweight="bold",
            arrowprops={"arrowstyle": "-", "color": "#777777", "linewidth": 0.55},
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 0.5},
            zorder=4,
        )
    coordinates = np.asarray(tuple(positions.values()))
    x_min, y_min = coordinates.min(axis=0)
    x_max, y_max = coordinates.max(axis=0)
    x_span = max(float(x_max - x_min), 1.0)
    y_span = max(float(y_max - y_min), 1.0)
    axis.set_xlim(x_min - 0.35 * x_span, x_max + 0.35 * x_span)
    axis.set_ylim(y_min - 0.10 * y_span, y_max + 0.10 * y_span)
    axis.set_aspect("equal")
    axis.axis("off")
    axis.set_title(f"{organelle.capitalize()} haplotype minimum-spanning network", fontweight="bold")
    add_taxon_legend(axis, {row["species"] for row in metadata.values()}, horizontal_bottom=True)
    axis.text(
        0.01,
        0.01,
        "Node area ∝ sample count; edge length reflects mutational distance; labels shown for n≥5",
        transform=axis.transAxes,
        fontsize=7.5,
    )
    return save_figure(figure, output, f"{organelle}.haplotype_network")


def render_fst(root: Path, run_id: str, organelle: str, output: Path) -> list[Path]:
    base = root / f"canonical_publication/results/popgen/{run_id}"
    populations = [row["population"] for row in read_tsv(base / f"{organelle}.population_summary.tsv")]
    if len(populations) != len(set(populations)):
        raise RuntimeError(f"Duplicate population in {organelle} population summary")
    index = {population: number for number, population in enumerate(populations)}
    matrix = np.full((len(populations), len(populations)), np.nan)
    np.fill_diagonal(matrix, 0.0)
    for row in read_tsv(base / f"{organelle}.pairwise_hudson_fst.tsv"):
        left, right = index[row["population_1"]], index[row["population_2"]]
        value = float(row["hudson_fst"])
        matrix[left, right] = value
        matrix[right, left] = value
    off_diagonal = matrix.copy()
    np.fill_diagonal(off_diagonal, np.nan)
    limit = signed_fst_limit(off_diagonal)
    cmap = FST_CMAP.copy()
    cmap.set_bad("#D9D9D9")
    figure_size = max(7.0, min(13.0, 4.5 + 0.22 * len(populations)))
    figure, axis = plt.subplots(figsize=(figure_size, figure_size), constrained_layout=True)
    image = axis.imshow(matrix, cmap=cmap, norm=TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit), interpolation="none")
    axis.set_xticks(range(len(populations)), labels=populations, rotation=90, fontsize=6)
    axis.set_yticks(range(len(populations)), labels=populations, fontsize=6)
    axis.set_title(f"{organelle.capitalize()} signed Hudson FST", fontweight="bold")
    colorbar = figure.colorbar(image, ax=axis, shrink=0.75)
    colorbar.set_label("Hudson FST (ratio of sums; negative estimates retained)")
    return save_figure(figure, output, f"{organelle}.signed_hudson_fst")


def render_tree(root: Path, run_id: str, organelle: str, output: Path) -> list[Path]:
    path = root / f"canonical_publication/results/trees/{run_id}/{organelle}.primary.treefile"
    tree = Phylo.read(path, "newick")
    metadata = metadata_for(root / f"canonical_publication/metadata/qc/{run_id}/{organelle}_samples.tsv")
    terminals = tree.get_terminals()
    observed = [terminal.name for terminal in terminals]
    if len(observed) != len(set(observed)) or set(observed) != set(metadata):
        raise RuntimeError(f"Tree/metadata sample mismatch for {organelle}")
    positions = unrooted_tree_layout(tree)
    figure, axis = plt.subplots(figsize=(10, 9), constrained_layout=True)
    for parent in tree.find_clades(order="level"):
        for child in parent.clades:
            x = (positions[id(parent)][0], positions[id(child)][0])
            y = (positions[id(parent)][1], positions[id(child)][1])
            support = child.name if not child.is_terminal() else None
            strong = support_is_strong(support)
            axis.plot(x, y, color="#222222" if strong else "#AAAAAA", linewidth=1.25 if strong else 0.45, zorder=1)
    present = set()
    for terminal in terminals:
        sample = terminal.name
        assert sample is not None
        taxon = metadata[sample]["species"]
        present.add(taxon)
        x, y = positions[id(terminal)]
        axis.scatter(x, y, s=14, color=TAXON_COLORS[taxon], edgecolor="black", linewidth=0.2, zorder=3)
    axis.set_aspect("equal")
    axis.axis("off")
    axis.set_title(f"{organelle.capitalize()} primary ML tree — unrooted", fontweight="bold")
    add_taxon_legend(axis, present)
    axis.text(
        0.01,
        0.01,
        "Dark branches: SH-aLRT≥80 and UFBoot≥95; tip IDs and exact supports are retained in the canonical treefile",
        transform=axis.transAxes,
        fontsize=7.5,
    )
    return save_figure(figure, output, f"{organelle}.primary_unrooted_tree")


def render_admixture(root: Path, run_id: str, organelle: str, output: Path) -> list[Path]:
    base = root / f"canonical_publication/results/supplement/{run_id}/admixture/{organelle}"
    replicates = read_tsv(base / "replicate_cv.tsv")
    summary = read_tsv(base / "k_summary.tsv")
    selected = read_tsv(base / "selected_solution.tsv")
    order = read_tsv(base / "sample_order.tsv")
    if len(selected) != 1:
        raise RuntimeError(f"Expected one selected ADMIXTURE solution for {organelle}")
    selected_k = int(selected[0]["selected_k"])
    q_path = root / selected[0]["q_path"]
    q = np.loadtxt(q_path, ndmin=2)
    sample_ids = [row["sample_id"] for row in order]
    metadata = metadata_for(root / f"canonical_publication/metadata/qc/{run_id}/{organelle}_samples.tsv")
    if len(sample_ids) != len(set(sample_ids)) or set(sample_ids) != set(metadata):
        raise RuntimeError(f"ADMIXTURE/metadata sample mismatch for {organelle}")
    if q.shape != (len(sample_ids), selected_k) or not np.isfinite(q).all() or not np.allclose(q.sum(axis=1), 1.0, atol=1e-5):
        raise RuntimeError(f"Malformed selected ADMIXTURE Q matrix for {organelle}")
    figure, axes = plt.subplots(2, 1, figsize=(11, 7.5), gridspec_kw={"height_ratios": [1, 1.25]}, constrained_layout=True)
    by_k: dict[int, list[float]] = defaultdict(list)
    for row in replicates:
        by_k[int(row["k"])].append(float(row["cv_error"]))
    for k, values in sorted(by_k.items()):
        offsets = np.linspace(-0.08, 0.08, len(values)) if len(values) > 1 else np.asarray([0.0])
        axes[0].scatter(k + offsets, values, s=16, color="#777777", alpha=0.65)
    means = {int(row["k"]): float(row["mean_cv_error"]) for row in summary}
    axes[0].plot(sorted(means), [means[k] for k in sorted(means)], marker="o", color="#332288", linewidth=1.5, label="replicate mean")
    axes[0].axvline(selected_k, color="#D55E00", linestyle="--", linewidth=1, label=f"selected K={selected_k}")
    axes[0].set(xlabel="K", ylabel="Cross-validation error")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].set_title("A  Replicate cross-validation sensitivity", loc="left", fontweight="bold")
    x = np.arange(len(sample_ids))
    bottom = np.zeros(len(sample_ids))
    for cluster in range(selected_k):
        axes[1].bar(x, q[:, cluster], bottom=bottom, width=1.0, color=CLUSTER_COLORS[cluster], linewidth=0)
        bottom += q[:, cluster]
    boundaries = []
    previous = None
    for number, sample in enumerate(sample_ids):
        population = metadata[sample]["popcode"]
        if previous is not None and population != previous:
            boundaries.append(number - 0.5)
        previous = population
    for boundary in boundaries:
        axes[1].axvline(boundary, color="white", linewidth=0.5)
    axes[1].set(xlim=(-0.5, len(sample_ids) - 0.5), ylim=(0, 1), ylabel="Q", xlabel="Samples in canonical Q-matrix order")
    axes[1].set_xticks([])
    axes[1].set_title("B  Selected descriptive clustering", loc="left", fontweight="bold")
    boundary_note = "boundary optimum; interpret sensitivity cautiously" if selected[0]["boundary_optimum"] == "yes" else "interior optimum"
    figure.suptitle(
        f"Supplementary {organelle} ADMIXTURE sensitivity — {boundary_note}\n"
        "Linked haploid organelle markers; Q values are not ancestry proportions",
        fontweight="bold",
        fontsize=11,
    )
    return save_figure(figure, output, f"{organelle}.supplementary_admixture")


def main() -> int:
    args = parse_args()
    validate_run_id(args.run_id)
    root = args.repository_root.resolve()
    config_path = root / repository_relative(args.config, root)
    config = tomllib.loads(config_path.read_text())
    run_state = root / "canonical_publication/provenance/runs" / args.run_id
    state_path = run_state / "figures.json"
    figure_dir = root / "canonical_publication/reports/figures" / args.run_id
    staging = root / "canonical_publication/work" / args.run_id / "figure_staging"
    state_paths = {
        "references": run_state / "references.json",
        "qc": run_state / "qc.json",
        "tree_reproducibility": run_state / "tree_reproducibility.json",
        **{f"pca:{organelle}": run_state / f"pca/{organelle}.json" for organelle in ORGANELLES},
        **{f"haplotypes:{organelle}": run_state / f"haplotypes/{organelle}.json" for organelle in ORGANELLES},
        **{f"popgen:{organelle}": run_state / f"popgen/{organelle}.json" for organelle in ORGANELLES},
        **{f"trees:{organelle}": run_state / f"trees/{organelle}.json" for organelle in ORGANELLES},
        **{f"admixture:{organelle}": run_state / f"admixture/{organelle}.json" for organelle in ORGANELLES},
    }
    states = {label: load_state(root, path) for label, path in state_paths.items()}
    input_paths: dict[str, Path] = {
        "config": config_path,
        "cp_fasta": root / "canonical_publication/references/selected/chloroplast.fa",
        "mt_fasta": root / "canonical_publication/references/selected/mitochondria.fa",
        "cp_annotation": root / "canonical_publication/references/annotations/chloroplast.projected.tsv",
        "mt_annotation": root / "canonical_publication/references/annotations/mitochondria.projected.tsv",
        "cp_ir": root / "canonical_publication/references/masks/chloroplast_ir_copies.bed",
        "cp_population": root / "canonical_publication/references/masks/chloroplast_population_sites.bed",
        "mt_repeat": root / "canonical_publication/references/masks/mitochondria_repeat_mask.bed",
        "mt_high": root / f"canonical_publication/references/masks/{args.run_id}/mitochondria_high_confidence_sites.bed",
        "breadth": root / f"canonical_publication/results/qc/{args.run_id}/sample_breadth.tsv",
        "preprocessing": root / f"canonical_publication/results/qc/{args.run_id}/read_preprocessing_summary.tsv",
    }
    for organelle in ORGANELLES:
        input_paths.update(
            {
                f"metadata:{organelle}": root / f"canonical_publication/metadata/qc/{args.run_id}/{organelle}_samples.tsv",
                f"pca_coordinates:{organelle}": root / f"canonical_publication/results/pca/{args.run_id}/{organelle}.coordinates.tsv",
                f"pca_variance:{organelle}": root / f"canonical_publication/results/pca/{args.run_id}/{organelle}.variance.tsv",
                f"hap_assignments:{organelle}": root
                / f"canonical_publication/results/haplotypes/{args.run_id}/{organelle}.sample_haplotypes.tsv",
                f"haplotypes:{organelle}": root / f"canonical_publication/results/haplotypes/{args.run_id}/{organelle}.haplotypes.tsv",
                f"hap_edges:{organelle}": root / f"canonical_publication/results/haplotypes/{args.run_id}/{organelle}.network_edges.tsv",
                f"population_summary:{organelle}": root
                / f"canonical_publication/results/popgen/{args.run_id}/{organelle}.population_summary.tsv",
                f"fst:{organelle}": root / f"canonical_publication/results/popgen/{args.run_id}/{organelle}.pairwise_hudson_fst.tsv",
                f"tree:{organelle}": root / f"canonical_publication/results/trees/{args.run_id}/{organelle}.primary.treefile",
                f"admixture_replicates:{organelle}": root
                / f"canonical_publication/results/supplement/{args.run_id}/admixture/{organelle}/replicate_cv.tsv",
                f"admixture_summary:{organelle}": root
                / f"canonical_publication/results/supplement/{args.run_id}/admixture/{organelle}/k_summary.tsv",
                f"admixture_order:{organelle}": root
                / f"canonical_publication/results/supplement/{args.run_id}/admixture/{organelle}/sample_order.tsv",
                f"admixture_selected:{organelle}": root
                / f"canonical_publication/results/supplement/{args.run_id}/admixture/{organelle}/selected_solution.tsv",
            }
        )
    for key, path in input_paths.items():
        if not path.is_file():
            raise RuntimeError(f"Required figure input is missing: {key} -> {path.relative_to(root)}")
    for path in input_paths.values():
        if path == config_path:
            continue
        relative = path.relative_to(root).as_posix()
        declaring = [state for state in states.values() if relative in state.get("outputs", {})]
        if len(declaring) != 1:
            raise RuntimeError(f"Figure input must be declared by exactly one upstream state: {relative}")
        require_declared(root, path, declaring[0])
    selected_q_paths = []
    for organelle in ORGANELLES:
        selected = read_tsv(input_paths[f"admixture_selected:{organelle}"])
        if len(selected) != 1:
            raise RuntimeError(f"Expected one selected ADMIXTURE row for {organelle}")
        q_path = root / selected[0]["q_path"]
        require_declared(root, q_path, states[f"admixture:{organelle}"])
        selected_q_paths.append(q_path)
    declared = {
        **runtime_provenance(
            root,
            {
                "biopython": ("python", "-c", "import Bio; print(Bio.__version__)"),
                "matplotlib": ("python", "-c", "import matplotlib; print(matplotlib.__version__)"),
                "networkx": ("python", "-c", "import networkx; print(networkx.__version__)"),
                "numpy": ("python", "-c", "import numpy; print(numpy.__version__)"),
            },
        ),
        **{path.relative_to(root).as_posix(): sha256_file(path) for path in (*input_paths.values(), *selected_q_paths)},
        **{path.relative_to(root).as_posix(): sha256_file(path) for path in state_paths.values()},
    }
    fingerprint = build_stage_fingerprint_from_hashes(
        "publication_figures",
        declared,
        {label: state["fingerprint"]["digest"] for label, state in states.items()},  # type: ignore[index]
        ["deterministic canonical renderer; PNG/PDF/SVG; unrooted tree layout; signed FST centered at zero"],
    )
    if args.resume and state_path.exists():
        saved = json.loads(state_path.read_text())
        validate_resume(saved["fingerprint"]["digest"], fingerprint)
        validate_saved_outputs(root, saved)
        print("resume-valid publication figures")
        return 0
    if state_path.exists() or figure_dir.exists() or staging.exists():
        raise RuntimeError("Existing unvalidated figure output; preserve it and use a new run ID")
    staging.mkdir(parents=True)
    matplotlib.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.hashsalt": "canonical-publication-figures",
        }
    )
    try:
        figure_paths = []
        figure_paths.extend(render_reference_maps(root, args.run_id, staging))
        figure_paths.extend(render_qc(root, args.run_id, config, staging))
        for organelle in ORGANELLES:
            figure_paths.extend(render_pca(root, args.run_id, organelle, staging))
            figure_paths.extend(render_haplotypes(root, args.run_id, organelle, staging))
            figure_paths.extend(render_fst(root, args.run_id, organelle, staging))
            figure_paths.extend(render_tree(root, args.run_id, organelle, staging))
            figure_paths.extend(render_admixture(root, args.run_id, organelle, staging))
        manifest = staging / "figure_manifest.tsv"
        with manifest.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["figure_id", "format", "path", "sha256"])
            for path in sorted(figure_paths):
                writer.writerow(
                    [
                        path.stem,
                        path.suffix.removeprefix("."),
                        f"canonical_publication/reports/figures/{args.run_id}/{path.name}",
                        sha256_file(path),
                    ]
                )
        figure_dir.parent.mkdir(parents=True, exist_ok=True)
        os.replace(staging, figure_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    final_paths = sorted(path for path in figure_dir.iterdir() if path.is_file())
    outputs = {path.relative_to(root).as_posix(): sha256_file(path) for path in final_paths}
    state_path.parent.mkdir(parents=True, exist_ok=True)
    state_path.write_text(
        json.dumps(
            {
                "status": "complete",
                "figure_count": len(figure_paths) // len(FORMATS),
                "formats": list(FORMATS),
                "fingerprint": asdict(fingerprint),
                "outputs": outputs,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(f"completed {len(figure_paths) // len(FORMATS)} publication figures in {len(FORMATS)} formats")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
