"""Deterministic visual semantics for canonical publication figures."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence

import networkx as nx
import numpy as np
from Bio.Phylo.BaseTree import Tree

TAXON_COLORS: dict[str, str] = {
    "D. abramsii ssp. abramsii": "#0072B2",
    "D. abramsii ssp. bettinae": "#E69F00",
    "D. abramsii ssp. murina": "#009E73",
    "D. cymosa": "#CC79A7",
    "D. setchellii": "#D55E00",
}


def signed_fst_limit(values: np.ndarray | Sequence[Sequence[float]]) -> float:
    """Return a finite symmetric color limit without discarding negative FST."""

    matrix = np.asarray(values, dtype=float)
    finite = np.abs(matrix[np.isfinite(matrix)])
    if finite.size == 0:
        raise ValueError("signed FST display requires at least one finite estimate")
    maximum = float(finite.max())
    return maximum if maximum > 0 else 1.0


def composition_counts(sample_ids: Iterable[str], sample_taxa: Mapping[str, str]) -> tuple[int, ...]:
    """Count samples in the fixed canonical taxon order."""

    counts = {taxon: 0 for taxon in TAXON_COLORS}
    for sample_id in sample_ids:
        if sample_id not in sample_taxa:
            raise ValueError(f"Missing taxon metadata for sample: {sample_id}")
        taxon = sample_taxa[sample_id]
        if taxon not in counts:
            raise ValueError(f"Unknown taxon for sample {sample_id}: {taxon}")
        counts[taxon] += 1
    return tuple(counts[taxon] for taxon in TAXON_COLORS)


def major_haplotype_labels(sample_counts: Mapping[str, int], minimum_count: int = 5) -> tuple[str, ...]:
    """Select haplotypes for labels using a declared, non-data-adaptive count rule."""

    if minimum_count < 1:
        raise ValueError("minimum haplotype label count must be positive")
    if any(count < 1 for count in sample_counts.values()):
        raise ValueError("haplotype sample counts must be positive")
    return tuple(sorted(haplotype for haplotype, count in sample_counts.items() if count >= minimum_count))


def side_label_layout(positions: Mapping[str, tuple[float, float]], labels: Iterable[str]) -> dict[str, tuple[float, float, str]]:
    """Place dense graph labels in balanced side columns with deterministic order."""

    selected = sorted(set(labels), key=lambda label: (positions[label][0], positions[label][1], label))
    if not selected:
        return {}
    coordinates = np.asarray(tuple(positions.values()), dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 2 or not np.isfinite(coordinates).all():
        raise ValueError("label layout requires finite two-dimensional positions")
    x_min, y_min = coordinates.min(axis=0)
    x_max, y_max = coordinates.max(axis=0)
    x_span = max(float(x_max - x_min), np.finfo(float).eps)
    split = (len(selected) + 1) // 2
    groups = ((selected[:split], float(x_min - 0.08 * x_span), "right"), (selected[split:], float(x_max + 0.08 * x_span), "left"))
    result: dict[str, tuple[float, float, str]] = {}
    for group, x_label, alignment in groups:
        ordered = sorted(group, key=lambda label: (positions[label][1], label))
        y_values = np.linspace(y_min, y_max, len(ordered)) if len(ordered) > 1 else np.asarray([positions[ordered[0]][1]])
        for label, y_label in zip(ordered, y_values, strict=True):
            result[label] = (x_label, float(y_label), alignment)
    return result


def support_is_strong(label: str | None, sh_alrt: float = 80.0, ufboot: float = 95.0) -> bool:
    """Return whether an IQ-TREE SH-aLRT/UFBoot label meets both thresholds."""

    if label is None:
        return False
    fields = label.split("/")
    if len(fields) != 2:
        return False
    try:
        sh_value, ufboot_value = (float(value) for value in fields)
    except ValueError:
        return False
    return sh_value >= sh_alrt and ufboot_value >= ufboot


def _scaled_kamada_kawai(graph: nx.Graph) -> dict[object, tuple[float, float]]:
    if not graph:
        return {}
    if not nx.is_connected(graph):
        raise ValueError("distance-aware layout requires a connected graph")
    distances = dict(nx.all_pairs_dijkstra_path_length(graph, weight="distance"))
    positions = nx.kamada_kawai_layout(graph, dist=distances, weight=None, dim=2, scale=1.0)
    return {node: (float(position[0]), float(position[1])) for node, position in sorted(positions.items(), key=lambda item: str(item[0]))}


def distance_aware_layout(
    nodes: Iterable[str],
    edges: Iterable[tuple[str, str, float]],
) -> dict[str, tuple[float, float]]:
    """Lay out a connected haplotype graph using mutational path distances."""

    graph = nx.Graph()
    graph.add_nodes_from(sorted(nodes))
    for left, right, distance in edges:
        value = float(distance)
        if not math.isfinite(value) or value <= 0:
            raise ValueError("mutational distances must be finite and positive")
        graph.add_edge(left, right, distance=value)
    return _scaled_kamada_kawai(graph)


def unrooted_tree_layout(tree: Tree) -> dict[int, tuple[float, float]]:
    """Return an orientation-free layout based only on tree branch distances."""

    graph = nx.Graph()
    for parent in tree.find_clades(order="level"):
        parent_id = id(parent)
        graph.add_node(parent_id)
        for child in parent.clades:
            length = 1.0 if child.branch_length is None else float(child.branch_length)
            if not math.isfinite(length) or length < 0:
                raise ValueError("tree branch lengths must be finite and nonnegative")
            graph.add_edge(parent_id, id(child), distance=max(length, np.finfo(float).eps))
    return _scaled_kamada_kawai(graph)
