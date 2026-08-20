"""Scientifically meaningful comparison of fixed-seed unrooted IQ-TREE results."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import Phylo

from organelle_pipeline.analysis import is_strong_iqtree_support, parse_iqtree_support


@dataclass(frozen=True)
class TreeReproducibility:
    taxa_equal: bool
    canonical_internal_splits: int
    replicate_internal_splits: int
    full_unrooted_rf: int
    canonical_strong_splits: int
    replicate_strong_splits: int
    strong_split_symmetric_difference: int
    max_shared_strong_branch_length_difference: float

    @property
    def strong_topology_reproduced(self) -> bool:
        return self.taxa_equal and self.strong_split_symmetric_difference == 0


def _split_maps(path: Path) -> tuple[frozenset[str], dict[frozenset[str], tuple[float, bool]]]:
    tree = Phylo.read(path, "newick")
    terminal_names = [terminal.name for terminal in tree.get_terminals()]
    if None in terminal_names or len(terminal_names) != len(set(terminal_names)):
        raise ValueError(f"Tree has missing or duplicate terminal names: {path}")
    taxa = frozenset(terminal_names)
    splits: dict[frozenset[str], tuple[float, bool]] = {}
    for clade in tree.get_nonterminals(order="postorder"):
        side = frozenset(terminal.name for terminal in clade.get_terminals())
        complement = taxa - side
        if len(side) < 2 or len(complement) < 2:
            continue
        canonical = side if (len(side), sorted(side)) <= (len(complement), sorted(complement)) else complement
        try:
            strong = is_strong_iqtree_support(parse_iqtree_support(clade.name, clade.confidence))
        except ValueError:
            strong = False
        splits[canonical] = (clade.branch_length or 0.0, strong)
    return taxa, splits


def compare_unrooted_trees(canonical_path: Path, replicate_path: Path) -> TreeReproducibility:
    """Compare full and strongly supported unrooted split sets."""

    canonical_taxa, canonical = _split_maps(canonical_path)
    replicate_taxa, replicate = _split_maps(replicate_path)
    taxa_equal = canonical_taxa == replicate_taxa
    if not taxa_equal:
        return TreeReproducibility(False, len(canonical), len(replicate), -1, 0, 0, -1, float("nan"))
    canonical_keys = set(canonical)
    replicate_keys = set(replicate)
    canonical_strong = {key for key, (_, strong) in canonical.items() if strong}
    replicate_strong = {key for key, (_, strong) in replicate.items() if strong}
    shared_strong = canonical_strong & replicate_strong
    branch_differences = [abs(canonical[key][0] - replicate[key][0]) for key in shared_strong]
    return TreeReproducibility(
        taxa_equal=True,
        canonical_internal_splits=len(canonical_keys),
        replicate_internal_splits=len(replicate_keys),
        full_unrooted_rf=len(canonical_keys ^ replicate_keys),
        canonical_strong_splits=len(canonical_strong),
        replicate_strong_splits=len(replicate_strong),
        strong_split_symmetric_difference=len(canonical_strong ^ replicate_strong),
        max_shared_strong_branch_length_difference=max(branch_differences, default=0.0),
    )
