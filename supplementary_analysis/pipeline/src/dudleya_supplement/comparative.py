"""Comparative-analysis invariants and tree operations."""

from __future__ import annotations

import copy
from pathlib import Path

from Bio import Phylo


def expected_pair_count(population_count: int) -> int:
    if population_count < 0:
        raise ValueError("population_count must be nonnegative")
    return population_count * (population_count - 1) // 2


def validate_resampling_spec(*, site_draws: int, site_seed: int, pi_draws: int, pi_seed: int, common_n: int) -> None:
    observed = (site_draws, site_seed, pi_draws, pi_seed, common_n)
    required = (1000, 424200, 1000, 424201, 4)
    if observed != required:
        raise ValueError(f"Resampling specification changed: observed={observed}, required={required}")


def _support_pair(name: str | None) -> tuple[float, float]:
    if not name or "/" not in name:
        return 0.0, 0.0
    try:
        left, right = name.split("/", 1)
        return float(left), float(right)
    except ValueError:
        return 0.0, 0.0


def supported_contracted_tree(path: Path, taxa: set[str]):
    """Read, prune, and contract an unrooted IQ-TREE topology."""
    tree = Phylo.read(path, "newick")
    for terminal in list(tree.get_terminals()):
        if terminal.name not in taxa:
            tree.prune(terminal)
    changed = True
    while changed:
        changed = False
        for clade in list(tree.get_nonterminals(order="postorder")):
            if clade is tree.root:
                continue
            sh_alrt, ufboot = _support_pair(clade.name)
            if sh_alrt < 80 or ufboot < 95:
                tree.collapse(clade)
                changed = True
                break
    tree.rooted = False
    return tree


def unrooted_splits(tree) -> set[frozenset[str]]:
    terminals = frozenset(tip.name for tip in tree.get_terminals())
    splits: set[frozenset[str]] = set()
    for clade in tree.get_nonterminals():
        side = frozenset(tip.name for tip in clade.get_terminals())
        if 1 < len(side) < len(terminals) - 1:
            complement = terminals - side
            splits.add(side if tuple(sorted(side)) <= tuple(sorted(complement)) else complement)
    return splits


def normalized_unrooted_rf(left, right) -> tuple[int, int, float]:
    left_taxa = {tip.name for tip in left.get_terminals()}
    right_taxa = {tip.name for tip in right.get_terminals()}
    if left_taxa != right_taxa:
        raise ValueError("RF trees do not have identical taxon sets")
    a, b = unrooted_splits(left), unrooted_splits(right)
    numerator = len(a - b) + len(b - a)
    denominator = len(a) + len(b)
    return numerator, denominator, numerator / denominator if denominator else 0.0


def prune_copy(tree, taxa: set[str]):
    result = copy.deepcopy(tree)
    for terminal in list(result.get_terminals()):
        if terminal.name not in taxa:
            result.prune(terminal)
    return result
