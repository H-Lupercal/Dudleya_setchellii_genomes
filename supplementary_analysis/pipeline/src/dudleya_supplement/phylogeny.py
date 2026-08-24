"""Likelihood-mapping decisions and bootstrap split handling."""

from __future__ import annotations

import re
from pathlib import Path


def parse_split_nexus(path: Path | str) -> tuple[tuple[str, ...], tuple[tuple[float, frozenset[str]], ...]]:
    text = Path(path).read_text()
    taxa_block = re.search(r"TAXLABELS(.*?);", text, re.IGNORECASE | re.DOTALL)
    matrix_block = re.search(r"BEGIN\s+Splits;.*?MATRIX(.*?);", text, re.IGNORECASE | re.DOTALL)
    if taxa_block is None or matrix_block is None:
        raise ValueError(f"Malformed split NEXUS: {path}")
    indexed = [(int(number), label) for number, label in re.findall(r"\[(\d+)\]\s+'([^']+)'", taxa_block.group(1))]
    indexed.sort()
    taxa = tuple(label for _, label in indexed)
    splits: list[tuple[float, frozenset[str]]] = []
    for raw in matrix_block.group(1).splitlines():
        line = raw.strip().rstrip(",")
        if not line:
            continue
        fields = line.split()
        weight = float(fields[0]) / 100.0
        members = frozenset(taxa[int(field) - 1] for field in fields[1:])
        if 1 < len(members) < len(taxa) - 1:
            splits.append((weight, members))
    return taxa, tuple(splits)


def splits_incompatible(left: frozenset[str], right: frozenset[str], taxa: frozenset[str]) -> bool:
    return bool(left & right) and bool(left - right) and bool(right - left) and bool(taxa - (left | right))


def supported_incompatible_pair(
    taxa: tuple[str, ...],
    splits: tuple[tuple[float, frozenset[str]], ...],
    *,
    minimum_frequency: float = 0.20,
) -> tuple[tuple[float, frozenset[str]], tuple[float, frozenset[str]]] | None:
    universe = frozenset(taxa)
    supported = [item for item in splits if item[0] >= minimum_frequency]
    for index, left in enumerate(supported):
        for right in supported[index + 1 :]:
            if splits_incompatible(left[1], right[1], universe):
                return left, right
    return None


def likelihood_decision(
    *,
    center_fraction: float,
    side_fraction: float,
    has_supported_conflict: bool,
    center_limit: float = 0.15,
    side_trigger: float = 0.20,
) -> str:
    if center_fraction > center_limit:
        return "INSUFFICIENT_INFORMATION"
    if side_fraction > side_trigger and has_supported_conflict:
        return "RUN_NEIGHBORNET"
    return "TREE_LIKE_NO_NETWORK"


def parse_identical_sequence_map(text: str) -> dict[str, str]:
    """Return only tips IQ-TREE collapsed and later restored at zero length.

    IQ-TREE also reports identical sequences that it explicitly keeps. Those
    tips are already in the representative tree and must not be labeled as
    restored collapsed samples.
    """
    return {duplicate: representative for duplicate, representative in re.findall(r"NOTE:\s+(\S+) \(identical to (\S+)\) is ignored", text)}
