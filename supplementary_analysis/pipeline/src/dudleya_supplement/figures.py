"""Supplementary figure-family contract and shared rendering helpers."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence

FIGURE_FAMILIES = (
    "robustness",
    "phylogenetic_information",
    "technical_confounders",
    "organelle_comparison",
    "population_diversity",
    "genome_coordinate",
)
FORMATS = ("png", "pdf", "svg")


def validate_figure_manifest(rows: Sequence[Mapping[str, str]]) -> None:
    grouped: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        grouped[row["family"]].add(row["format"])
    if tuple(grouped) != FIGURE_FAMILIES and set(grouped) != set(FIGURE_FAMILIES):
        raise ValueError(f"Figure manifest must contain exactly {FIGURE_FAMILIES}")
    for family in FIGURE_FAMILIES:
        if grouped[family] != set(FORMATS):
            raise ValueError(f"Figure family {family} does not have PNG, PDF, and SVG")
