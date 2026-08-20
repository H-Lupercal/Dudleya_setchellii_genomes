"""Utilities for evidence-backed external feature projection."""

from __future__ import annotations

from .reference_evidence import BlastHit


def gff_phase(feature_type: str, codon_start: str | None) -> str:
    """Convert GenBank's one-based CDS codon_start to GFF3 phase."""

    if feature_type != "CDS":
        return "."
    value = int(codon_start or "1")
    if value not in {1, 2, 3}:
        raise ValueError(f"Invalid GenBank codon_start: {codon_start!r}")
    return str(value - 1)


def projected_interval(hit: BlastHit) -> tuple[int, int, str]:
    """Return a 0-based half-open target interval and projected strand."""

    return (
        min(hit.subject_start, hit.subject_end) - 1,
        max(hit.subject_start, hit.subject_end),
        "+" if hit.subject_start <= hit.subject_end else "-",
    )
