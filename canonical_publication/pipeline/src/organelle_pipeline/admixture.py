"""Explicit pseudo-diploid encoding used only for supplementary ADMIXTURE."""

from __future__ import annotations

import math


def pseudo_diploid_alleles(genotype: str, ref: str, alt: str) -> tuple[str, str]:
    if genotype == "0":
        return ref, ref
    if genotype == "1" and len(alt) == 1:
        return alt, alt
    return "0", "0"


def validate_q_matrix(lines: list[str], sample_count: int, k: int, tolerance: float = 1e-5) -> None:
    """Reject truncated, malformed, or non-probability ADMIXTURE Q output."""

    if len(lines) != sample_count:
        raise ValueError(f"ADMIXTURE Q row count {len(lines)} != sample count {sample_count}")
    for row_number, line in enumerate(lines, 1):
        values = [float(value) for value in line.split()]
        if len(values) != k:
            raise ValueError(f"ADMIXTURE Q row {row_number} has {len(values)} columns, expected K={k}")
        if any(not math.isfinite(value) or value < 0 or value > 1 for value in values):
            raise ValueError(f"ADMIXTURE Q row {row_number} contains an invalid probability")
        if not math.isclose(sum(values), 1.0, rel_tol=0, abs_tol=tolerance):
            raise ValueError(f"ADMIXTURE Q row {row_number} probabilities do not sum to one")
