"""Organelle haplotypes retaining all accepted variable positions."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass


@dataclass(frozen=True)
class HaplotypeSummary:
    positions: tuple[int, ...]
    sample_haplotypes: dict[str, str]
    sequences: dict[str, str]
    counts: dict[str, int]


def summarize_haplotypes(records: dict[str, str]) -> HaplotypeSummary:
    if not records:
        raise ValueError("haplotype summary requires sequences")
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) != 1:
        raise ValueError("haplotype alignment sequences must have equal lengths")
    length = next(iter(lengths))
    positions = tuple(
        index for index in range(length) if len({sequence[index] for sequence in records.values() if sequence[index] in "ACGT"}) > 1
    )
    sample_strings = {sample: "".join(sequence[index] for index in positions) for sample, sequence in records.items()}
    unique = sorted({sequence for sequence in sample_strings.values() if set(sequence) <= set("ACGT")})
    labels = {sequence: f"H{index}" for index, sequence in enumerate(unique, 1)}
    sample_haplotypes = {sample: labels[sequence] if sequence in labels else "AMBIGUOUS" for sample, sequence in sample_strings.items()}
    counts = Counter(label for label in sample_haplotypes.values() if label != "AMBIGUOUS")
    return HaplotypeSummary(
        positions=positions,
        sample_haplotypes=sample_haplotypes,
        sequences={label: sequence for sequence, label in labels.items()},
        counts=dict(sorted(counts.items())),
    )
