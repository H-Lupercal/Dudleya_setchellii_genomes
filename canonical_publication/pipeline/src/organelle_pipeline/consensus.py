"""Construct full callable-site haploid alignments from all-site calls."""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class ReferenceConcordance:
    callable_bases: int
    reference_matches: int
    nonreference_bases: int
    identity: float


def analysis_mask_length(intervals: list[tuple[int, int]], reference_length: int) -> int:
    """Return the union length of validated 0-based half-open mask intervals."""

    if reference_length < 1:
        raise ValueError("reference length must be positive")
    if not intervals:
        raise ValueError("analysis mask contains no intervals")
    for start, end in intervals:
        if start < 0 or end <= start or end > reference_length:
            raise ValueError(f"Analysis-mask interval {start}:{end} lies outside reference")
    merged_length = 0
    current_start, current_end = sorted(intervals)[0]
    for start, end in sorted(intervals)[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged_length += current_end - current_start
            current_start, current_end = start, end
    return merged_length + current_end - current_start


def build_callable_sequence(
    reference: str,
    rows: list[tuple[int, str, str, str]],
    passing_variant_positions: set[int],
    invariant_rows: list[tuple[int, str, str, int, tuple[int, ...]]] | None = None,
    minimum_depth: int = 5,
    minimum_genotype_quality: int = 20,
) -> str:
    """Keep confident invariant/accepted variant calls; represent all else as N.

    Bcftools' multiallelic caller does not emit GQ at invariant sites. Reference
    confidence is therefore calculated from the raw diploid-layout mpileup PL
    vector. For haploid data, only homozygous likelihoods are relevant: index 0
    for reference and indices k*(k+3)/2 for each alternative allele k. All
    heterozygous likelihoods are ignored.
    """

    sequence = ["N"] * len(reference)
    for position, ref, alt, depth, likelihoods in invariant_rows or []:
        if position < 1 or position > len(reference):
            raise ValueError(f"Position {position} lies outside reference")
        ref = ref.upper()
        start = position - 1
        end = start + len(ref)
        if end > len(reference) or reference[start:end].upper() != ref:
            raise ValueError(f"Reference mismatch at position {position}")
        alternatives = [] if alt == "." else alt.split(",")
        allele_count = 1 + len(alternatives)
        expected_likelihood_count = allele_count * (allele_count + 1) // 2
        homozygous_alt_indices = [allele * (allele + 3) // 2 for allele in range(1, allele_count)]
        if len(ref) != 1 or depth < minimum_depth or len(likelihoods) != expected_likelihood_count or not homozygous_alt_indices:
            continue
        reference_likelihood = likelihoods[0]
        best_alt_likelihood = min(likelihoods[index] for index in homozygous_alt_indices)
        if reference_likelihood < best_alt_likelihood and best_alt_likelihood - reference_likelihood >= minimum_genotype_quality:
            sequence[start] = ref
    for position, ref, alt, genotype in rows:
        if position < 1 or position > len(reference):
            raise ValueError(f"Position {position} lies outside reference")
        ref = ref.upper()
        start = position - 1
        end = start + len(ref)
        if end > len(reference) or reference[start:end].upper() != ref:
            raise ValueError(f"Reference mismatch at position {position}")
        if genotype in {".", "./.", ".|."}:
            continue
        if alt in {".", "<*>", "<NON_REF>"}:
            continue
        if position not in passing_variant_positions:
            continue
        alleles = [ref, *(value.upper() for value in alt.split(","))]
        try:
            allele_index = int(genotype)
        except ValueError:
            continue
        if allele_index >= len(alleles) or len(alleles[allele_index]) != 1:
            continue
        sequence[position - 1] = alleles[allele_index]
    return "".join(sequence)


def reference_concordance(reference: str, sequence: str) -> ReferenceConcordance:
    """Summarize read-backed callable consensus agreement with its reference."""

    if len(reference) != len(sequence):
        raise ValueError("Reference and callable sequence must have equal lengths")
    pairs = [(ref.upper(), observed.upper()) for ref, observed in zip(reference, sequence, strict=True) if observed.upper() in "ACGT"]
    matches = sum(reference_base == observed for reference_base, observed in pairs)
    callable_bases = len(pairs)
    return ReferenceConcordance(
        callable_bases=callable_bases,
        reference_matches=matches,
        nonreference_bases=callable_bases - matches,
        identity=matches / callable_bases if callable_bases else math.nan,
    )
