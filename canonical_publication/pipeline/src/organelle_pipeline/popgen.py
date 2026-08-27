"""Callable-site diversity and Hudson FST for haploid organelle data."""

from __future__ import annotations

import math
import random
from collections import Counter
from dataclasses import dataclass
from itertools import combinations

import numpy as np

BASES = frozenset("ACGT")


@dataclass(frozen=True)
class DiversityResult:
    differences: int
    compared_sites: int
    jointly_callable_sites: int
    pi: float


@dataclass(frozen=True)
class HudsonResult:
    numerator: float
    denominator: float
    callable_sites: int
    fst: float


@dataclass(frozen=True)
class HaplotypeDiversityResult:
    haplotype_count: int
    diversity: float
    assigned_samples: int
    ambiguous_samples: int


@dataclass(frozen=True)
class PairwiseDistanceResult:
    differences: int
    sites_compared: int
    p_distance: float


@dataclass(frozen=True)
class PackedSequence:
    length: int
    base_masks: tuple[int, int, int, int]
    callable_mask: int


def pack_sequence(sequence: str) -> PackedSequence:
    """Pack A/C/G/T positions into bitsets for fast all-pairs comparisons."""

    encoded = np.frombuffer(sequence.upper().encode("ascii"), dtype=np.uint8)
    masks = tuple(
        int.from_bytes(np.packbits(encoded == ord(base), bitorder="little").tobytes(), "little")
        for base in "ACGT"
    )
    return PackedSequence(
        length=len(sequence),
        base_masks=masks,  # type: ignore[arg-type]
        callable_mask=masks[0] | masks[1] | masks[2] | masks[3],
    )


def packed_pairwise_distance(left: PackedSequence, right: PackedSequence) -> PairwiseDistanceResult:
    """Count substitutions at positions callable in both packed sequences."""

    if left.length != right.length:
        raise ValueError("Pairwise-distance sequences must have equal lengths")
    jointly_callable = left.callable_mask & right.callable_mask
    sites_compared = jointly_callable.bit_count()
    matches = sum((left_mask & right_mask).bit_count() for left_mask, right_mask in zip(left.base_masks, right.base_masks, strict=True))
    differences = sites_compared - matches
    return PairwiseDistanceResult(
        differences=differences,
        sites_compared=sites_compared,
        p_distance=(differences / sites_compared) if sites_compared else math.nan,
    )


def pairwise_sequence_distance(left: str, right: str) -> PairwiseDistanceResult:
    """Count raw nucleotide differences, excluding non-ACGT calls pairwise."""

    if len(left) != len(right):
        raise ValueError("Pairwise-distance sequences must have equal lengths")
    return packed_pairwise_distance(pack_sequence(left), pack_sequence(right))


def haplotype_diversity_from_assignments(
    assignments: tuple[str, ...] | list[str],
    ambiguous_label: str = "AMBIGUOUS",
) -> HaplotypeDiversityResult:
    """Calculate unbiased diversity from organelle-wide haplotype labels."""

    counts = Counter(label for label in assignments if label != ambiguous_label)
    assigned = sum(counts.values())
    if assigned > 1:
        diversity = (assigned / (assigned - 1)) * (1 - sum((count / assigned) ** 2 for count in counts.values()))
    elif assigned == 1:
        diversity = 0.0
    else:
        diversity = math.nan
    return HaplotypeDiversityResult(
        haplotype_count=len(counts),
        diversity=diversity,
        assigned_samples=assigned,
        ambiguous_samples=len(assignments) - assigned,
    )


def _equal_sequence_length(*populations: tuple[str, ...] | list[str]) -> int:
    lengths = {len(sequence) for population in populations for sequence in population}
    if len(lengths) != 1:
        raise ValueError("Hudson FST sequences must have equal lengths")
    return next(iter(lengths))


def callable_nucleotide_diversity(sequences: tuple[str, ...] | list[str]) -> DiversityResult:
    """Mean pairwise differences across positions callable in every sample."""

    if not sequences:
        return DiversityResult(0, 0, 0, math.nan)
    lengths = {len(sequence) for sequence in sequences}
    if len(lengths) != 1:
        raise ValueError("Nucleotide-diversity sequences must have equal lengths")
    joint_positions = [index for index in range(next(iter(lengths))) if all(sequence[index].upper() in BASES for sequence in sequences)]
    differences = 0
    for left, right in combinations(sequences, 2):
        differences += sum(left[index].upper() != right[index].upper() for index in joint_positions)
    compared = len(joint_positions) * (len(sequences) * (len(sequences) - 1) // 2)
    return DiversityResult(
        differences=differences,
        compared_sites=compared,
        jointly_callable_sites=len(joint_positions),
        pi=(differences / compared) if compared else math.nan,
    )


def _site_hudson(alleles1: list[str], alleles2: list[str]) -> tuple[float, float] | None:
    if len(alleles1) < 2 or len(alleles2) < 2:
        return None
    observed = sorted(set(alleles1) | set(alleles2))
    if len(observed) != 2:
        return None
    alternate = observed[1]
    n1, n2 = len(alleles1), len(alleles2)
    p1 = alleles1.count(alternate) / n1
    p2 = alleles2.count(alternate) / n2
    numerator = (p1 - p2) ** 2 - (p1 * (1 - p1) / (n1 - 1)) - (p2 * (1 - p2) / (n2 - 1))
    denominator = p1 * (1 - p2) + p2 * (1 - p1)
    if denominator == 0:
        return None
    return numerator, denominator


def hudson_fst(
    population1: tuple[str, ...] | list[str],
    population2: tuple[str, ...] | list[str],
) -> HudsonResult:
    """Calculate the ratio-of-sums Hudson FST without zero clamping."""

    if not population1 or not population2:
        return HudsonResult(0.0, 0.0, 0, math.nan)
    site_count = _equal_sequence_length(population1, population2)
    numerator = 0.0
    denominator = 0.0
    callable_sites = 0
    for index in range(site_count):
        alleles1 = [seq[index].upper() for seq in population1 if seq[index].upper() in BASES]
        alleles2 = [seq[index].upper() for seq in population2 if seq[index].upper() in BASES]
        if len(alleles1) >= 2 and len(alleles2) >= 2:
            callable_sites += 1
        components = _site_hudson(alleles1, alleles2)
        if components is None:
            continue
        site_numerator, site_denominator = components
        numerator += site_numerator
        denominator += site_denominator
    return HudsonResult(
        numerator=numerator,
        denominator=denominator,
        callable_sites=callable_sites,
        fst=(numerator / denominator) if denominator else math.nan,
    )


def _hudson_blocks(
    population1: tuple[str, ...] | list[str],
    population2: tuple[str, ...] | list[str],
    block_size: int,
) -> tuple[tuple[float, float], ...]:
    """Sum Hudson components within every physically callable sequence block."""

    if not population1 or not population2:
        return ()
    site_count = _equal_sequence_length(population1, population2)
    blocks: list[tuple[float, float]] = []
    for start in range(0, site_count, block_size):
        numerator = 0.0
        denominator = 0.0
        callable_block = False
        for index in range(start, min(start + block_size, site_count)):
            alleles1 = [seq[index].upper() for seq in population1 if seq[index].upper() in BASES]
            alleles2 = [seq[index].upper() for seq in population2 if seq[index].upper() in BASES]
            if len(alleles1) >= 2 and len(alleles2) >= 2:
                callable_block = True
            components = _site_hudson(alleles1, alleles2)
            if components is not None:
                numerator += components[0]
                denominator += components[1]
        if callable_block:
            blocks.append((numerator, denominator))
    return tuple(blocks)


def block_bootstrap_hudson_fst(
    population1: tuple[str, ...] | list[str],
    population2: tuple[str, ...] | list[str],
    block_size: int = 1000,
    replicates: int = 1000,
    seed: int = 1729,
) -> tuple[float, float]:
    """Return deterministic percentile confidence limits from site blocks."""

    if block_size <= 0 or replicates <= 0:
        raise ValueError("block_size and replicates must be positive")
    blocks = _hudson_blocks(population1, population2, block_size)
    if not blocks:
        return math.nan, math.nan
    rng = random.Random(seed)
    estimates: list[float] = []
    for _ in range(replicates):
        chosen = [rng.choice(blocks) for _ in blocks]
        numerator = sum(block[0] for block in chosen)
        denominator = sum(block[1] for block in chosen)
        estimate = numerator / denominator if denominator else math.nan
        if math.isfinite(estimate):
            estimates.append(estimate)
    if not estimates:
        return math.nan, math.nan
    estimates.sort()
    lower = estimates[int(0.025 * (len(estimates) - 1))]
    upper = estimates[int(0.975 * (len(estimates) - 1))]
    return lower, upper


def private_variant_sites(
    population: str,
    groups: dict[str, tuple[str, ...] | list[str]],
    reference: str,
) -> tuple[int, ...]:
    """Return zero-based sites with any allele private to a population."""

    return private_variant_sites_all(groups, reference)[population]


def private_variant_sites_all(
    groups: dict[str, tuple[str, ...] | list[str]],
    reference: str,
    require_joint_callability: bool = False,
) -> dict[str, tuple[int, ...]]:
    """Find private non-reference alleles for all populations in one pass."""

    all_sequences = [seq for values in groups.values() for seq in values]
    result: dict[str, list[int]] = {population: [] for population in groups}
    if not all_sequences:
        return {population: () for population in groups}
    sequence_lengths = {len(sequence) for sequence in all_sequences}
    if len(sequence_lengths) != 1:
        raise ValueError("Private-variant sequences must have equal lengths")
    site_count = next(iter(sequence_lengths))
    if len(reference) < site_count:
        raise ValueError("Reference is shorter than population sequences")
    for index in range(site_count):
        if require_joint_callability and any(sequence[index].upper() not in BASES for sequence in all_sequences):
            continue
        allele_populations: dict[str, set[str]] = {}
        for population, sequences in groups.items():
            for sequence in sequences:
                allele = sequence[index].upper()
                if allele in BASES and allele != reference[index].upper():
                    allele_populations.setdefault(allele, set()).add(population)
        private_populations = {next(iter(populations)) for populations in allele_populations.values() if len(populations) == 1}
        for population in private_populations:
            result[population].append(index)
    return {population: tuple(sites) for population, sites in result.items()}
