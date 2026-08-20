"""Deterministic commands and interpretation for supplementary analyses."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class KChoice:
    k: int
    cv_error: float
    is_boundary: bool


@dataclass(frozen=True)
class IQTreeSupport:
    sh_alrt: float
    ultrafast_bootstrap: float


@dataclass(frozen=True)
class AlignmentCallabilityCounts:
    coordinate_span_sites: int
    sites_with_any_callable_sample: int
    sites_with_at_least_two_callable_samples: int
    jointly_callable_sites: int


def alignment_callability_counts(records: dict[str, str]) -> AlignmentCallabilityCounts:
    """Distinguish padded coordinate span from sites carrying sequence data."""

    if not records:
        raise ValueError("Alignment callability counts require sequences")
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) != 1:
        raise ValueError("Tree alignment sequences must have equal lengths")
    span = next(iter(lengths))
    callable_counts = [sum(sequence[index].upper() in "ACGT" for sequence in records.values()) for index in range(span)]
    return AlignmentCallabilityCounts(
        coordinate_span_sites=span,
        sites_with_any_callable_sample=sum(count >= 1 for count in callable_counts),
        sites_with_at_least_two_callable_samples=sum(count >= 2 for count in callable_counts),
        jointly_callable_sites=sum(count == len(records) for count in callable_counts),
    )


def alignment_site_counts(records: dict[str, str]) -> tuple[int, int]:
    """Return variable and parsimony-informative alignment-site counts."""

    if not records:
        raise ValueError("Alignment site counts require sequences")
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) != 1:
        raise ValueError("Tree alignment sequences must have equal lengths")
    variable = 0
    parsimony_informative = 0
    for index in range(next(iter(lengths))):
        counts: dict[str, int] = {}
        for sequence in records.values():
            allele = sequence[index].upper()
            if allele in "ACGT":
                counts[allele] = counts.get(allele, 0) + 1
        variable += len(counts) > 1
        parsimony_informative += sum(count >= 2 for count in counts.values()) >= 2
    return variable, parsimony_informative


def parse_iqtree_support(name: str | None, confidence: float | None) -> IQTreeSupport:
    """Parse IQ-TREE's SH-aLRT/UFBoot internal-node label."""

    if name and "/" in name:
        parts = name.split("/")
        if len(parts) == 2:
            return IQTreeSupport(float(parts[0]), float(parts[1]))
    raise ValueError(f"Expected dual SH-aLRT/UFBoot support label, found name={name!r}, confidence={confidence!r}")


def is_strong_iqtree_support(
    support: IQTreeSupport,
    minimum_sh_alrt: float = 80.0,
    minimum_ultrafast_bootstrap: float = 95.0,
) -> bool:
    """Use conventional joint support cutoffs for conflict screening."""

    return support.sh_alrt >= minimum_sh_alrt and support.ultrafast_bootstrap >= minimum_ultrafast_bootstrap


def select_best_k(cv_errors: dict[int, float], tested_min: int, tested_max: int) -> KChoice:
    if not cv_errors:
        raise ValueError("At least one cross-validation result is required")
    k = min(cv_errors, key=lambda value: (cv_errors[value], value))
    return KChoice(k=k, cv_error=cv_errors[k], is_boundary=k in {tested_min, tested_max})


def build_iqtree_command(
    alignment: str | Path,
    prefix: str | Path,
    seed: int,
    partition_file: str | Path | None = None,
    model: str = "MFP",
    sh_alrt_replicates: int = 1000,
    ultrafast_bootstrap_replicates: int = 1000,
    bnni: bool = True,
) -> str:
    partition = f" -p {partition_file}" if partition_file is not None else ""
    bnni_flag = " -bnni" if bnni else ""
    return (
        f"iqtree3 -s {alignment}{partition} -pre {prefix} -st DNA -m {model} "
        f"-alrt {sh_alrt_replicates} -B {ultrafast_bootstrap_replicates}{bnni_flag} -seed {seed}"
    )
