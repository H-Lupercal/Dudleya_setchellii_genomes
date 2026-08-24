"""Conservative sample-identity and mixed-allele screening."""

from __future__ import annotations

import re
import statistics
from collections.abc import Mapping, Sequence
from dataclasses import dataclass


@dataclass(frozen=True)
class MixedAlleleCall:
    depth: int
    reference_depth: int
    alternate_depth: int


@dataclass(frozen=True)
class MixedAlleleSummary:
    sample_id: str
    evaluated_site_count: int
    mixed_site_count: int
    mixed_site_fraction: float
    status: str


def parse_structured_id(filename: str) -> dict[str, str]:
    patterns = {
        "plate_well": r"(LP_\d+)",
        "specimen": r"(Du-\d+[A-Za-z]?)",
        "demultiplex_sample": r"_(S\d+)_",
        "lane": r"_(L\d{3})_",
    }
    return {key: match.group(1) for key, pattern in patterns.items() if (match := re.search(pattern, filename, re.IGNORECASE))}


def index_hopping_status(index_sequences: Sequence[str], demultiplex_metrics: Sequence[Mapping[str, object]]) -> str:
    return "testable" if index_sequences and demultiplex_metrics else "untestable"


def _is_mixed(call: MixedAlleleCall) -> bool:
    if call.depth < 20 or call.reference_depth < 3 or call.alternate_depth < 3:
        return False
    total = call.reference_depth + call.alternate_depth
    fraction = min(call.reference_depth, call.alternate_depth) / total
    return 0.20 <= fraction <= 0.80


def classify_mixed_allele_samples(calls: Mapping[str, Sequence[MixedAlleleCall]]) -> list[MixedAlleleSummary]:
    raw: list[tuple[str, int, int, float]] = []
    for sample, values in calls.items():
        evaluated = sum(value.depth >= 20 for value in values)
        mixed = sum(_is_mixed(value) for value in values)
        raw.append((sample, evaluated, mixed, mixed / evaluated if evaluated else 0.0))
    rates = [value[3] for value in raw]
    median = statistics.median(rates) if rates else 0.0
    mad = statistics.median(abs(value - median) for value in rates) if rates else 0.0
    threshold = median + 5 * mad
    return [
        MixedAlleleSummary(sample, evaluated, mixed, rate, "suspected" if mixed >= 10 and rate > threshold else "not_flagged")
        for sample, evaluated, mixed, rate in raw
    ]
