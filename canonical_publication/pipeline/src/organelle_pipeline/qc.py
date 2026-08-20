"""Coverage-based organelle sample eligibility."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class SampleBreadth:
    sample_id: str
    cp_dp5: float
    mt_dp5: float


@dataclass(frozen=True)
class OrganelleSampleSets:
    cp: tuple[str, ...]
    mt: tuple[str, ...]
    shared: tuple[str, ...]


@dataclass(frozen=True)
class DepthSummary:
    reference_length: int
    mean_depth: float
    breadth: dict[int, float]


@dataclass(frozen=True)
class ReadPreprocessingSummary:
    input_reads: int
    input_bases: int
    input_q20_rate: float
    input_q30_rate: float
    passing_reads: int
    passing_bases: int
    passing_q20_rate: float
    passing_q30_rate: float
    read_retention: float
    low_quality_reads: int
    too_many_n_reads: int
    adapter_dimer_reads: int
    too_short_reads: int
    adapter_trimmed_reads: int
    adapter_trimmed_bases: int
    duplication_rate: float


def summarize_fastp_report(payload: Mapping[str, object]) -> ReadPreprocessingSummary:
    """Validate and extract acceptance-relevant metrics from fastp JSON."""

    summary = payload.get("summary")
    filtering = payload.get("filtering_result")
    adapter = payload.get("adapter_cutting")
    duplication = payload.get("duplication")
    if not all(isinstance(value, Mapping) for value in (summary, filtering, adapter, duplication)):
        raise ValueError("fastp report is missing a required metrics section")
    assert isinstance(summary, Mapping)
    assert isinstance(filtering, Mapping)
    assert isinstance(adapter, Mapping)
    assert isinstance(duplication, Mapping)
    before = summary.get("before_filtering")
    after = summary.get("after_filtering")
    if not isinstance(before, Mapping) or not isinstance(after, Mapping):
        raise ValueError("fastp report is missing before/after filtering summaries")

    input_reads = int(before["total_reads"])
    passing_reads = int(after["total_reads"])
    reported_passing = int(filtering["passed_filter_reads"])
    if input_reads <= 0 or passing_reads < 0 or passing_reads > input_reads:
        raise ValueError("fastp report has impossible input/passing read counts")
    if reported_passing != passing_reads:
        raise ValueError("fastp passed-filter count disagrees with its after-filtering summary")
    rate_fields = (
        float(before["q20_rate"]),
        float(before["q30_rate"]),
        float(after["q20_rate"]),
        float(after["q30_rate"]),
        float(duplication["rate"]),
    )
    if any(not 0 <= value <= 1 for value in rate_fields):
        raise ValueError("fastp report contains a rate outside [0, 1]")
    adapter_trimmed_reads = int(adapter["adapter_trimmed_reads"])
    if not 0 <= adapter_trimmed_reads <= input_reads:
        raise ValueError("fastp adapter-trimmed read count is impossible")
    return ReadPreprocessingSummary(
        input_reads=input_reads,
        input_bases=int(before["total_bases"]),
        input_q20_rate=rate_fields[0],
        input_q30_rate=rate_fields[1],
        passing_reads=passing_reads,
        passing_bases=int(after["total_bases"]),
        passing_q20_rate=rate_fields[2],
        passing_q30_rate=rate_fields[3],
        read_retention=passing_reads / input_reads,
        low_quality_reads=int(filtering["low_quality_reads"]),
        too_many_n_reads=int(filtering["too_many_N_reads"]),
        adapter_dimer_reads=int(filtering["adapter_dimer_reads"]),
        too_short_reads=int(filtering["too_short_reads"]),
        adapter_trimmed_reads=adapter_trimmed_reads,
        adapter_trimmed_bases=int(adapter["adapter_trimmed_bases"]),
        duplication_rate=rate_fields[4],
    )


def build_depth_command(
    bam: Path | str,
    record: str,
    minimum_mapping_quality: int,
    minimum_base_quality: int,
) -> list[str]:
    """Build `samtools depth` arguments with its BQ/MQ option semantics."""

    return [
        "samtools",
        "depth",
        "-aa",
        "-s",
        "-q",
        str(minimum_base_quality),
        "-Q",
        str(minimum_mapping_quality),
        "-r",
        record,
        str(bam),
    ]


def select_organelle_samples(
    rows: list[SampleBreadth] | tuple[SampleBreadth, ...],
    minimum_dp5_breadth: float = 0.80,
) -> OrganelleSampleSets:
    """Select cpDNA and mtDNA samples independently using an inclusive cutoff."""

    if not 0 <= minimum_dp5_breadth <= 1:
        raise ValueError("minimum_dp5_breadth must be between zero and one")
    cp = tuple(row.sample_id for row in rows if row.cp_dp5 >= minimum_dp5_breadth)
    mt = tuple(row.sample_id for row in rows if row.mt_dp5 >= minimum_dp5_breadth)
    mt_set = set(mt)
    shared = tuple(sample_id for sample_id in cp if sample_id in mt_set)
    return OrganelleSampleSets(cp=cp, mt=mt, shared=shared)


def summarize_depths(depths: list[int] | tuple[int, ...], thresholds: tuple[int, ...] = (1, 3, 5, 10)) -> DepthSummary:
    if not depths:
        raise ValueError("depth vector cannot be empty")
    if any(threshold <= 0 for threshold in thresholds):
        raise ValueError("depth thresholds must be positive")
    length = len(depths)
    return DepthSummary(
        reference_length=length,
        mean_depth=sum(depths) / length,
        breadth={threshold: sum(depth >= threshold for depth in depths) / length for threshold in thresholds},
    )


def summarize_masked_depths(
    depths: list[int] | tuple[int, ...],
    included: list[bool] | tuple[bool, ...],
    thresholds: tuple[int, ...] = (1, 3, 5, 10),
) -> DepthSummary:
    """Summarize breadth over explicitly mappable analysis coordinates."""

    if len(depths) != len(included):
        raise ValueError("depth and analysis-mask vectors must have equal lengths")
    selected = [depth for depth, keep in zip(depths, included, strict=True) if keep]
    if not selected:
        raise ValueError("analysis mask contains no included positions")
    return summarize_depths(selected, thresholds)


def support_intervals(
    support_counts: list[int] | tuple[int, ...],
    sample_count: int,
    minimum_fraction: float,
    minimum_length: int,
) -> tuple[tuple[int, int], ...]:
    """Convert per-position sample support to 0-based half-open intervals."""

    if sample_count <= 0 or not 0 <= minimum_fraction <= 1 or minimum_length <= 0:
        raise ValueError("invalid support interval parameters")
    intervals: list[tuple[int, int]] = []
    start: int | None = None
    for index, support in enumerate(support_counts):
        passing = support / sample_count >= minimum_fraction
        if passing and start is None:
            start = index
        elif not passing and start is not None:
            if index - start >= minimum_length:
                intervals.append((start, index))
            start = None
    if start is not None and len(support_counts) - start >= minimum_length:
        intervals.append((start, len(support_counts)))
    return tuple(intervals)
