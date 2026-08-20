"""Reference-alignment metrics and inverted-repeat detection."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class BlastHit:
    query: str
    subject: str
    identity_percent: float
    alignment_length: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    bitscore: float
    query_length: int
    subject_length: int


@dataclass(frozen=True)
class BlastSummary:
    query_covered_bp: int
    subject_covered_bp: int
    query_coverage: float
    subject_coverage: float
    weighted_identity_percent: float
    hit_count: int


@dataclass(frozen=True)
class OrientationSummary:
    same_orientation_query_bp: int
    reverse_complement_query_bp: int
    same_orientation_fraction_of_covered: float
    reverse_complement_fraction_of_covered: float


def read_blast_hits(lines: list[str] | tuple[str, ...]) -> list[BlastHit]:
    hits: list[BlastHit] = []
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip().split("\t")
        if len(fields) != 11:
            raise ValueError(f"Expected 11 BLAST columns, found {len(fields)}")
        hits.append(
            BlastHit(
                query=fields[0],
                subject=fields[1],
                identity_percent=float(fields[2]),
                alignment_length=int(fields[3]),
                query_start=int(fields[4]),
                query_end=int(fields[5]),
                subject_start=int(fields[6]),
                subject_end=int(fields[7]),
                bitscore=float(fields[8]),
                query_length=int(fields[9]),
                subject_length=int(fields[10]),
            )
        )
    return hits


def _union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    merged = 0
    current_start, current_end = sorted(intervals)[0]
    for start, end in sorted(intervals)[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            merged += current_end - current_start + 1
            current_start, current_end = start, end
    return merged + current_end - current_start + 1


def summarize_blast_hits(hits: list[BlastHit]) -> BlastSummary:
    if not hits:
        raise ValueError("Cannot summarize zero BLAST hits")
    query_covered = _union_length([(min(hit.query_start, hit.query_end), max(hit.query_start, hit.query_end)) for hit in hits])
    subject_covered = _union_length([(min(hit.subject_start, hit.subject_end), max(hit.subject_start, hit.subject_end)) for hit in hits])
    # Attribute each query position to only its highest-bitscore HSP. Repeats
    # and fragmented secondary hits otherwise inflate and bias the identity.
    identity_by_query_position: list[float | None] = [None] * hits[0].query_length
    for hit in sorted(hits, key=lambda value: value.bitscore, reverse=True):
        start = min(hit.query_start, hit.query_end) - 1
        end = max(hit.query_start, hit.query_end)
        for index in range(start, end):
            if identity_by_query_position[index] is None:
                identity_by_query_position[index] = hit.identity_percent
    assigned_identity = [value for value in identity_by_query_position if value is not None]
    weighted_identity = sum(assigned_identity) / len(assigned_identity)
    return BlastSummary(
        query_covered_bp=query_covered,
        subject_covered_bp=subject_covered,
        query_coverage=query_covered / hits[0].query_length,
        subject_coverage=subject_covered / hits[0].subject_length,
        weighted_identity_percent=weighted_identity,
        hit_count=len(hits),
    )


def summarize_query_orientation(hits: list[BlastHit]) -> OrientationSummary:
    """Assign each covered query base to its highest-bitscore HSP orientation."""

    if not hits:
        raise ValueError("Cannot summarize orientation for zero BLAST hits")
    assigned: list[bool | None] = [None] * hits[0].query_length
    for hit in sorted(hits, key=lambda value: value.bitscore, reverse=True):
        same_orientation = (hit.query_end - hit.query_start) * (hit.subject_end - hit.subject_start) > 0
        start = min(hit.query_start, hit.query_end) - 1
        end = max(hit.query_start, hit.query_end)
        for index in range(start, end):
            if assigned[index] is None:
                assigned[index] = same_orientation
    same = sum(value is True for value in assigned)
    reverse = sum(value is False for value in assigned)
    covered = same + reverse
    return OrientationSummary(
        same_orientation_query_bp=same,
        reverse_complement_query_bp=reverse,
        same_orientation_fraction_of_covered=same / covered,
        reverse_complement_fraction_of_covered=reverse / covered,
    )


def select_inverted_repeat_pair(
    hits: list[BlastHit], minimum_length: int = 10_000, minimum_identity: float = 99.0
) -> tuple[tuple[int, int], tuple[int, int]]:
    candidates = [
        hit
        for hit in hits
        if hit.query == hit.subject
        and hit.subject_start > hit.subject_end
        and hit.alignment_length >= minimum_length
        and hit.identity_percent >= minimum_identity
    ]
    if not candidates:
        raise ValueError("No qualifying inverted-repeat pair found")
    hit = max(candidates, key=lambda value: (value.alignment_length, value.bitscore))
    query_interval = (min(hit.query_start, hit.query_end) - 1, max(hit.query_start, hit.query_end))
    subject_interval = (
        min(hit.subject_start, hit.subject_end) - 1,
        max(hit.subject_start, hit.subject_end),
    )
    first, second = sorted((query_interval, subject_interval))
    return first, second


def select_terminal_direct_repeat(
    hits: list[BlastHit],
    retained_length: int,
    sequence_length: int,
    minimum_identity: float = 99.0,
    boundary_tolerance: int = 25,
) -> tuple[tuple[int, int], tuple[int, int]]:
    """Validate a configured circular trim using a start-to-end direct repeat."""

    candidates: list[tuple[BlastHit, tuple[int, int], tuple[int, int]]] = []
    for hit in hits:
        if hit.query != hit.subject or hit.identity_percent < minimum_identity:
            continue
        if (hit.query_end - hit.query_start) * (hit.subject_end - hit.subject_start) <= 0:
            continue
        first, second = sorted(
            (
                (min(hit.query_start, hit.query_end) - 1, max(hit.query_start, hit.query_end)),
                (min(hit.subject_start, hit.subject_end) - 1, max(hit.subject_start, hit.subject_end)),
            )
        )
        if (
            first[0] <= boundary_tolerance
            and abs(second[0] - retained_length) <= boundary_tolerance
            and abs(second[1] - sequence_length) <= boundary_tolerance
        ):
            candidates.append((hit, first, second))
    if not candidates:
        raise ValueError("No high-identity start-to-end repeat validates the chloroplast trim boundary")
    _, first, second = max(candidates, key=lambda value: (value[0].alignment_length, value[0].bitscore))
    return first, second


def self_repeat_intervals(hits: list[BlastHit], minimum_length: int = 100, minimum_identity: float = 95.0) -> tuple[tuple[int, int], ...]:
    """Merge both copies of qualifying non-diagonal self-alignment hits."""

    intervals: list[tuple[int, int]] = []
    for hit in hits:
        query = (min(hit.query_start, hit.query_end) - 1, max(hit.query_start, hit.query_end))
        subject = (
            min(hit.subject_start, hit.subject_end) - 1,
            max(hit.subject_start, hit.subject_end),
        )
        if hit.query != hit.subject or hit.alignment_length < minimum_length or hit.identity_percent < minimum_identity or query == subject:
            continue
        intervals.extend((query, subject))
    if not intervals:
        return ()
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return tuple(merged)
