from organelle_pipeline.reference_evidence import (
    BlastHit,
    select_inverted_repeat_pair,
    select_terminal_direct_repeat,
    self_repeat_intervals,
    summarize_blast_hits,
    summarize_query_orientation,
)


def test_blast_summary_uses_interval_unions_for_coverage() -> None:
    hits = [
        BlastHit("q", "s", 99.0, 60, 1, 60, 1, 60, 100.0, 100, 100),
        BlastHit("q", "s", 98.0, 51, 50, 100, 50, 100, 90.0, 100, 100),
    ]
    summary = summarize_blast_hits(hits)
    assert summary.query_covered_bp == 100
    assert summary.subject_covered_bp == 100
    assert summary.query_coverage == 1.0
    assert summary.subject_coverage == 1.0


def test_blast_identity_does_not_double_count_overlapping_secondary_hits() -> None:
    hits = [
        BlastHit("q", "s", 99.0, 60, 1, 60, 1, 60, 100.0, 100, 100),
        BlastHit("q", "s", 0.0, 51, 50, 100, 50, 100, 90.0, 100, 100),
    ]
    assert summarize_blast_hits(hits).weighted_identity_percent == 59.4


def test_orientation_summary_assigns_overlaps_to_highest_bitscore_hsp() -> None:
    hits = [
        BlastHit("q", "s", 99.0, 60, 1, 60, 1, 60, 100.0, 100, 100),
        BlastHit("q", "s", 98.0, 51, 50, 100, 100, 50, 90.0, 100, 100),
    ]
    summary = summarize_query_orientation(hits)
    assert summary.same_orientation_query_bp == 60
    assert summary.reverse_complement_query_bp == 40
    assert summary.same_orientation_fraction_of_covered == 0.6
    assert summary.reverse_complement_fraction_of_covered == 0.4


def test_longest_reciprocal_reverse_hits_define_ir_copies() -> None:
    hits = [
        BlastHit("cp", "cp", 99.9, 20_000, 10, 20_009, 50_000, 30_001, 1000, 100_000, 100_000),
        BlastHit("cp", "cp", 99.9, 20_000, 30_001, 50_000, 20_009, 10, 1000, 100_000, 100_000),
        BlastHit("cp", "cp", 100.0, 100_000, 1, 100_000, 1, 100_000, 2000, 100_000, 100_000),
    ]
    first, second = select_inverted_repeat_pair(hits, minimum_length=10_000)
    assert first == (9, 20_009)
    assert second == (30_000, 50_000)


def test_self_repeat_mask_excludes_non_diagonal_copies() -> None:
    hits = [
        BlastHit("mt", "mt", 100, 1000, 1, 1000, 1, 1000, 2000, 1000, 1000),
        BlastHit("mt", "mt", 99, 100, 10, 109, 500, 599, 200, 1000, 1000),
    ]
    assert self_repeat_intervals(hits, minimum_length=100) == ((9, 109), (499, 599))


def test_terminal_trim_boundary_requires_a_high_identity_end_to_start_repeat() -> None:
    hits = [
        BlastHit("cp", "cp", 99.9, 100, 1, 100, 901, 1000, 500, 1000, 1000),
        BlastHit("cp", "cp", 99.9, 100, 901, 1000, 1, 100, 500, 1000, 1000),
    ]
    assert select_terminal_direct_repeat(hits, retained_length=900, sequence_length=1000) == (
        (0, 100),
        (900, 1000),
    )
