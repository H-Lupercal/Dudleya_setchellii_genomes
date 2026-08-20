import pytest
from organelle_pipeline.consensus import analysis_mask_length, build_callable_sequence, reference_concordance


def test_callable_sequence_keeps_callable_invariants_and_masks_failed_variants() -> None:
    rows = [(2, "C", "T", "1"), (3, "G", "A", "1")]
    invariant_rows = [
        (1, "A", "<*>", 5, (0, 15, 121)),
        (4, "T", "<*>", 4, (0, 12, 99)),
    ]
    sequence = build_callable_sequence(
        "ACGT",
        rows,
        passing_variant_positions={2},
        invariant_rows=invariant_rows,
    )
    assert sequence == "ATNN"


def test_invariant_haploid_gq_uses_homozygous_pls_and_requires_reference_best() -> None:
    invariant_rows = [
        (1, "A", "<*>", 5, (0, 8, 19)),
        (2, "C", "<*>", 5, (0, 8, 20)),
        (3, "G", "<*>", 10, (30, 8, 0)),
    ]

    assert build_callable_sequence("ACG", [], set(), invariant_rows=invariant_rows) == "NCN"


def test_haploid_reference_gq_supports_multiallelic_mpileup_pl_layout() -> None:
    likelihood_rows = [
        # REF=A, ALT=C,<*>; diploid PL order is 0/0, 0/1, 1/1,
        # 0/2, 1/2, 2/2. Haploid comparison uses indices 0, 2, and 5.
        (1, "A", "C,<*>", 9, (0, 27, 229, 27, 229, 229)),
        # The heterozygous PL is deliberately best but is biologically
        # irrelevant for the haploid organelle model.
        (2, "C", "T,A,<*>", 11, (0, 0, 247, 33, 247, 247, 33, 247, 247, 247)),
        # A homozygous alternative is better than reference, so do not call ref.
        (3, "G", "A,<*>", 10, (30, 8, 0, 20, 8, 40)),
    ]

    assert build_callable_sequence("ACG", [], set(), invariant_rows=likelihood_rows) == "ACN"


def test_reference_genotype_at_passing_variant_is_retained() -> None:
    rows = [(1, "A", "T", "0")]
    assert build_callable_sequence("A", rows, {1}) == "A"


def test_overlapping_excluded_indel_validates_full_reference_span_without_becoming_callable() -> None:
    rows = [
        (1, "GAAAAAAAA", "GAAAAAAAAA,GAAAAAAA", "0"),
    ]
    invariant_rows = [
        (1, "G", "<*>", 10, (0, 10, 100)),
        (2, "A", "<*>", 10, (0, 10, 100)),
    ]

    assert build_callable_sequence("GAAAAAAAAA", rows, set(), invariant_rows=invariant_rows) == "GANNNNNNNN"


def test_reference_concordance_ignores_uncertain_positions() -> None:
    summary = reference_concordance("ACGT", "ATNT")
    assert summary.callable_bases == 3
    assert summary.reference_matches == 2
    assert summary.nonreference_bases == 1
    assert summary.identity == 2 / 3


def test_analysis_mask_length_uses_interval_union_and_validates_reference_bounds() -> None:
    assert analysis_mask_length([(0, 5), (3, 8), (10, 12)], reference_length=12) == 10
    with pytest.raises(ValueError, match="outside reference"):
        analysis_mask_length([(0, 13)], reference_length=12)
