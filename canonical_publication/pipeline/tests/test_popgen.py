import math

from organelle_pipeline.popgen import (
    _hudson_blocks,
    block_bootstrap_hudson_fst,
    callable_nucleotide_diversity,
    haplotype_diversity_from_assignments,
    hudson_fst,
    private_variant_sites,
    private_variant_sites_all,
)


def test_haplotype_diversity_uses_global_assignments_and_excludes_ambiguous_samples() -> None:
    result = haplotype_diversity_from_assignments(("H1", "H1", "H2", "AMBIGUOUS"))

    assert result.haplotype_count == 2
    assert result.assigned_samples == 3
    assert result.ambiguous_samples == 1
    assert math.isclose(result.diversity, 2 / 3)


def test_haplotype_diversity_is_undefined_when_every_sample_is_ambiguous() -> None:
    result = haplotype_diversity_from_assignments(("AMBIGUOUS", "AMBIGUOUS"))

    assert result.haplotype_count == 0
    assert result.assigned_samples == 0
    assert result.ambiguous_samples == 2
    assert math.isnan(result.diversity)


def test_pi_uses_all_pairwise_callable_positions() -> None:
    result = callable_nucleotide_diversity(("AAAA", "AAAT"))
    assert result.differences == 1
    assert result.compared_sites == 4
    assert result.jointly_callable_sites == 4
    assert result.pi == 0.25


def test_pi_excludes_missing_positions_from_denominator() -> None:
    result = callable_nucleotide_diversity(("AANA", "AATA"))
    assert result.differences == 0
    assert result.compared_sites == 3
    assert result.jointly_callable_sites == 3
    assert result.pi == 0.0


def test_pi_excludes_a_position_from_every_pair_if_not_jointly_callable() -> None:
    result = callable_nucleotide_diversity(("AN", "AT", "AT"))
    assert result.jointly_callable_sites == 1
    assert result.compared_sites == 3


def test_hudson_fst_reports_components_and_does_not_clamp_negative_values() -> None:
    fixed = hudson_fst(("AA", "AA"), ("TT", "TT"))
    assert fixed.fst == 1.0
    assert fixed.callable_sites == 2

    same_frequencies = hudson_fst(("A", "T"), ("A", "T"))
    assert same_frequencies.fst < 0
    assert math.isclose(same_frequencies.fst, -1.0)


def test_hudson_fst_callable_count_includes_jointly_callable_invariant_sites() -> None:
    result = hudson_fst(("AA", "AT"), ("AA", "AA"))
    assert result.callable_sites == 2


def test_block_bootstrap_is_deterministic_for_fixed_seed() -> None:
    first = block_bootstrap_hudson_fst(
        ("AAAATTTT", "AAAATTTT"),
        ("TTTTAAAA", "TTTTAAAA"),
        block_size=2,
        replicates=20,
        seed=1729,
    )
    second = block_bootstrap_hudson_fst(
        ("AAAATTTT", "AAAATTTT"),
        ("TTTTAAAA", "TTTTAAAA"),
        block_size=2,
        replicates=20,
        seed=1729,
    )
    assert first == second


def test_hudson_bootstrap_retains_callable_invariant_physical_blocks() -> None:
    blocks = _hudson_blocks(
        ("AATT", "AATT"),
        ("AATT", "AATT"),
        block_size=2,
    )

    assert blocks == ((0.0, 0.0), (0.0, 0.0))


def test_hudson_rejects_unequal_sequence_lengths() -> None:
    try:
        hudson_fst(("AA", "AA"), ("A", "A"))
    except ValueError as error:
        assert "equal lengths" in str(error)
    else:
        raise AssertionError("unequal sequence lengths were silently truncated")


def test_private_variant_count_includes_singletons() -> None:
    groups = {"p1": ("AAAA", "AAAT"), "p2": ("AAAA", "AAAA")}
    assert private_variant_sites("p1", groups, reference="AAAA") == (3,)
    assert private_variant_sites_all(groups, reference="AAAA") == {"p1": (3,), "p2": ()}


def test_private_variants_count_only_nonreference_alleles() -> None:
    groups = {"reference_population": ("A", "A"), "alternate_population": ("T", "T")}
    assert private_variant_sites_all(groups, reference="A") == {
        "reference_population": (),
        "alternate_population": (0,),
    }


def test_strict_private_variant_summary_does_not_treat_missing_as_absence() -> None:
    groups = {"p1": ("T",), "p2": ("N",)}
    assert private_variant_sites_all(groups, reference="A") == {"p1": (0,), "p2": ()}
    assert private_variant_sites_all(groups, reference="A", require_joint_callability=True) == {
        "p1": (),
        "p2": (),
    }
