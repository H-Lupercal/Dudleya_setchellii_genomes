from organelle_pipeline.haplotypes import summarize_haplotypes


def test_haplotype_summary_uses_jointly_callable_variable_sites() -> None:
    summary = summarize_haplotypes({"s1": "AANA", "s2": "AATA", "s3": "AATT"})
    assert summary.positions == (3,)
    assert summary.sample_haplotypes == {"s1": "H1", "s2": "H1", "s3": "H2"}
    assert summary.counts == {"H1": 2, "H2": 1}


def test_haplotype_summary_retains_variable_site_and_marks_missing_sample_ambiguous() -> None:
    summary = summarize_haplotypes({"s1": "AN", "s2": "AA", "s3": "AT"})
    assert summary.positions == (1,)
    assert summary.sample_haplotypes == {"s1": "AMBIGUOUS", "s2": "H1", "s3": "H2"}
    assert summary.counts == {"H1": 1, "H2": 1}


def test_haplotype_summary_rejects_unequal_sequence_lengths() -> None:
    try:
        summarize_haplotypes({"short": "A", "long": "AA"})
    except ValueError as error:
        assert "equal lengths" in str(error)
    else:
        raise AssertionError("unequal haplotype alignments were silently truncated")
