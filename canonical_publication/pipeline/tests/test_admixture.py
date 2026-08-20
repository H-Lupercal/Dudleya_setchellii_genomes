import pytest
from organelle_pipeline.admixture import pseudo_diploid_alleles, validate_q_matrix


def test_haploid_alleles_are_explicitly_duplicated_for_supplementary_admixture() -> None:
    assert pseudo_diploid_alleles("0", "A", "T") == ("A", "A")
    assert pseudo_diploid_alleles("1", "A", "T") == ("T", "T")
    assert pseudo_diploid_alleles(".", "A", "T") == ("0", "0")


def test_admixture_q_matrix_must_match_samples_k_and_unit_sum() -> None:
    validate_q_matrix(["0.25 0.75", "1 0"], sample_count=2, k=2)
    with pytest.raises(ValueError, match="row count"):
        validate_q_matrix(["0.25 0.75"], sample_count=2, k=2)
    with pytest.raises(ValueError, match="sum"):
        validate_q_matrix(["0.25 0.25"], sample_count=1, k=2)
