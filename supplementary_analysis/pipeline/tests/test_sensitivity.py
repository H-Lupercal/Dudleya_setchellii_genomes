import numpy as np
import pytest
from dudleya_supplement.comparative_analysis import finite_pair_spearman, summarize_population_resampling
from dudleya_supplement.sensitivity import (
    classify_fst,
    classify_pi,
    compare_pi,
    procrustes_permutation_test,
    rank_extreme_cases,
    select_eligible,
)
from dudleya_supplement.technical_sensitivity import (
    EXPECTED_FULLY_CALLED_MARKERS,
    classify_pca_sensitivity,
    select_fully_called_markers,
)


def test_scenario_eligibility_uses_its_own_depth_and_breadth() -> None:
    rows = [
        {"sample_id": "a", "breadth_dp3": "0.75", "breadth_dp5": "0.79", "breadth_dp10": "0.91"},
        {"sample_id": "b", "breadth_dp3": "0.69", "breadth_dp5": "0.85", "breadth_dp10": "0.89"},
    ]
    assert select_eligible(rows, eligibility_dp=3, minimum_breadth=0.70) == ["a"]
    assert select_eligible(rows, eligibility_dp=10, minimum_breadth=0.90) == ["a"]


def test_pi_comparison_reports_zero_transitions_without_dividing_by_zero() -> None:
    result = compare_pi({"A": 0.0, "B": 0.1, "C": 0.2}, {"A": 0.01, "B": 0.11, "C": 0.18})
    assert result.zero_to_nonzero == 1
    assert result.nonzero_proportional_changes == (0.1, 0.1)
    assert classify_pi(0.96, 0.10) == "PASS"
    assert classify_pi(0.92, 0.20) == "PASS_WITH_CAVEAT"
    assert classify_pi(0.89, 0.05) == "FAIL"
    assert classify_fst(0.96, 0.05) == "PASS"


def test_extreme_pi_cases_rank_transitions_first_and_limit_to_ten() -> None:
    canonical = {f"P{index}": float(index) for index in range(1, 13)} | {"ZERO": 0.0}
    scenario = {key: value + 1.0 for key, value in canonical.items()}
    rows = rank_extreme_cases(canonical, scenario, scenario="strict", organelle="mitochondria", metric="pi")
    assert len(rows) == 10
    assert rows[0]["population_1"] == "ZERO"
    assert rows[0]["transition_type"] == "zero_to_nonzero"
    assert rows[0]["proportional_change"] == "NA"
    assert [int(row["rank"]) for row in rows] == list(range(1, 11))


def test_extreme_fst_cases_rank_by_absolute_change() -> None:
    canonical = {("A", "B"): 0.1, ("A", "C"): 0.2, ("B", "C"): -0.1}
    scenario = {("A", "B"): 0.11, ("A", "C"): -0.2, ("B", "C"): -0.15}
    rows = rank_extreme_cases(canonical, scenario, scenario="strict", organelle="chloroplast", metric="fst")
    assert [(row["population_1"], row["population_2"]) for row in rows] == [("A", "C"), ("B", "C"), ("A", "B")]
    assert float(rows[0]["absolute_change"]) == pytest.approx(0.4)
    assert all(row["proportional_change"] == "NA" for row in rows)


def test_procrustes_permutation_is_fixed_seed_reproducible() -> None:
    left = np.asarray([[0.0, 0.0, 0.0], [1.0, 0.5, 0.0], [0.0, 1.0, 0.5], [1.0, 1.0, 1.0]])
    right = left @ np.asarray([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    first = procrustes_permutation_test(left, right, permutations=99, seed=424210)
    second = procrustes_permutation_test(left, right, permutations=99, seed=424210)
    assert first == second
    assert first.correlation > 0.999
    assert first.permutations == 99


def test_fst_agreement_excludes_nonfinite_pairs_without_hiding_them() -> None:
    rho, finite_count = finite_pair_spearman([0.1, np.nan, 0.9], [0.2, 0.4, 0.8])
    assert rho == pytest.approx(1.0)
    assert finite_count == 2


def test_fully_called_marker_filter_removes_every_marker_with_missing_genotypes() -> None:
    matrix = np.asarray([[0.0, 1.0, np.nan], [1.0, np.nan, 0.0], [0.0, 1.0, 1.0]])
    filtered = select_fully_called_markers(matrix)
    assert filtered.shape == (3, 1)
    assert filtered[:, 0].tolist() == [0.0, 1.0, 0.0]


def test_complete_site_pca_status_uses_predeclared_thresholds() -> None:
    assert EXPECTED_FULLY_CALLED_MARKERS == {"chloroplast": 1111, "mitochondria": 31}
    assert classify_pca_sensitivity(0.90, 0.0009) == "PASS"
    assert classify_pca_sensitivity(0.90, 0.001) == "PASS_WITH_CAVEAT"
    assert classify_pca_sensitivity(0.80, 0.20) == "PASS_WITH_CAVEAT"
    assert classify_pca_sensitivity(0.799, 0.0001) == "FAIL"


def test_population_resampling_summary_compares_named_outliers_at_common_n() -> None:
    site_rows = [
        {
            "cp_multi_population_haplotypes": str(value),
            "observed_mt_multi_population_haplotypes": "11",
        }
        for value in (8, 10, 11, 12, 14)
    ]
    pi_rows = [
        {"population": population, "sample_size": "4", "nucleotide_diversity": str(value)}
        for population, values in {
            "CY_SIE": (0.90, 1.00),
            "CY_CAS": (0.80, 0.85),
            "OTHER": (0.10, 0.20),
        }.items()
        for value in values
    ]
    summary = summarize_population_resampling(site_rows, pi_rows, named_outliers=("CY_SIE", "CY_CAS"))
    assert summary["marker_count_result"] == "observed_within_cp_distribution"
    assert summary["sample_size_result"] == "named_outlier_medians_remain_top_ranked"
    assert summary["common_sample_size"] == 4
    assert summary["CY_SIE_median_pi"] == pytest.approx(0.95)
    assert summary["CY_CAS_median_pi"] == pytest.approx(0.825)
