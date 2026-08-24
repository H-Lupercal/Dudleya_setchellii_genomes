import numpy as np
from dudleya_supplement.sensitivity import (
    classify_fst,
    classify_pi,
    compare_pi,
    procrustes_permutation_test,
    select_eligible,
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


def test_procrustes_permutation_is_fixed_seed_reproducible() -> None:
    left = np.asarray([[0.0, 0.0, 0.0], [1.0, 0.5, 0.0], [0.0, 1.0, 0.5], [1.0, 1.0, 1.0]])
    right = left @ np.asarray([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    first = procrustes_permutation_test(left, right, permutations=99, seed=424210)
    second = procrustes_permutation_test(left, right, permutations=99, seed=424210)
    assert first == second
    assert first.correlation > 0.999
    assert first.permutations == 99
