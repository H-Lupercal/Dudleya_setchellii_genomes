"""Sensitivity comparison and deterministic Procrustes testing."""

from __future__ import annotations

import statistics
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
from scipy.stats import spearmanr


@dataclass(frozen=True)
class PiComparison:
    spearman_rho: float
    median_proportional_change: float
    maximum_proportional_change: float
    zero_to_nonzero: int
    nonzero_to_zero: int
    nonzero_proportional_changes: tuple[float, ...]


@dataclass(frozen=True)
class ProcrustesResult:
    correlation: float
    p_value: float
    permutations: int
    seed: int


def select_eligible(rows: Sequence[Mapping[str, str]], eligibility_dp: int, minimum_breadth: float) -> list[str]:
    key = f"breadth_dp{eligibility_dp}"
    return [row["sample_id"] for row in rows if float(row[key]) >= minimum_breadth]


def compare_pi(canonical: Mapping[str, float], scenario: Mapping[str, float]) -> PiComparison:
    common = sorted(set(canonical) & set(scenario))
    left = np.asarray([canonical[key] for key in common], dtype=float)
    right = np.asarray([scenario[key] for key in common], dtype=float)
    rho = float(spearmanr(left, right).statistic) if len(common) >= 2 else float("nan")
    proportional = tuple(round(abs(scenario[key] - canonical[key]) / abs(canonical[key]), 12) for key in common if canonical[key] != 0)
    return PiComparison(
        rho,
        statistics.median(proportional) if proportional else float("nan"),
        max(proportional, default=float("nan")),
        sum(canonical[key] == 0 and scenario[key] != 0 for key in common),
        sum(canonical[key] != 0 and scenario[key] == 0 for key in common),
        proportional,
    )


def classify_pi(rho: float, median_proportional_change: float) -> str:
    if rho >= 0.95 and median_proportional_change <= 0.10:
        return "PASS"
    if rho >= 0.90 and median_proportional_change <= 0.25:
        return "PASS_WITH_CAVEAT"
    return "FAIL"


def classify_fst(rho: float, median_absolute_change: float) -> str:
    if rho >= 0.95 and median_absolute_change <= 0.05:
        return "PASS"
    if rho >= 0.90 and median_absolute_change <= 0.10:
        return "PASS_WITH_CAVEAT"
    return "FAIL"


def _procrustes_correlation(left: np.ndarray, right: np.ndarray) -> float:
    if left.shape != right.shape or left.ndim != 2:
        raise ValueError("Procrustes matrices must have the same two-dimensional shape")
    x = left - left.mean(axis=0)
    y = right - right.mean(axis=0)
    x_norm = np.linalg.norm(x)
    y_norm = np.linalg.norm(y)
    if x_norm == 0 or y_norm == 0:
        raise ValueError("Procrustes matrices must contain variation")
    x /= x_norm
    y /= y_norm
    left_singular, _, right_singular = np.linalg.svd(y.T @ x)
    rotation = left_singular @ right_singular
    aligned = y @ rotation
    return float(np.sum(x * aligned))


def procrustes_permutation_test(
    left: np.ndarray,
    right: np.ndarray,
    *,
    permutations: int = 9999,
    seed: int,
) -> ProcrustesResult:
    observed = _procrustes_correlation(np.asarray(left, dtype=float), np.asarray(right, dtype=float))
    rng = np.random.default_rng(seed)
    exceedances = 0
    for _ in range(permutations):
        permuted = right[rng.permutation(right.shape[0])]
        exceedances += _procrustes_correlation(left, permuted) >= observed - 1e-12
    return ProcrustesResult(observed, (exceedances + 1) / (permutations + 1), permutations, seed)
