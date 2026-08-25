"""Fully called-site PCA sensitivity helpers."""

from __future__ import annotations

import numpy as np

EXPECTED_FULLY_CALLED_MARKERS = {"chloroplast": 1111, "mitochondria": 31}


def select_fully_called_markers(genotypes: np.ndarray) -> np.ndarray:
    """Retain sample-by-marker columns with a genotype in every sample."""
    matrix = np.asarray(genotypes, dtype=float)
    if matrix.ndim != 2:
        raise ValueError("Genotype matrix must be two-dimensional")
    return matrix[:, ~np.isnan(matrix).any(axis=0)]


def classify_pca_sensitivity(correlation: float, p_value: float) -> str:
    if correlation >= 0.90 and p_value < 0.001:
        return "PASS"
    if correlation >= 0.80:
        return "PASS_WITH_CAVEAT"
    return "FAIL"
