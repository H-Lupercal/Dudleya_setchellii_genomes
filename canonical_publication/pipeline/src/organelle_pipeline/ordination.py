"""Haploid genotype preparation for organelle-specific PCA."""

from __future__ import annotations

import numpy as np


def prepare_haploid_pca_matrix(genotypes: np.ndarray) -> np.ndarray:
    """Mean-impute missing calls and standardize polymorphic markers."""

    matrix = np.asarray(genotypes, dtype=float).copy()
    if matrix.ndim != 2 or matrix.shape[0] < 2 or matrix.shape[1] < 1:
        raise ValueError("PCA requires a two-dimensional sample-by-marker matrix")
    means = np.nanmean(matrix, axis=0)
    if np.isnan(means).any():
        raise ValueError("PCA contains a marker with no called genotypes")
    missing = np.where(np.isnan(matrix))
    matrix[missing] = means[missing[1]]
    standard_deviations = matrix.std(axis=0)
    polymorphic = standard_deviations > 0
    if not polymorphic.any():
        raise ValueError("PCA requires at least one polymorphic marker")
    matrix = matrix[:, polymorphic]
    means = means[polymorphic]
    standard_deviations = standard_deviations[polymorphic]
    return (matrix - means) / standard_deviations
