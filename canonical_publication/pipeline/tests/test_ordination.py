import numpy as np
from organelle_pipeline.ordination import prepare_haploid_pca_matrix


def test_pca_matrix_mean_imputes_and_standardizes_markers() -> None:
    genotypes = np.array([[0.0, 0.0], [1.0, np.nan], [1.0, 1.0]])
    matrix = prepare_haploid_pca_matrix(genotypes)
    assert matrix.shape == (3, 2)
    assert np.allclose(matrix.mean(axis=0), 0.0)
    assert np.allclose(matrix.std(axis=0), 1.0)
