import numpy as np
from scipy import sparse

from lncspacemap.pipeline.week2_evaluation import (
    evaluate_target,
    normalize_cutar_truth,
)


def test_truth_normalization_uses_gene_plus_all_cutar_library():
    cutar = sparse.csr_matrix([[10], [5]], dtype=np.float32)
    all_cutar = sparse.csr_matrix([[10, 10], [5, 5]], dtype=np.float32)
    genes = sparse.csr_matrix([[80, 0], [40, 0]], dtype=np.float32)
    observed = normalize_cutar_truth(
        cutar,
        genes,
        all_cutar_counts=all_cutar,
        scale_factor=100.0,
    ).toarray()
    expected = np.log1p([[10.0], [10.0]])
    np.testing.assert_allclose(observed, expected)


def test_perfect_prediction_has_perfect_six_metrics():
    truth = np.asarray([0.0, 1.0, 2.0, 0.0, 3.0])
    result = evaluate_target(
        truth,
        truth,
        truth,
        n_permutations=20,
        rng=np.random.default_rng(0),
    )
    assert result["pearson"] == 1.0
    assert result["spearman"] == 1.0
    assert result["z_nrmse"] == 0.0
    assert result["detection_auroc"] == 1.0
    assert result["detection_auprc"] == 1.0
    assert result["topk_recall"] == 1.0


def test_topk_denominator_penalizes_abstained_positive():
    truth = np.asarray([1.0, 1.0, 0.0, 0.0])
    estimate = np.asarray([np.nan, 2.0, 1.0, 0.0])
    result = evaluate_target(
        truth,
        truth,
        estimate,
        n_permutations=0,
        rng=np.random.default_rng(0),
    )
    assert result["truth_positive_coverage"] == 0.5
    assert result["topk_recall"] == 0.5
