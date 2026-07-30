import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from lncspacemap.evaluation.minimal import evaluate_predictions
from lncspacemap.io.spatial import validate_spatial_counts
from lncspacemap.preprocessing.anchors import select_shared_anchors
from lncspacemap.preprocessing.proxies import build_spatial_proxy_folds


def _adata(matrix, genes, prefix):
    obj = AnnData(
        X=sparse.csr_matrix(matrix, dtype=np.float32),
        obs=pd.DataFrame(index=[f"{prefix}{i}" for i in range(len(matrix))]),
        var=pd.DataFrame(index=genes),
    )
    obj.layers["counts"] = obj.X.copy()
    return obj


def test_anchor_selection_excludes_masked_target():
    genes = [f"ENSG{i}" for i in range(110)] + ["MASKED"]
    rng = np.random.default_rng(0)
    reference = _adata(rng.poisson(2, (30, len(genes))), genes, "c")
    spatial = _adata(rng.poisson(2, (12, len(genes))), genes, "s")
    annotated = pd.DataFrame(
        {"feature_type": "gene", "gene_type": "protein_coding"}, index=genes
    )
    anchors = select_shared_anchors(
        reference,
        spatial,
        annotated,
        ["MASKED"],
        n_anchors=100,
        min_reference_fraction=0,
        min_spatial_spots=1,
    )
    assert len(anchors) == 100
    assert "MASKED" not in anchors.index


def test_spatial_contract_and_minimal_metrics():
    genes = ["A", "B"]
    spatial = _adata([[1, 0], [2, 1], [3, 2]], genes, "s")
    spatial.obsm["spatial"] = np.array([[0, 0], [1, 1], [2, 2]], dtype=float)
    assert validate_spatial_counts(spatial)["status"] == "PASS"
    truth = pd.DataFrame({"A": [1, 2, 3]}, index=spatial.obs_names)
    predicted = truth.copy()
    per_gene, summary = evaluate_predictions(predicted, truth)
    assert per_gene.loc["A", "spearman"] == 1
    assert per_gene.loc["A", "z_nrmse"] == 0
    assert summary["targets"] == 1


def test_spatial_proxy_fallback_is_evaluable_balanced_and_deterministic():
    genes = [f"ENSG{i:03d}" for i in range(30)]
    rng = np.random.default_rng(1)
    reference = _adata(rng.poisson(1, (40, len(genes))), genes, "c")
    spatial_matrix = np.zeros((20, len(genes)), dtype=np.float32)
    for index in range(len(genes)):
        spatial_matrix[: 3 + index % 2, index] = 1
    spatial = _adata(spatial_matrix, genes, "s")
    annotated = pd.DataFrame(
        {
            "feature_type": "gene",
            "gene_type": "protein_coding",
            "total_counts": np.arange(200, 230),
            "detected_cells": np.arange(180, 210),
        },
        index=genes,
    )
    original = annotated.iloc[:10].copy()
    original["fold"] = np.arange(10) % 5
    first = build_spatial_proxy_folds(
        reference,
        spatial,
        annotated,
        original,
        n_folds=5,
        genes_per_fold=4,
        min_detected_spots=3,
        max_detected_fraction=0.25,
    )
    second = build_spatial_proxy_folds(
        reference,
        spatial,
        annotated,
        original,
        n_folds=5,
        genes_per_fold=4,
        min_detected_spots=3,
        max_detected_fraction=0.25,
    )
    assert first.index.tolist() == second.index.tolist()
    assert first.groupby("fold").size().tolist() == [4, 4, 4, 4, 4]
    assert first["spatial_detected_spots"].ge(3).all()
