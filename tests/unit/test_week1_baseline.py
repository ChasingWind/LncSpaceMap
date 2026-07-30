import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from lncspacemap.evaluation.minimal import evaluate_predictions
from lncspacemap.io.spatial import validate_spatial_counts
from lncspacemap.preprocessing.anchors import select_shared_anchors


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
