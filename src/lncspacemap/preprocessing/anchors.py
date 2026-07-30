"""Leakage-free shared anchor selection."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse


def _moments(matrix) -> tuple[np.ndarray, np.ndarray]:
    matrix = sparse.csr_matrix(matrix, dtype=np.float64)
    mean = np.asarray(matrix.mean(axis=0)).ravel()
    mean_sq = np.asarray(matrix.power(2).mean(axis=0)).ravel()
    variance = np.maximum(mean_sq - mean**2, 0)
    detected = np.asarray((matrix > 0).mean(axis=0)).ravel()
    return variance, detected


def select_shared_anchors(
    reference,
    spatial,
    annotated_feature_qc: pd.DataFrame,
    masked_targets: list[str],
    *,
    n_anchors: int = 2000,
    min_reference_fraction: float = 0.002,
    min_spatial_spots: int = 5,
) -> pd.DataFrame:
    """Rank shared expressed protein-coding genes and exclude every target."""
    masked = set(map(str, masked_targets))
    protein_coding = set(
        annotated_feature_qc.index[
            annotated_feature_qc["feature_type"].eq("gene")
            & annotated_feature_qc["gene_type"].eq("protein_coding")
        ].astype(str)
    )
    shared = sorted(
        (set(reference.var_names.astype(str)) & set(spatial.var_names.astype(str)))
        & protein_coding
        - masked
    )
    if not shared:
        raise ValueError("no shared protein-coding anchors")

    ref_idx = reference.var_names.get_indexer(shared)
    spa_idx = spatial.var_names.get_indexer(shared)
    ref_var, ref_detect = _moments(reference.layers["counts"][:, ref_idx])
    spa_var, spa_detect = _moments(spatial.layers["counts"][:, spa_idx])
    table = pd.DataFrame(
        {
            "gene_id": shared,
            "reference_variance": ref_var,
            "reference_detected_fraction": ref_detect,
            "spatial_variance": spa_var,
            "spatial_detected_spots": np.rint(spa_detect * spatial.n_obs).astype(int),
        }
    ).set_index("gene_id")
    table = table[
        table["reference_detected_fraction"].ge(min_reference_fraction)
        & table["spatial_detected_spots"].ge(min_spatial_spots)
        & table["reference_variance"].gt(0)
        & table["spatial_variance"].gt(0)
    ].copy()
    table["rank_score"] = (
        table["reference_variance"].rank(pct=True)
        + table["spatial_variance"].rank(pct=True)
        + table["reference_detected_fraction"].rank(pct=True)
        + table["spatial_detected_spots"].rank(pct=True)
    )
    table = table.sort_values(
        ["rank_score", "gene_id"], ascending=[False, True]
    ).head(n_anchors)
    if len(table) < min(100, n_anchors):
        raise ValueError(f"only {len(table)} usable shared anchors")
    if masked & set(table.index):
        raise AssertionError("masked targets leaked into anchors")
    table["anchor_rank"] = np.arange(1, len(table) + 1)
    return table
