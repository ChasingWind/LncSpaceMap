"""Spatially evaluable proxy panels for masked-gene smoke benchmarks."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse


def build_spatial_proxy_folds(
    reference,
    spatial,
    annotated_feature_qc: pd.DataFrame,
    reference_proxy_folds: pd.DataFrame,
    *,
    n_folds: int = 5,
    genes_per_fold: int = 25,
    min_detected_spots: int = 3,
    max_detected_fraction: float = 0.25,
) -> pd.DataFrame:
    """Build deterministic low-expression proxy folds observable in spatial data.

    Spatial expression is used only to establish whether a gene has measurable
    ground truth and to exclude broadly detected genes. It is never supplied to
    a mapping backend for a masked target.
    """
    required = {"feature_type", "gene_type", "total_counts", "detected_cells"}
    missing = required - set(annotated_feature_qc.columns)
    if missing:
        raise ValueError(f"annotated feature QC is missing: {sorted(missing)}")
    if not {"total_counts", "detected_cells"}.issubset(reference_proxy_folds.columns):
        raise ValueError("reference proxy folds lack expression difficulty columns")
    if n_folds < 1 or genes_per_fold < 1:
        raise ValueError("n_folds and genes_per_fold must be positive")

    shared = (
        set(reference.var_names.astype(str))
        & set(spatial.var_names.astype(str))
        & set(annotated_feature_qc.index.astype(str))
    )
    candidates = annotated_feature_qc.loc[
        annotated_feature_qc.index.astype(str).isin(shared)
        & annotated_feature_qc["feature_type"].eq("gene")
        & annotated_feature_qc["gene_type"].eq("protein_coding")
        & annotated_feature_qc["total_counts"].gt(0)
        & annotated_feature_qc["detected_cells"].gt(0)
    ].copy()
    candidates.index = candidates.index.astype(str)
    candidates = candidates[~candidates.index.duplicated(keep=False)]
    if candidates.empty:
        raise ValueError("no shared protein-coding proxy candidates")

    spatial_idx = spatial.var_names.get_indexer(candidates.index)
    if (spatial_idx < 0).any():
        raise AssertionError("shared proxy candidate indexing failed")
    counts = sparse.csr_matrix(spatial.layers["counts"][:, spatial_idx])
    detected_spots = np.asarray((counts > 0).sum(axis=0)).ravel().astype(int)
    spatial_totals = np.asarray(counts.sum(axis=0)).ravel()
    candidates["spatial_detected_spots"] = detected_spots
    candidates["spatial_total_counts"] = spatial_totals
    candidates = candidates[
        candidates["spatial_detected_spots"].ge(min_detected_spots)
        & candidates["spatial_detected_spots"].le(
            int(np.floor(spatial.n_obs * max_detected_fraction))
        )
    ].copy()

    target_count = n_folds * genes_per_fold
    if len(candidates) < target_count:
        raise ValueError(
            "only "
            f"{len(candidates)} spatially evaluable low-expression proxy candidates; "
            f"{target_count} required"
        )

    center_total = float(reference_proxy_folds["total_counts"].median())
    center_detected = float(reference_proxy_folds["detected_cells"].median())
    scale_total = max(
        float(reference_proxy_folds["total_counts"].quantile(0.75)
              - reference_proxy_folds["total_counts"].quantile(0.25)),
        1.0,
    )
    scale_detected = max(
        float(reference_proxy_folds["detected_cells"].quantile(0.75)
              - reference_proxy_folds["detected_cells"].quantile(0.25)),
        1.0,
    )
    candidates["reference_difficulty_distance"] = np.sqrt(
        (
            (candidates["total_counts"].astype(float) - center_total)
            / scale_total
        )
        ** 2
        + (
            (candidates["detected_cells"].astype(float) - center_detected)
            / scale_detected
        )
        ** 2
    )
    candidates["spatial_detection_rank"] = candidates[
        "spatial_detected_spots"
    ].rank(method="average", pct=True)
    candidates["proxy_score"] = (
        candidates["reference_difficulty_distance"]
        + 0.25 * candidates["spatial_detection_rank"]
    )
    selected = candidates.sort_values(
        ["proxy_score", "spatial_detected_spots", "spatial_total_counts"],
        kind="mergesort",
    ).head(target_count)
    selected["fold"] = np.arange(target_count) % n_folds
    selected["proxy_role"] = "meld_spatial_evaluable_low_expression"
    selected["proxy_selection"] = "spatial_fallback"
    return selected
