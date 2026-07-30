"""Minimal leakage-free baseline metrics."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def evaluate_predictions(predicted: pd.DataFrame, truth: pd.DataFrame):
    """Return per-target Spearman and z-scored NRMSE plus a compact summary."""
    if not predicted.index.equals(truth.index) or not predicted.columns.equals(
        truth.columns
    ):
        raise ValueError("prediction and truth axes are not identical")
    values = predicted.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("predictions contain NaN, infinity, or negative values")
    rows = []
    for gene in truth.columns:
        observed = truth[gene].to_numpy(dtype=float)
        estimate = predicted[gene].to_numpy(dtype=float)
        truth_sd = float(np.std(observed))
        pred_sd = float(np.std(estimate))
        rho = (
            float(spearmanr(observed, estimate).statistic)
            if truth_sd > 0 and pred_sd > 0
            else np.nan
        )
        if truth_sd > 0:
            z_truth = (observed - observed.mean()) / truth_sd
            z_pred = (
                (estimate - estimate.mean()) / pred_sd
                if pred_sd > 0
                else np.zeros_like(estimate)
            )
            z_nrmse = float(np.sqrt(np.mean((z_pred - z_truth) ** 2)))
        else:
            z_nrmse = np.nan
        rows.append(
            {
                "gene_id": gene,
                "spearman": rho,
                "z_nrmse": z_nrmse,
                "truth_detected_spots": int(np.count_nonzero(observed)),
                "truth_total": float(observed.sum()),
            }
        )
    per_gene = pd.DataFrame(rows).set_index("gene_id")
    summary = {
        "targets": int(len(per_gene)),
        "evaluable_spearman": int(per_gene["spearman"].notna().sum()),
        "median_spearman": float(per_gene["spearman"].median()),
        "median_z_nrmse": float(per_gene["z_nrmse"].median()),
    }
    return per_gene, summary
