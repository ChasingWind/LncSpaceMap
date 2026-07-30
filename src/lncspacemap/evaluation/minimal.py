"""Leakage-free per-target metrics with spatial permutation calibration."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score


METRIC_DIRECTIONS = {
    "spearman": "higher",
    "pearson": "higher",
    "z_nrmse": "lower",
    "detection_auroc": "higher",
    "detection_auprc": "higher",
    "topk_recall": "higher",
}


def _metric_values(observed: np.ndarray, estimate: np.ndarray) -> dict[str, float]:
    truth_sd = float(np.std(observed))
    pred_sd = float(np.std(estimate))
    spearman = (
        float(spearmanr(observed, estimate).statistic)
        if truth_sd > 0 and pred_sd > 0
        else np.nan
    )
    pearson = (
        float(pearsonr(observed, estimate).statistic)
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

    detected = observed > 0
    detected_count = int(detected.sum())
    if 0 < detected_count < len(detected):
        auroc = float(roc_auc_score(detected, estimate))
        auprc = float(average_precision_score(detected, estimate))
        top_indices = np.argsort(-estimate, kind="stable")[:detected_count]
        topk_recall = float(detected[top_indices].mean())
    else:
        auroc = np.nan
        auprc = np.nan
        topk_recall = np.nan
    return {
        "spearman": spearman,
        "pearson": pearson,
        "z_nrmse": z_nrmse,
        "detection_auroc": auroc,
        "detection_auprc": auprc,
        "topk_recall": topk_recall,
    }


def _permutation_calibration(
    observed: np.ndarray,
    estimate: np.ndarray,
    observed_metrics: dict[str, float],
    *,
    n_permutations: int,
    rng: np.random.Generator,
) -> dict[str, float]:
    if n_permutations < 1:
        return {}
    null = {metric: [] for metric in METRIC_DIRECTIONS}
    for _ in range(n_permutations):
        permuted = rng.permutation(estimate)
        values = _metric_values(observed, permuted)
        for metric, value in values.items():
            null[metric].append(value)

    result = {}
    for metric, direction in METRIC_DIRECTIONS.items():
        observed_value = observed_metrics[metric]
        values = np.asarray(null[metric], dtype=float)
        values = values[np.isfinite(values)]
        result[f"{metric}_null_median"] = (
            float(np.median(values)) if values.size else np.nan
        )
        if not np.isfinite(observed_value) or not values.size:
            result[f"{metric}_permutation_p"] = np.nan
        elif direction == "higher":
            result[f"{metric}_permutation_p"] = float(
                (1 + np.count_nonzero(values >= observed_value))
                / (1 + len(values))
            )
        else:
            result[f"{metric}_permutation_p"] = float(
                (1 + np.count_nonzero(values <= observed_value))
                / (1 + len(values))
            )
    return result


def evaluate_predictions(
    predicted: pd.DataFrame,
    truth: pd.DataFrame,
    *,
    n_permutations: int = 0,
    seed: int = 0,
):
    """Return six per-target metrics, permutation nulls, and a summary."""
    if not predicted.index.equals(truth.index) or not predicted.columns.equals(
        truth.columns
    ):
        raise ValueError("prediction and truth axes are not identical")
    values = predicted.to_numpy(dtype=float)
    truth_values = truth.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("predictions contain NaN, infinity, or negative values")
    if not np.isfinite(truth_values).all() or (truth_values < 0).any():
        raise ValueError("truth contains NaN, infinity, or negative values")
    if n_permutations < 0:
        raise ValueError("n_permutations must be non-negative")

    rng = np.random.default_rng(seed)
    rows = []
    for gene in truth.columns:
        observed = truth[gene].to_numpy(dtype=float)
        estimate = predicted[gene].to_numpy(dtype=float)
        metrics = _metric_values(observed, estimate)
        calibration = _permutation_calibration(
            observed,
            estimate,
            metrics,
            n_permutations=n_permutations,
            rng=rng,
        )
        rows.append(
            {
                "gene_id": gene,
                **metrics,
                **calibration,
                "truth_detected_spots": int(np.count_nonzero(observed)),
                "truth_prevalence": float(np.count_nonzero(observed) / len(observed)),
                "truth_total": float(observed.sum()),
            }
        )
    per_gene = pd.DataFrame(rows).set_index("gene_id")
    summary: dict[str, float | int] = {
        "targets": int(len(per_gene)),
        "permutations": int(n_permutations),
    }
    for metric in METRIC_DIRECTIONS:
        summary[f"evaluable_{metric}"] = int(per_gene[metric].notna().sum())
        summary[f"median_{metric}"] = float(per_gene[metric].median())
        p_column = f"{metric}_permutation_p"
        if p_column in per_gene:
            summary[f"significant_{metric}_p05"] = int(
                per_gene[p_column].le(0.05).sum()
            )
    return per_gene, summary
