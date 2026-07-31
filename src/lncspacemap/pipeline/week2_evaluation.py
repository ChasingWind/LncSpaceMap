"""Week 2C evaluation of real MelD cuTAR projections against held-out truth."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score

LOG = logging.getLogger("lncspacemap.week2_evaluation")

METRIC_DIRECTIONS = {
    "pearson": "higher",
    "spearman": "higher",
    "z_nrmse": "lower",
    "detection_auroc": "higher",
    "detection_auprc": "higher",
    "topk_recall": "higher",
}


def normalize_cutar_truth(
    cutar_counts,
    gene_counts,
    *,
    all_cutar_counts=None,
    scale_factor: float = 1e4,
):
    """Normalize cuTARs by the combined gene-plus-cuTAR spot library."""
    cutar = sp.csr_matrix(cutar_counts, dtype=np.float32)
    genes = sp.csr_matrix(gene_counts, dtype=np.float32)
    if cutar.shape[0] != genes.shape[0]:
        raise ValueError("gene and cuTAR truth spot axes differ")
    all_cutar = (
        cutar
        if all_cutar_counts is None
        else sp.csr_matrix(all_cutar_counts, dtype=np.float32)
    )
    if all_cutar.shape[0] != cutar.shape[0]:
        raise ValueError("target and all-cutar truth spot axes differ")
    totals = np.asarray(all_cutar.sum(axis=1)).ravel()
    totals += np.asarray(genes.sum(axis=1)).ravel()
    factors = np.divide(
        float(scale_factor),
        totals,
        out=np.zeros_like(totals, dtype=np.float32),
        where=totals > 0,
    )
    normalized = sp.diags(factors) @ cutar
    normalized.data = np.log1p(normalized.data)
    return normalized.tocsr()


def _six_metrics(
    truth_continuous: np.ndarray,
    truth_detected: np.ndarray,
    estimate: np.ndarray,
    valid: np.ndarray,
) -> dict[str, float]:
    """Compute metrics while retaining all positives in top-k's denominator."""
    total_positives = int(truth_detected.sum())
    valid_positives = int((truth_detected & valid).sum())
    n_valid = int(valid.sum())
    result = {
        "valid_spots": n_valid,
        "prediction_coverage": float(valid.mean()),
        "truth_detected_spots": total_positives,
        "truth_positive_coverage": (
            float(valid_positives / total_positives) if total_positives else np.nan
        ),
        "valid_truth_prevalence": (
            float(valid_positives / n_valid) if n_valid else np.nan
        ),
        "topk_expected_random": np.nan,
    }
    if n_valid < 3:
        return {**result, **{metric: np.nan for metric in METRIC_DIRECTIONS}}

    observed = truth_continuous[valid].astype(float, copy=False)
    detected = truth_detected[valid]
    predicted = estimate[valid].astype(float, copy=False)
    truth_sd = float(np.std(observed))
    pred_sd = float(np.std(predicted))
    result["pearson"] = (
        float(pearsonr(observed, predicted).statistic)
        if truth_sd > 0 and pred_sd > 0
        else np.nan
    )
    result["spearman"] = (
        float(spearmanr(observed, predicted).statistic)
        if truth_sd > 0 and pred_sd > 0
        else np.nan
    )
    if truth_sd > 0:
        z_truth = (observed - observed.mean()) / truth_sd
        z_pred = (
            (predicted - predicted.mean()) / pred_sd
            if pred_sd > 0
            else np.zeros_like(predicted)
        )
        result["z_nrmse"] = float(np.sqrt(np.mean((z_pred - z_truth) ** 2)))
    else:
        result["z_nrmse"] = np.nan

    if 0 < valid_positives < n_valid:
        result["detection_auroc"] = float(roc_auc_score(detected, predicted))
        result["detection_auprc"] = float(
            average_precision_score(detected, predicted)
        )
    else:
        result["detection_auroc"] = np.nan
        result["detection_auprc"] = np.nan

    if total_positives and n_valid:
        selected = min(total_positives, n_valid)
        top = np.argsort(-predicted, kind="stable")[:selected]
        result["topk_recall"] = float(detected[top].sum() / total_positives)
        result["topk_expected_random"] = float(
            selected * valid_positives / n_valid / total_positives
        )
    else:
        result["topk_recall"] = np.nan
    return result


def evaluate_target(
    truth_continuous: np.ndarray,
    truth_counts: np.ndarray,
    estimate: np.ndarray,
    *,
    n_permutations: int,
    rng: np.random.Generator,
) -> dict[str, float]:
    """Evaluate one target with abstention-aware permutation calibration."""
    truth_continuous = np.asarray(truth_continuous, dtype=float)
    truth_counts = np.asarray(truth_counts, dtype=float)
    estimate = np.asarray(estimate, dtype=float)
    if not (
        truth_continuous.shape == truth_counts.shape == estimate.shape
        and truth_continuous.ndim == 1
    ):
        raise ValueError("truth and prediction target vectors must share one axis")
    if (
        not np.isfinite(truth_continuous).all()
        or not np.isfinite(truth_counts).all()
        or (truth_continuous < 0).any()
        or (truth_counts < 0).any()
    ):
        raise ValueError("truth values are invalid")
    valid = np.isfinite(estimate) & (estimate >= 0)
    detected = truth_counts > 0
    observed = _six_metrics(truth_continuous, detected, estimate, valid)
    prevalence = observed["valid_truth_prevalence"]
    observed["detection_auprc_lift"] = (
        observed["detection_auprc"] / prevalence
        if np.isfinite(observed["detection_auprc"])
        and np.isfinite(prevalence)
        and prevalence > 0
        else np.nan
    )

    null = {metric: [] for metric in METRIC_DIRECTIONS}
    valid_estimate = estimate[valid].copy()
    for _ in range(int(n_permutations)):
        permuted = estimate.copy()
        permuted[valid] = rng.permutation(valid_estimate)
        values = _six_metrics(truth_continuous, detected, permuted, valid)
        for metric in METRIC_DIRECTIONS:
            null[metric].append(values[metric])
    for metric, direction in METRIC_DIRECTIONS.items():
        values = np.asarray(null[metric], dtype=float)
        values = values[np.isfinite(values)]
        null_median = float(np.median(values)) if values.size else np.nan
        observed[f"{metric}_null_median"] = null_median
        value = observed[metric]
        if not np.isfinite(value) or not values.size:
            p_value = np.nan
        elif direction == "higher":
            p_value = float((1 + np.count_nonzero(values >= value)) / (1 + len(values)))
        else:
            p_value = float((1 + np.count_nonzero(values <= value)) / (1 + len(values)))
        observed[f"{metric}_permutation_p"] = p_value
    topk_null = observed["topk_expected_random"]
    observed["topk_recall_lift"] = (
        observed["topk_recall"] / topk_null
        if np.isfinite(observed["topk_recall"])
        and np.isfinite(topk_null)
        and topk_null > 0
        else np.nan
    )
    return observed


def _panel_labels(detection_bin: str) -> list[str]:
    labels = ["all_frozen"]
    if detection_bin == ">10 spots":
        labels.extend(["primary_ge6", "high_confidence_gt10"])
    elif detection_bin == "6-10 spots":
        labels.append("primary_ge6")
    elif detection_bin == "3-5 spots":
        labels.append("stress_3_5")
    return labels


def summarize_metrics(per_target: pd.DataFrame, *, alpha: float) -> pd.DataFrame:
    """Summarize primary spatial and reference-support strata."""
    groups: list[tuple[str, pd.DataFrame]] = []
    for panel in (
        "all_frozen",
        "primary_ge6",
        "high_confidence_gt10",
        "stress_3_5",
    ):
        groups.append(
            (
                panel,
                per_target.loc[
                    per_target["evaluation_panels"].str.split(";").apply(
                        lambda values: panel in values
                    )
                ],
            )
        )
    groups.extend(
        [
            (
                "reference_samples_2",
                per_target.loc[per_target["reference_quantified_samples"].eq(2)],
            ),
            (
                "reference_samples_ge3",
                per_target.loc[per_target["reference_quantified_samples"].ge(3)],
            ),
        ]
    )
    rows = []
    for name, frame in groups:
        row: dict[str, float | int | str] = {
            "stratum": name,
            "targets": len(frame),
            "median_prediction_coverage": float(frame["prediction_coverage"].median()),
            "median_truth_positive_coverage": float(
                frame["truth_positive_coverage"].median()
            ),
        }
        for metric in METRIC_DIRECTIONS:
            row[f"evaluable_{metric}"] = int(frame[metric].notna().sum())
            row[f"median_{metric}"] = float(frame[metric].median())
            row[f"significant_{metric}"] = int(
                frame[f"{metric}_permutation_p"].le(alpha).sum()
            )
        row["median_detection_auprc_lift"] = float(
            frame["detection_auprc_lift"].median()
        )
        row["median_topk_recall_lift"] = float(
            frame["topk_recall_lift"].median()
        )
        rows.append(row)
    return pd.DataFrame(rows)


def scientific_gate(summary: pd.DataFrame, config: dict) -> tuple[str, dict]:
    """Apply a small, predeclared primary-panel decision gate."""
    gate = config["scientific_gate"]
    row = summary.set_index("stratum").loc[str(gate["panel"])]
    checks = {
        "targets": int(row["targets"]) >= int(gate["min_targets"]),
        "spearman": float(row["median_spearman"])
        >= float(gate["min_median_spearman"]),
        "detection_auroc": float(row["median_detection_auroc"])
        >= float(gate["min_median_detection_auroc"]),
        "detection_auprc_lift": float(row["median_detection_auprc_lift"])
        >= float(gate["min_median_detection_auprc_lift"]),
        "topk_recall_lift": float(row["median_topk_recall_lift"])
        >= float(gate["min_median_topk_recall_lift"]),
    }
    directional = sum(
        checks[key]
        for key in (
            "spearman",
            "detection_auroc",
            "detection_auprc_lift",
            "topk_recall_lift",
        )
    )
    passed = checks["targets"] and directional >= int(gate["min_directional_passes"])
    decision = (
        "ADMIT_W2D_CALIBRATION"
        if passed
        else "HOLD_W2D_RUN_BOUNDED_COMPARATOR"
    )
    return decision, {**checks, "directional_passes": directional}


def run_meld_cutar_evaluation(
    prediction_path: Path,
    spatial_gene_path: Path,
    spatial_cutar_truth_path: Path,
    config_path: Path,
    output_dir: Path,
    review_dir: Path,
) -> None:
    """Evaluate the frozen W2B predictions against separately loaded W2A truth."""
    import anndata as ad

    for path in (
        prediction_path,
        spatial_gene_path,
        spatial_cutar_truth_path,
        config_path,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)
    config = yaml.safe_load(config_path.read_text())
    eval_cfg = config["evaluation"]
    prediction = ad.read_h5ad(prediction_path)
    spatial_gene = ad.read_h5ad(spatial_gene_path)
    truth = ad.read_h5ad(spatial_cutar_truth_path)
    if not (
        prediction.obs_names.equals(spatial_gene.obs_names)
        and prediction.obs_names.equals(truth.obs_names)
    ):
        raise ValueError("prediction, gene, and cuTAR truth spot axes differ")
    targets = prediction.var_names.astype(str).tolist()
    truth_idx = truth.var_names.get_indexer(targets)
    if (truth_idx < 0).any():
        raise ValueError("prediction targets are absent from cuTAR truth")
    if "counts" not in spatial_gene.layers or "counts" not in truth.layers:
        raise ValueError("gene and cuTAR truth require layers['counts']")
    predicted = np.asarray(prediction.layers["relative_expression"], dtype=float)
    support = np.asarray(prediction.layers["reference_support"], dtype=float)
    target_counts_sparse = sp.csr_matrix(truth.layers["counts"][:, truth_idx])
    truth_normalized = normalize_cutar_truth(
        target_counts_sparse,
        spatial_gene.layers["counts"],
        all_cutar_counts=truth.layers["counts"],
        scale_factor=float(eval_cfg["scale_factor"]),
    )
    truth_continuous = truth_normalized.toarray()
    truth_counts = target_counts_sparse.toarray()

    rng = np.random.default_rng(int(eval_cfg["random_state"]))
    rows = []
    for index, target in enumerate(targets):
        values = evaluate_target(
            truth_continuous[:, index],
            truth_counts[:, index],
            predicted[:, index],
            n_permutations=int(eval_cfg["permutations"]),
            rng=rng,
        )
        metadata = prediction.var.loc[target]
        detection_bin = str(metadata["spatial_detection_bin"])
        rows.append(
            {
                "target_id": target,
                "panel_rank": int(metadata["panel_rank"]),
                "spatial_detection_bin": detection_bin,
                "reference_selection_tier": str(
                    metadata["reference_selection_tier"]
                ),
                "reference_quantified_samples": int(
                    metadata["reference_quantified_samples"]
                ),
                "median_reference_support": float(np.median(support[:, index])),
                "evaluation_panels": ";".join(_panel_labels(detection_bin)),
                **values,
            }
        )
    per_target = pd.DataFrame(rows).set_index("target_id")
    summary = summarize_metrics(per_target, alpha=float(eval_cfg["alpha"]))
    decision, gate_checks = scientific_gate(summary, eval_cfg)

    output_dir.mkdir(parents=True, exist_ok=True)
    (review_dir / "metrics").mkdir(parents=True, exist_ok=True)
    (review_dir / "manifests").mkdir(parents=True, exist_ok=True)
    per_target.to_csv(
        review_dir / "metrics/w2c_meld_per_target_metrics.tsv", sep="\t"
    )
    summary.to_csv(
        review_dir / "metrics/w2c_meld_stratified_summary.tsv",
        sep="\t",
        index=False,
    )
    bundle = ad.AnnData(
        X=predicted.astype(np.float32),
        obs=prediction.obs.copy(),
        var=prediction.var.copy(),
    )
    bundle.layers["prediction_relative_expression"] = predicted.astype(np.float32)
    bundle.layers["truth_log1p_combined_cp10k"] = truth_continuous.astype(np.float32)
    bundle.layers["truth_counts"] = truth_counts.astype(np.float32)
    bundle.layers["reference_support"] = support.astype(np.float32)
    bundle.obsm["spatial"] = np.asarray(prediction.obsm["spatial"]).copy()
    bundle_path = output_dir / "meld_w2c_evaluation_bundle.h5ad"
    bundle.write_h5ad(bundle_path, compression="gzip")

    manifest = {
        "schema_version": "0.1",
        "stage": "w2c_meld_real_cutar_evaluation",
        "status": "PASS",
        "decision": decision,
        "spots": int(prediction.n_obs),
        "targets": int(prediction.n_vars),
        "metrics": list(METRIC_DIRECTIONS),
        "permutations": int(eval_cfg["permutations"]),
        "prediction_coverage": float(np.isfinite(predicted).mean()),
        "truth_normalization": "log1p(cuTAR/(gene+all_cuTAR)*10000)",
        "detection_truth": "raw_cuTAR_count_gt_0",
        "abstention_policy": "exclude_from_five_metrics_penalize_topk_denominator",
        "spot_axis_check": "PASS",
        "target_axis_check": "PASS",
        "truth_separation_check": "PASS",
        "scientific_gate": gate_checks,
        "large_outputs": {"evaluation_bundle_h5ad": str(bundle_path)},
    }
    (review_dir / "manifests/w2c_meld_evaluation_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    LOG.info(
        "PASS_W2C_MELD_REAL_CUTAR_EVALUATION targets=%d decision=%s",
        prediction.n_vars,
        decision,
    )
