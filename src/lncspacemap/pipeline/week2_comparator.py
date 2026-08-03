"""Week 2C2 bounded comparison of raw and mask-aware Tangram projections."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from scipy.stats import wilcoxon

from lncspacemap.pipeline.week2_evaluation import (
    METRIC_DIRECTIONS,
    _panel_labels,
    evaluate_target,
    normalize_cutar_truth,
    scientific_gate,
    summarize_metrics,
)

LOG = logging.getLogger("lncspacemap.week2_comparator")

COMPARISON_METRICS = [
    "pearson",
    "spearman",
    "z_nrmse",
    "detection_auroc",
    "detection_auprc",
    "detection_auprc_lift",
    "topk_recall",
    "topk_recall_lift",
]


def paired_method_comparison(
    relative: pd.DataFrame,
    raw: pd.DataFrame,
    *,
    panel: str,
    material_delta: dict[str, float],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create target-level directional deltas and primary-panel statistics."""
    if not relative.index.equals(raw.index):
        raise ValueError("raw and relative metric target axes differ")
    required = {
        "evaluation_panels",
        "spatial_detection_bin",
        "reference_quantified_samples",
        *COMPARISON_METRICS,
    }
    for label, frame in (("relative", relative), ("raw", raw)):
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"{label} metrics are missing columns: {sorted(missing)}")

    paired = relative[
        [
            "spatial_detection_bin",
            "reference_quantified_samples",
            "evaluation_panels",
        ]
    ].copy()
    for metric in COMPARISON_METRICS:
        relative_values = pd.to_numeric(relative[metric], errors="coerce")
        raw_values = pd.to_numeric(raw[metric], errors="coerce")
        improvement = (
            relative_values - raw_values
            if metric == "z_nrmse"
            else raw_values - relative_values
        )
        paired[f"relative_{metric}"] = relative_values
        paired[f"raw_{metric}"] = raw_values
        paired[f"improvement_{metric}"] = improvement

    in_panel = paired["evaluation_panels"].str.split(";").apply(
        lambda values: panel in values
    )
    rows = []
    for metric in COMPARISON_METRICS:
        frame = paired.loc[in_panel]
        valid = frame[
            [
                f"relative_{metric}",
                f"raw_{metric}",
                f"improvement_{metric}",
            ]
        ].dropna()
        improvements = valid[f"improvement_{metric}"].to_numpy(dtype=float)
        if improvements.size and np.any(improvements != 0):
            try:
                p_value = float(
                    wilcoxon(improvements, alternative="greater").pvalue
                )
            except ValueError:
                p_value = np.nan
        else:
            p_value = np.nan
        threshold = float(material_delta.get(metric, 0.0))
        rows.append(
            {
                "panel": panel,
                "metric": metric,
                "targets_paired": len(valid),
                "relative_median": float(valid[f"relative_{metric}"].median()),
                "raw_median": float(valid[f"raw_{metric}"].median()),
                "median_directional_improvement": float(
                    valid[f"improvement_{metric}"].median()
                ),
                "targets_improved": int((improvements > 0).sum()),
                "improved_fraction": (
                    float((improvements > 0).mean()) if improvements.size else np.nan
                ),
                "material_delta": threshold,
                "material_improvement": bool(
                    improvements.size
                    and float(np.median(improvements)) >= threshold
                ),
                "wilcoxon_greater_p": p_value,
            }
        )
    return paired, pd.DataFrame(rows)


def _evaluate_raw(
    prediction,
    spatial_gene,
    truth,
    eval_cfg: dict,
) -> pd.DataFrame:
    targets = prediction.var_names.astype(str).tolist()
    truth_idx = truth.var_names.get_indexer(targets)
    if (truth_idx < 0).any():
        raise ValueError("prediction targets are absent from cuTAR truth")
    target_counts_sparse = sp.csr_matrix(truth.layers["counts"][:, truth_idx])
    truth_continuous = normalize_cutar_truth(
        target_counts_sparse,
        spatial_gene.layers["counts"],
        all_cutar_counts=truth.layers["counts"],
        scale_factor=float(eval_cfg["scale_factor"]),
    ).toarray()
    truth_counts = target_counts_sparse.toarray()
    raw = np.asarray(prediction.layers["raw_projection"], dtype=float)
    support = np.asarray(prediction.layers["reference_support"], dtype=float)
    if not np.isfinite(raw).all() or (raw < 0).any():
        raise ValueError("raw projection contains invalid values")

    rng = np.random.default_rng(int(eval_cfg["random_state"]))
    rows = []
    for index, target in enumerate(targets):
        values = evaluate_target(
            truth_continuous[:, index],
            truth_counts[:, index],
            raw[:, index],
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
    return pd.DataFrame(rows).set_index("target_id")


def run_meld_projection_comparator(
    prediction_path: Path,
    spatial_gene_path: Path,
    spatial_cutar_truth_path: Path,
    relative_metrics_path: Path,
    config_path: Path,
    review_dir: Path,
) -> None:
    """Run the single permitted raw-versus-relative Tangram comparison."""
    import anndata as ad

    for path in (
        prediction_path,
        spatial_gene_path,
        spatial_cutar_truth_path,
        relative_metrics_path,
        config_path,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)
    config = yaml.safe_load(config_path.read_text())
    eval_cfg = config["evaluation"]
    comparator_cfg = config["comparator"]
    prediction = ad.read_h5ad(prediction_path)
    spatial_gene = ad.read_h5ad(spatial_gene_path)
    truth = ad.read_h5ad(spatial_cutar_truth_path)
    if not (
        prediction.obs_names.equals(spatial_gene.obs_names)
        and prediction.obs_names.equals(truth.obs_names)
    ):
        raise ValueError("prediction, gene, and cuTAR truth spot axes differ")
    for layer in ("raw_projection", "relative_expression", "reference_support"):
        if layer not in prediction.layers:
            raise ValueError(f"prediction is missing layer {layer}")

    relative = pd.read_csv(relative_metrics_path, sep="\t", index_col=0)
    relative.index = relative.index.astype(str)
    if not relative.index.equals(prediction.var_names.astype(str)):
        raise ValueError("relative metric and prediction target axes differ")
    raw = _evaluate_raw(prediction, spatial_gene, truth, eval_cfg)
    raw_summary = summarize_metrics(raw, alpha=float(eval_cfg["alpha"]))
    relative_summary = summarize_metrics(relative, alpha=float(eval_cfg["alpha"]))
    raw_summary.insert(0, "method", "raw_projection")
    relative_summary.insert(0, "method", "relative_expression")
    method_summary = pd.concat(
        [relative_summary, raw_summary], ignore_index=True
    )
    paired, paired_summary = paired_method_comparison(
        relative,
        raw,
        panel=str(comparator_cfg["panel"]),
        material_delta=comparator_cfg["material_delta"],
    )

    gate_decision, gate_checks = scientific_gate(raw_summary, eval_cfg)
    if gate_decision == "ADMIT_W2D_CALIBRATION":
        decision = "ADMIT_W2D_CALIBRATION_RAW_PROJECTION"
    else:
        decision = "FREEZE_TANGRAM_NEGATIVE_BASELINE_ENTER_W3_RELATION_MODEL"
    gate_metrics = {
        "spearman",
        "detection_auroc",
        "detection_auprc_lift",
        "topk_recall_lift",
    }
    paired_gate = paired_summary.loc[
        paired_summary["metric"].isin(gate_metrics)
    ]
    comparison_checks = {
        "directional_wins": int(
            paired_gate["median_directional_improvement"].gt(0).sum()
        ),
        "material_wins": int(paired_gate["material_improvement"].sum()),
        "significant_wins": int(paired_gate["wilcoxon_greater_p"].le(0.05).sum()),
    }

    (review_dir / "metrics").mkdir(parents=True, exist_ok=True)
    (review_dir / "manifests").mkdir(parents=True, exist_ok=True)
    raw.to_csv(
        review_dir / "metrics/w2c2_meld_raw_per_target_metrics.tsv", sep="\t"
    )
    paired.to_csv(
        review_dir / "metrics/w2c2_meld_paired_target_deltas.tsv", sep="\t"
    )
    paired_summary.to_csv(
        review_dir / "metrics/w2c2_meld_paired_summary.tsv",
        sep="\t",
        index=False,
    )
    method_summary.to_csv(
        review_dir / "metrics/w2c2_meld_method_summary.tsv",
        sep="\t",
        index=False,
    )
    manifest = {
        "schema_version": "0.1",
        "stage": "w2c2_meld_bounded_projection_comparator",
        "status": "PASS",
        "decision": decision,
        "comparison": "raw_projection_vs_relative_expression",
        "mapping_refit": False,
        "targets": int(prediction.n_vars),
        "spots": int(prediction.n_obs),
        "permutations": int(eval_cfg["permutations"]),
        "primary_panel": str(comparator_cfg["panel"]),
        "raw_scientific_gate": gate_checks,
        "paired_comparison": comparison_checks,
        "terminal_policy": (
            "No further Tangram tuning. Raw must pass the original scientific "
            "gate; otherwise enter Week 3 relation modeling."
        ),
    }
    (review_dir / "manifests/w2c2_meld_comparator_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    LOG.info(
        "PASS_W2C2_BOUNDED_COMPARATOR targets=%d decision=%s",
        prediction.n_vars,
        decision,
    )
