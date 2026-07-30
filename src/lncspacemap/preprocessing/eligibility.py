"""Target eligibility and leakage-free proxy-gene selection."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class EligibilityPolicy:
    min_total_counts: int = 20
    min_detected_cells: int = 10
    min_detected_fraction: float = 0.002
    min_supported_samples: int = 2


def annotate_gene_features(
    feature_qc: pd.DataFrame, gene_table: pd.DataFrame
) -> pd.DataFrame:
    out = feature_qc.copy()
    for column in (
        "gene_id_versioned",
        "gene_name",
        "gene_type",
        "chrom",
        "start",
        "end",
        "strand",
        "annotation_source",
    ):
        mapped = gene_table[column].reindex(out.index)
        if column in out:
            out[column] = out[column].where(out[column].notna(), mapped)
        else:
            out[column] = mapped
    return out


def build_target_catalog(
    feature_qc: pd.DataFrame, policy: EligibilityPolicy
) -> pd.DataFrame:
    required = {
        "feature_type",
        "total_counts",
        "detected_cells",
        "quantified_cells",
        "quantified_sample_count",
        "chrom",
        "start",
        "end",
        "strand",
    }
    missing = sorted(required - set(feature_qc.columns))
    if missing:
        raise ValueError(f"feature QC missing columns: {missing}")
    targets = feature_qc.loc[feature_qc["feature_type"].eq("cuTAR")].copy()
    targets.index.name = "target_id"
    dynamic_min = np.maximum(
        policy.min_detected_cells,
        np.ceil(policy.min_detected_fraction * targets["quantified_cells"]).astype(int),
    )
    enough_counts = targets["total_counts"] >= policy.min_total_counts
    enough_cells = targets["detected_cells"] >= dynamic_min
    enough_samples = (
        targets["quantified_sample_count"] >= policy.min_supported_samples
    )
    coordinates = targets[["chrom", "start", "end", "strand"]].notna().all(axis=1)
    eligible = enough_counts & enough_cells & enough_samples & coordinates
    any_signal = (targets["total_counts"] > 0) & (targets["detected_cells"] > 0)
    targets["eligibility"] = np.where(
        eligible,
        "eligible",
        np.where(any_signal, "exploratory", "insufficient_reference_signal"),
    )
    reasons = []
    for index in targets.index:
        failed = []
        if not enough_counts.loc[index]:
            failed.append("low_total_counts")
        if not enough_cells.loc[index]:
            failed.append("low_cell_detection")
        if not enough_samples.loc[index]:
            failed.append("low_sample_support")
        if not coordinates.loc[index]:
            failed.append("missing_coordinates")
        reasons.append(";".join(failed))
    targets["eligibility_reason"] = reasons
    targets["required_detected_cells"] = dynamic_min
    return targets


def _robust_scale(values: np.ndarray) -> tuple[float, float]:
    center = float(np.median(values))
    mad = float(np.median(np.abs(values - center)))
    return center, max(1.4826 * mad, 1e-6)


def select_masked_gene_proxies(
    annotated_features: pd.DataFrame,
    targets: pd.DataFrame,
    n_genes: int = 250,
    n_folds: int = 5,
) -> pd.DataFrame:
    eligible_targets = targets.loc[targets["eligibility"].eq("eligible")]
    if eligible_targets.empty:
        raise ValueError("no eligible targets available for difficulty matching")
    coding = annotated_features.loc[
        annotated_features["feature_type"].eq("gene")
        & annotated_features["gene_type"].eq("protein_coding")
        & (annotated_features["total_counts"] > 0)
        & (annotated_features["detected_cells"] > 0)
    ].copy()
    if len(coding) < n_folds:
        raise ValueError("insufficient annotated protein-coding proxy genes")

    target_log_counts = np.log1p(eligible_targets["total_counts"].to_numpy(float))
    target_detection = eligible_targets["detected_cell_fraction"].to_numpy(float)
    count_center, count_scale = _robust_scale(target_log_counts)
    detect_center, detect_scale = _robust_scale(target_detection)
    coding["difficulty_distance"] = np.sqrt(
        (
            (np.log1p(coding["total_counts"].to_numpy(float)) - count_center)
            / count_scale
        )
        ** 2
        + (
            (coding["detected_cell_fraction"].to_numpy(float) - detect_center)
            / detect_scale
        )
        ** 2
    )
    selected = coding.sort_values(
        ["difficulty_distance", "detected_cells", "total_counts"],
        kind="mergesort",
    ).head(min(n_genes, len(coding)))
    selected = selected.copy()
    selected.index.name = "gene_id"
    selected["fold"] = np.arange(len(selected), dtype=int) % n_folds
    selected["proxy_role"] = "masked_low_expression_protein_coding"
    return selected
