"""Aggregate the five-fold MelD backend evaluation."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from lncspacemap.evaluation.minimal import METRIC_DIRECTIONS


def _read_backend(review_dir: Path, backend: str, fold_ids: list[int]):
    tables = []
    required_columns = {
        "gene_id",
        "backend",
        "fold",
        "truth_detected_spots",
        *METRIC_DIRECTIONS,
        *(f"{metric}_permutation_p" for metric in METRIC_DIRECTIONS),
    }
    for fold in fold_ids:
        path = review_dir / f"metrics/week1_meld_fold{fold}_{backend}.tsv"
        if not path.is_file():
            raise FileNotFoundError(path)
        table = pd.read_csv(path, sep="\t")
        missing = required_columns.difference(table.columns)
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        if not table["backend"].eq(backend).all() or not table["fold"].eq(fold).all():
            raise ValueError(f"backend/fold labels do not match {path}")
        if table["gene_id"].duplicated().any():
            raise ValueError(f"duplicate targets in {path}")
        manifest_path = (
            review_dir
            / f"manifests/week1_meld_fold{fold}_{backend}_manifest.json"
        )
        if not manifest_path.is_file():
            raise FileNotFoundError(manifest_path)
        manifest = json.loads(manifest_path.read_text())
        for key in ("status", "leakage_check", "spot_order_check"):
            if manifest.get(key) != "PASS":
                raise ValueError(f"{manifest_path} failed {key}")
        tables.append(table)
    result = pd.concat(tables, ignore_index=True)
    if result[["fold", "gene_id"]].duplicated().any():
        raise ValueError(f"{backend} has duplicate fold/target pairs")
    if result["gene_id"].duplicated().any():
        raise ValueError(f"{backend} reuses a target across folds")
    return result


def aggregate_multifold(
    review_dir: Path,
    *,
    fold_ids: list[int],
    backends: tuple[str, ...] = ("spage", "tangram"),
    alpha: float = 0.05,
) -> dict[str, object]:
    """Validate, aggregate, and write lightweight five-fold review outputs."""
    if not fold_ids or len(set(fold_ids)) != len(fold_ids):
        raise ValueError("fold_ids must be non-empty and unique")
    if len(backends) < 2 or len(set(backends)) != len(backends):
        raise ValueError("backends must contain at least two unique names")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between zero and one")
    by_backend = {
        backend: _read_backend(review_dir, backend, fold_ids)
        for backend in backends
    }
    reference_backend = backends[0]
    reference = by_backend[reference_backend]
    reference_keys = set(zip(reference["fold"], reference["gene_id"]))
    for backend, table in by_backend.items():
        if set(zip(table["fold"], table["gene_id"])) != reference_keys:
            raise ValueError(f"{backend} does not have the same fold/target panel")

    combined = pd.concat(by_backend.values(), ignore_index=True)
    combined["detection_bin"] = pd.cut(
        combined["truth_detected_spots"],
        bins=[0, 5, 10, np.inf],
        labels=["3-5 spots", "6-10 spots", ">10 spots"],
        include_lowest=True,
    )
    summary_rows = []
    for backend, table in combined.groupby("backend", sort=True):
        row: dict[str, object] = {
            "backend": backend,
            "folds": table["fold"].nunique(),
            "targets": len(table),
            "unique_targets": table["gene_id"].nunique(),
            "median_truth_detected_spots": table["truth_detected_spots"].median(),
        }
        for metric in METRIC_DIRECTIONS:
            row[f"median_{metric}"] = table[metric].median()
            p_column = f"{metric}_permutation_p"
            row[f"significant_{metric}_p{alpha:g}"] = int(
                table[p_column].le(alpha).sum()
            )
        summary_rows.append(row)
    backend_summary = pd.DataFrame(summary_rows)

    fold_summary = (
        combined.groupby(["backend", "fold"], observed=True)
        .agg(
            targets=("gene_id", "size"),
            median_detected_spots=("truth_detected_spots", "median"),
            median_spearman=("spearman", "median"),
            median_pearson=("pearson", "median"),
            median_z_nrmse=("z_nrmse", "median"),
            median_detection_auroc=("detection_auroc", "median"),
            median_detection_auprc=("detection_auprc", "median"),
            median_topk_recall=("topk_recall", "median"),
        )
        .reset_index()
    )
    detection_summary = (
        combined.groupby(["backend", "detection_bin"], observed=True)
        .agg(
            targets=("gene_id", "size"),
            median_detected_spots=("truth_detected_spots", "median"),
            median_spearman=("spearman", "median"),
            median_pearson=("pearson", "median"),
            median_z_nrmse=("z_nrmse", "median"),
            median_detection_auroc=("detection_auroc", "median"),
            median_detection_auprc=("detection_auprc", "median"),
            median_topk_recall=("topk_recall", "median"),
        )
        .reset_index()
    )

    summary_indexed = backend_summary.set_index("backend")
    winners = {}
    for metric, direction in METRIC_DIRECTIONS.items():
        values = summary_indexed[f"median_{metric}"]
        if values.isna().any():
            winner = "not_evaluable"
        else:
            best_value = values.max() if direction == "higher" else values.min()
            tied = values.index[np.isclose(values, best_value, rtol=1e-12, atol=1e-12)]
            winner = str(tied[0]) if len(tied) == 1 else "tie"
        winners[metric] = winner
    eligible_winners = [
        winner for winner in winners.values() if winner in set(backends)
    ]
    winner_counts = pd.Series(eligible_winners, dtype="object").value_counts()
    provisional_winner = (
        str(winner_counts.index[0])
        if not winner_counts.empty and int(winner_counts.iloc[0]) >= 4
        else "inconclusive"
    )

    metrics_dir = review_dir / "metrics"
    manifests_dir = review_dir / "manifests"
    metrics_dir.mkdir(parents=True, exist_ok=True)
    manifests_dir.mkdir(parents=True, exist_ok=True)
    backend_summary.to_csv(
        metrics_dir / "week1_meld_multifold_backend_summary.tsv",
        sep="\t",
        index=False,
    )
    fold_summary.to_csv(
        metrics_dir / "week1_meld_multifold_fold_summary.tsv",
        sep="\t",
        index=False,
    )
    detection_summary.to_csv(
        metrics_dir / "week1_meld_multifold_detection_summary.tsv",
        sep="\t",
        index=False,
    )
    gate = {
        "schema_version": "0.1",
        "stage": "week1_meld_multifold_evaluation",
        "status": "PASS",
        "folds": fold_ids,
        "backends": list(backends),
        "targets_per_backend": {
            backend: int(len(table)) for backend, table in by_backend.items()
        },
        "target_panel_match": "PASS",
        "leakage_checks": "PASS",
        "spot_order_checks": "PASS",
        "per_metric_winner": winners,
        "winner_metric_counts": {
            str(key): int(value) for key, value in winner_counts.items()
        },
        "provisional_winner": provisional_winner,
        "decision": (
            "ADMIT_TARGET_SPECIFIC_DEVELOPMENT"
            if provisional_winner != "inconclusive"
            else "RETAIN_BOTH_BACKENDS"
        ),
    }
    (manifests_dir / "week1_meld_multifold_gate.json").write_text(
        json.dumps(gate, indent=2, sort_keys=True) + "\n"
    )
    return gate
