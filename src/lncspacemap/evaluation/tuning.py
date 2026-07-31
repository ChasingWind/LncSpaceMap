"""Paired, leakage-aware comparison for the bounded Week 1 Tangram tuning."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from lncspacemap.evaluation.minimal import METRIC_DIRECTIONS


def read_fold_metrics(
    review_dir: Path,
    *,
    folds: list[int],
    backend: str = "tangram",
) -> pd.DataFrame:
    """Read per-target metrics and enforce a unique fold/target contract."""
    tables = []
    for fold in folds:
        path = review_dir / f"metrics/week1_meld_fold{fold}_{backend}.tsv"
        if not path.is_file():
            raise FileNotFoundError(path)
        table = pd.read_csv(path, sep="\t")
        required = {"gene_id", "backend", "fold", *METRIC_DIRECTIONS}
        missing = required.difference(table.columns)
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        if not table["backend"].eq(backend).all() or not table["fold"].eq(fold).all():
            raise ValueError(f"backend/fold labels do not match {path}")
        tables.append(table)
    result = pd.concat(tables, ignore_index=True)
    if result[["fold", "gene_id"]].duplicated().any():
        raise ValueError("duplicate fold/target rows")
    return result


def compare_tuning_candidate(
    baseline: pd.DataFrame,
    candidate: pd.DataFrame,
    *,
    candidate_name: str,
    material_delta: dict[str, float],
    min_directional_wins: int,
    min_material_wins: int,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Compare identical targets and return metric rows plus an eligibility gate."""
    keys = ["fold", "gene_id"]
    if set(map(tuple, baseline[keys].to_numpy())) != set(
        map(tuple, candidate[keys].to_numpy())
    ):
        raise ValueError(f"{candidate_name} does not use the baseline target panel")
    paired = baseline.merge(
        candidate,
        on=keys,
        validate="one_to_one",
        suffixes=("_baseline", "_candidate"),
    )
    rows = []
    for metric, direction in METRIC_DIRECTIONS.items():
        threshold = float(material_delta[metric])
        baseline_median = float(paired[f"{metric}_baseline"].median())
        candidate_median = float(paired[f"{metric}_candidate"].median())
        advantage = (
            candidate_median - baseline_median
            if direction == "higher"
            else baseline_median - candidate_median
        )
        rows.append(
            {
                "candidate": candidate_name,
                "metric": metric,
                "direction": direction,
                "baseline_median": baseline_median,
                "candidate_median": candidate_median,
                "candidate_advantage": advantage,
                "material_threshold": threshold,
                "directional_win": advantage > 0,
                "material_win": advantage >= threshold,
                "material_regression": advantage <= -threshold,
                "normalized_gain": float(np.clip(advantage / threshold, -2, 2)),
            }
        )
    metric_table = pd.DataFrame(rows)
    directional_wins = int(metric_table["directional_win"].sum())
    material_wins = int(metric_table["material_win"].sum())
    material_regressions = int(metric_table["material_regression"].sum())
    core = metric_table.set_index("metric")
    core_nonregression = not bool(
        core.loc[["spearman", "detection_auroc"], "material_regression"].any()
    )
    eligible = (
        directional_wins >= int(min_directional_wins)
        and material_wins >= int(min_material_wins)
        and core_nonregression
    )
    summary = {
        "candidate": candidate_name,
        "targets": int(len(paired)),
        "folds": sorted(map(int, paired["fold"].unique())),
        "directional_wins": directional_wins,
        "material_wins": material_wins,
        "material_regressions": material_regressions,
        "core_nonregression": core_nonregression,
        "normalized_gain_sum": float(metric_table["normalized_gain"].sum()),
        "eligible": eligible,
    }
    return metric_table, summary


def select_candidate(summaries: list[dict[str, object]]) -> str | None:
    """Select only among candidates that cleared the material-improvement gate."""
    eligible = [row for row in summaries if bool(row["eligible"])]
    if not eligible:
        return None
    eligible.sort(
        key=lambda row: (
            int(row["material_wins"]),
            int(row["directional_wins"]),
            float(row["normalized_gain_sum"]),
            str(row["candidate"]),
        ),
        reverse=True,
    )
    return str(eligible[0]["candidate"])


def write_tuning_gate(path: Path, gate: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(gate, indent=2, sort_keys=True) + "\n")
