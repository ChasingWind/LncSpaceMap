import pandas as pd

from lncspacemap.evaluation.tuning import (
    compare_tuning_candidate,
    select_candidate,
)
from lncspacemap.evaluation.minimal import METRIC_DIRECTIONS


THRESHOLDS = {
    "spearman": 0.01,
    "pearson": 0.01,
    "z_nrmse": 0.01,
    "detection_auroc": 0.02,
    "detection_auprc": 0.002,
    "topk_recall": 0.02,
}


def _metrics(offset):
    rows = []
    for index in range(3):
        row = {"fold": 0, "gene_id": f"g{index}", "backend": "tangram"}
        for metric, direction in METRIC_DIRECTIONS.items():
            baseline = 0.5 if direction == "higher" else 1.0
            row[metric] = baseline + offset if direction == "higher" else baseline - offset
        rows.append(row)
    return pd.DataFrame(rows)


def test_materially_better_candidate_is_eligible():
    table, summary = compare_tuning_candidate(
        _metrics(0),
        _metrics(0.03),
        candidate_name="better",
        material_delta=THRESHOLDS,
        min_directional_wins=4,
        min_material_wins=2,
    )
    assert summary["eligible"]
    assert summary["directional_wins"] == 6
    assert summary["material_wins"] == 6
    assert table["candidate_advantage"].gt(0).all()
    assert select_candidate([summary]) == "better"


def test_small_or_core_regressing_candidate_is_rejected():
    candidate = _metrics(0.001)
    candidate["spearman"] = _metrics(0)["spearman"] - 0.02
    _, summary = compare_tuning_candidate(
        _metrics(0),
        candidate,
        candidate_name="weak",
        material_delta=THRESHOLDS,
        min_directional_wins=4,
        min_material_wins=2,
    )
    assert not summary["eligible"]
    assert not summary["core_nonregression"]
    assert select_candidate([summary]) is None
