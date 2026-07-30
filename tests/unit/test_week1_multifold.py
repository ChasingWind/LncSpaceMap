import json

import pandas as pd

from lncspacemap.evaluation.multifold import aggregate_multifold
from lncspacemap.evaluation.minimal import METRIC_DIRECTIONS


def _write_backend_result(review_dir, backend, fold, offset):
    rows = []
    for index in range(2):
        row = {
            "gene_id": f"gene-{fold}-{index}",
            "backend": backend,
            "fold": fold,
            "truth_detected_spots": 3 + fold + index,
        }
        for metric, direction in METRIC_DIRECTIONS.items():
            base = 0.8 if direction == "higher" else 0.2
            row[metric] = base + offset if direction == "higher" else base - offset
            row[f"{metric}_permutation_p"] = 0.01
        rows.append(row)
    pd.DataFrame(rows).to_csv(
        review_dir / f"metrics/week1_meld_fold{fold}_{backend}.tsv",
        sep="\t",
        index=False,
    )
    manifest = {
        "status": "PASS",
        "leakage_check": "PASS",
        "spot_order_check": "PASS",
    }
    (
        review_dir
        / f"manifests/week1_meld_fold{fold}_{backend}_manifest.json"
    ).write_text(json.dumps(manifest))


def test_multifold_aggregation_selects_consistent_winner(tmp_path):
    (tmp_path / "metrics").mkdir()
    (tmp_path / "manifests").mkdir()
    for fold in (0, 1):
        _write_backend_result(tmp_path, "spage", fold, 0.0)
        _write_backend_result(tmp_path, "tangram", fold, 0.1)

    gate = aggregate_multifold(
        tmp_path,
        fold_ids=[0, 1],
        backends=("spage", "tangram"),
    )

    assert gate["status"] == "PASS"
    assert gate["provisional_winner"] == "tangram"
    assert gate["winner_metric_counts"]["tangram"] == 6
    assert (
        tmp_path / "metrics/week1_meld_multifold_detection_summary.tsv"
    ).is_file()


def test_multifold_aggregation_does_not_break_metric_ties(tmp_path):
    (tmp_path / "metrics").mkdir()
    (tmp_path / "manifests").mkdir()
    _write_backend_result(tmp_path, "spage", 0, 0.0)
    _write_backend_result(tmp_path, "tangram", 0, 0.0)

    gate = aggregate_multifold(
        tmp_path,
        fold_ids=[0],
        backends=("spage", "tangram"),
    )

    assert gate["provisional_winner"] == "inconclusive"
    assert set(gate["per_metric_winner"].values()) == {"tie"}
