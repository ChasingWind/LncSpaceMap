#!/usr/bin/env python3
"""Run one Tangram selection round and one untouched-fold confirmation round."""

from __future__ import annotations

import argparse
import copy
import logging
import sys
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lncspacemap.evaluation.tuning import (  # noqa: E402
    compare_tuning_candidate,
    read_fold_metrics,
    select_candidate,
    write_tuning_gate,
)
from lncspacemap.pipeline.week1_baseline import run_fold  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--spatial", type=Path, required=True)
    parser.add_argument("--annotated-feature-qc", type=Path, required=True)
    parser.add_argument("--folds", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "configs/week1_meld_baseline.yaml",
    )
    parser.add_argument("--baseline-review-dir", type=Path, default=ROOT / "git_eval")
    parser.add_argument(
        "--tuning-review-dir",
        type=Path,
        default=ROOT / "git_eval/tuning/week1_tangram",
    )
    return parser.parse_args()


def _candidate_config(
    base: dict,
    candidate: dict,
    destination: Path,
) -> Path:
    config = copy.deepcopy(base)
    config["anchors"]["n_genes"] = int(candidate["anchors"])
    config["tangram"]["num_epochs"] = int(candidate["num_epochs"])
    config["tangram"]["density_prior"] = str(candidate["density_prior"])
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(yaml.safe_dump(config, sort_keys=False))
    return destination


def main() -> int:
    args = parse_args()
    for path in (
        args.reference,
        args.spatial,
        args.annotated_feature_qc,
        args.folds,
        args.config,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)
    config = yaml.safe_load(args.config.read_text())
    tuning = config["tuning"]
    selection_fold = int(tuning["selection_fold"])
    confirmation_folds = list(map(int, tuning["confirmation_folds"]))
    thresholds = {
        key: float(value) for key, value in tuning["material_delta"].items()
    }
    args.tuning_review_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.baseline_review_dir / "logs").mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(
                args.baseline_review_dir / "logs/week1_tangram_tuning.log"
            ),
        ],
    )

    baseline_selection = read_fold_metrics(
        args.baseline_review_dir,
        folds=[selection_fold],
    )
    round1_tables = []
    round1_summaries = []
    candidate_configs = {}
    for name, candidate in tuning["candidates"].items():
        logging.info("START_TUNING_ROUND1 candidate=%s fold=%d", name, selection_fold)
        candidate_dir = args.tuning_review_dir / name
        config_path = _candidate_config(
            config,
            candidate,
            args.output_dir / "configs" / f"{name}.yaml",
        )
        candidate_configs[name] = config_path
        run_fold(
            args.reference,
            args.spatial,
            args.annotated_feature_qc,
            args.folds,
            config_path,
            args.output_dir / name,
            candidate_dir,
            backend="tangram",
            fold=selection_fold,
        )
        candidate_metrics = read_fold_metrics(
            candidate_dir,
            folds=[selection_fold],
        )
        table, summary = compare_tuning_candidate(
            baseline_selection,
            candidate_metrics,
            candidate_name=name,
            material_delta=thresholds,
            min_directional_wins=int(tuning["selection_min_directional_wins"]),
            min_material_wins=int(tuning["selection_min_material_wins"]),
        )
        round1_tables.append(table)
        round1_summaries.append(summary)

    pd.concat(round1_tables, ignore_index=True).to_csv(
        args.tuning_review_dir / "tuning_round1_metrics.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(round1_summaries).to_csv(
        args.tuning_review_dir / "tuning_round1_summary.tsv",
        sep="\t",
        index=False,
    )
    selected = select_candidate(round1_summaries)
    gate: dict[str, object] = {
        "schema_version": "0.1",
        "stage": "week1_tangram_bounded_tuning",
        "status": "PASS",
        "selection_fold": selection_fold,
        "confirmation_folds": confirmation_folds,
        "round1_candidates": round1_summaries,
        "selected_candidate": selected,
    }
    if selected is None:
        gate["decision"] = "STOP_TUNING_ENTER_W2"
        gate["reason"] = "no candidate passed the fold-0 material-improvement gate"
        write_tuning_gate(
            args.tuning_review_dir / "week1_tangram_tuning_gate.json",
            gate,
        )
        logging.info("PASS_TUNING_STOP decision=STOP_TUNING_ENTER_W2")
        return 0

    logging.info(
        "START_TUNING_CONFIRMATION candidate=%s folds=%s",
        selected,
        confirmation_folds,
    )
    selected_dir = args.tuning_review_dir / selected
    for fold in confirmation_folds:
        run_fold(
            args.reference,
            args.spatial,
            args.annotated_feature_qc,
            args.folds,
            candidate_configs[selected],
            args.output_dir / selected,
            selected_dir,
            backend="tangram",
            fold=fold,
        )
    baseline_confirmation = read_fold_metrics(
        args.baseline_review_dir,
        folds=confirmation_folds,
    )
    candidate_confirmation = read_fold_metrics(
        selected_dir,
        folds=confirmation_folds,
    )
    confirmation_table, confirmation_summary = compare_tuning_candidate(
        baseline_confirmation,
        candidate_confirmation,
        candidate_name=selected,
        material_delta=thresholds,
        min_directional_wins=int(tuning["confirmation_min_directional_wins"]),
        min_material_wins=int(tuning["confirmation_min_material_wins"]),
    )
    confirmation_table.to_csv(
        args.tuning_review_dir / "tuning_confirmation_metrics.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame([confirmation_summary]).to_csv(
        args.tuning_review_dir / "tuning_confirmation_summary.tsv",
        sep="\t",
        index=False,
    )
    gate["confirmation"] = confirmation_summary
    if bool(confirmation_summary["eligible"]):
        gate["decision"] = "ACCEPT_TUNED_TANGRAM_ENTER_W2"
        gate["reason"] = "candidate passed untouched-fold material-improvement gate"
    else:
        gate["decision"] = "STOP_TUNING_ENTER_W2"
        gate["reason"] = "candidate did not reproduce material gains on untouched folds"
    write_tuning_gate(
        args.tuning_review_dir / "week1_tangram_tuning_gate.json",
        gate,
    )
    logging.info(
        "PASS_TUNING_COMPLETE selected=%s decision=%s",
        selected,
        gate["decision"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
