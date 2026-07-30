#!/usr/bin/env python3
"""Run and aggregate all MelD masked-gene folds for SpaGE and Tangram."""

from __future__ import annotations

import argparse
import gc
import logging
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lncspacemap.evaluation import aggregate_multifold  # noqa: E402
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
    parser.add_argument("--review-dir", type=Path, default=ROOT / "git_eval")
    parser.add_argument(
        "--backends",
        nargs="+",
        choices=["spage", "tangram"],
        default=["spage", "tangram"],
    )
    parser.add_argument(
        "--fold-ids",
        nargs="+",
        type=int,
        default=[0, 1, 2, 3, 4],
    )
    return parser.parse_args()


def _release_memory() -> None:
    gc.collect()
    try:
        import torch

        if torch.cuda.is_available():
            torch.cuda.empty_cache()
    except ImportError:
        pass


def main() -> int:
    args = parse_args()
    args.review_dir.mkdir(parents=True, exist_ok=True)
    (args.review_dir / "logs").mkdir(exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(
                args.review_dir / "logs/week1_meld_multifold.log"
            ),
        ],
    )
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
    for backend in args.backends:
        for fold in args.fold_ids:
            logging.info("START_MULTIFOLD backend=%s fold=%d", backend, fold)
            run_fold(
                args.reference,
                args.spatial,
                args.annotated_feature_qc,
                args.folds,
                args.config,
                args.output_dir,
                args.review_dir,
                backend=backend,
                fold=fold,
            )
            _release_memory()
    gate = aggregate_multifold(
        args.review_dir,
        fold_ids=args.fold_ids,
        backends=tuple(args.backends),
        alpha=float(config["evaluation"]["alpha"]),
    )
    logging.info(
        "PASS_WEEK1_MULTIFOLD provisional_winner=%s decision=%s",
        gate["provisional_winner"],
        gate["decision"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
