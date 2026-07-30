#!/usr/bin/env python3
"""Prepare MelD spatial input and run leakage-free Week 1 baselines."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lncspacemap.pipeline.week1_baseline import prepare_spatial, run_fold  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="mode", required=True)
    prepare = sub.add_parser("prepare-spatial")
    prepare.add_argument("--matrix", type=Path, required=True)
    prepare.add_argument("--positions", type=Path, required=True)
    prepare.add_argument("--output", type=Path, required=True)
    run = sub.add_parser("run-fold")
    run.add_argument("--reference", type=Path, required=True)
    run.add_argument("--spatial", type=Path, required=True)
    run.add_argument("--annotated-feature-qc", type=Path, required=True)
    run.add_argument("--folds", type=Path, required=True)
    run.add_argument("--backend", choices=["spage", "tangram"], required=True)
    run.add_argument("--fold", type=int, default=0)
    run.add_argument("--output-dir", type=Path, required=True)
    for command in (prepare, run):
        command.add_argument(
            "--config", type=Path, default=ROOT / "configs/week1_meld_baseline.yaml"
        )
        command.add_argument("--review-dir", type=Path, default=ROOT / "git_eval")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.review_dir.mkdir(parents=True, exist_ok=True)
    (args.review_dir / "logs").mkdir(exist_ok=True)
    log_name = f"week1_meld_{args.mode.replace('-', '_')}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(args.review_dir / "logs" / log_name),
        ],
    )
    if args.mode == "prepare-spatial":
        prepare_spatial(args.matrix, args.positions, args.output, args.review_dir)
    else:
        run_fold(
            args.reference,
            args.spatial,
            args.annotated_feature_qc,
            args.folds,
            args.config,
            args.output_dir,
            args.review_dir,
            backend=args.backend,
            fold=args.fold,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
