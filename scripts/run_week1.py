#!/usr/bin/env python3
"""Run LncSpaceMap Week 1 contract, annotation, and masking preparation."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lncspacemap.pipeline.week1 import run_week1_prepare  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["prepare"])
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--feature-qc", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/week1_spanc_lnc.yaml")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("~/Spatial/data/Spanc-Lnc/processed/week1").expanduser(),
    )
    parser.add_argument("--review-dir", type=Path, default=ROOT / "git_eval")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.review_dir.mkdir(parents=True, exist_ok=True)
    (args.review_dir / "logs").mkdir(exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(args.review_dir / "logs/week1_prepare.log"),
        ],
    )
    for path in (args.reference, args.feature_qc, args.gtf, args.config):
        if not path.is_file():
            raise FileNotFoundError(path)
    run_week1_prepare(
        args.reference,
        args.feature_qc,
        args.gtf,
        args.config,
        args.output_dir,
        args.review_dir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
