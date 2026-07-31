#!/usr/bin/env python3
"""Audit and align observed MelD cuTAR counts for Week 2."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lncspacemap.pipeline.week2_cutar import audit_meld_cutar  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["audit"])
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--reference-feature-qc", type=Path, required=True)
    parser.add_argument("--spatial-gene", type=Path, required=True)
    parser.add_argument("--spatial-cutar", type=Path, required=True)
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "configs/week2_cutar.yaml",
    )
    parser.add_argument("--review-dir", type=Path, default=ROOT / "git_eval")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    (args.review_dir / "logs").mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(args.review_dir / "logs/week2a_meld_cutar.log"),
        ],
    )
    audit_meld_cutar(
        args.reference,
        args.reference_feature_qc,
        args.spatial_gene,
        args.spatial_cutar,
        args.bed,
        args.config,
        args.output_dir,
        args.review_dir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
