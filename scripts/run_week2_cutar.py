#!/usr/bin/env python3
"""Run Week 2A cuTAR admission or Week 2B real-target mapping."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lncspacemap.pipeline.week2_cutar import audit_meld_cutar  # noqa: E402
from lncspacemap.pipeline.week2_mapping import run_meld_cutar_mapping  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["audit", "map"])
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--reference-feature-qc", type=Path)
    parser.add_argument("--annotated-feature-qc", type=Path)
    parser.add_argument("--spatial-gene", type=Path, required=True)
    parser.add_argument("--spatial-cutar", type=Path)
    parser.add_argument("--bed", type=Path)
    parser.add_argument("--frozen-targets", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "configs/week2_cutar.yaml",
    )
    parser.add_argument("--review-dir", type=Path, default=ROOT / "git_eval")
    return parser.parse_args()


def _require(args: argparse.Namespace, names: tuple[str, ...]) -> None:
    missing = [name.replace("_", "-") for name in names if getattr(args, name) is None]
    if missing:
        raise SystemExit(
            f"{args.mode}: missing required arguments: "
            + ", ".join(f"--{name}" for name in missing)
        )


def main() -> int:
    args = parse_args()
    (args.review_dir / "logs").mkdir(parents=True, exist_ok=True)
    stage = "week2a_meld_cutar" if args.mode == "audit" else "week2b_meld_mapping"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(args.review_dir / f"logs/{stage}.log"),
        ],
    )
    if args.mode == "audit":
        _require(args, ("reference_feature_qc", "spatial_cutar", "bed"))
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
    else:
        _require(args, ("annotated_feature_qc", "frozen_targets"))
        run_meld_cutar_mapping(
            args.reference,
            args.annotated_feature_qc,
            args.spatial_gene,
            args.frozen_targets,
            args.config,
            args.output_dir,
            args.review_dir,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
