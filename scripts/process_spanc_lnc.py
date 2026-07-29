#!/usr/bin/env python3
"""Audit SPanC-Lnc inputs and build matched gene+cuTAR references."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lncspacemap.io.spanc_lnc import (  # noqa: E402
    build_references,
    load_sample_pairs,
    run_audit,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["audit", "build", "all"])
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("~/Spatial/data/Spanc-Lnc").expanduser(),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("~/Spatial/data/Spanc-Lnc/processed/reference").expanduser(),
    )
    parser.add_argument(
        "--review-dir", type=Path, default=ROOT / "git_eval"
    )
    parser.add_argument(
        "--config", type=Path, default=ROOT / "configs/spanc_lnc.yaml"
    )
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
            logging.FileHandler(args.review_dir / "logs/spanc_lnc_preprocess.log"),
        ],
    )
    if not args.data_dir.is_dir():
        raise FileNotFoundError(args.data_dir)
    pairs = load_sample_pairs(args.config)
    required = {
        args.data_dir / name
        for pair in pairs
        for name in (pair.gene_file, pair.cutar_file)
    }
    missing = sorted(str(path) for path in required if not path.is_file())
    if missing:
        raise FileNotFoundError("Missing required inputs:\n" + "\n".join(missing))
    if args.mode in {"audit", "all"}:
        run_audit(args.data_dir, args.review_dir, pairs)
    if args.mode in {"build", "all"}:
        build_references(args.data_dir, args.output_dir, args.review_dir, pairs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
