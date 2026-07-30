"""Week 1 P0-P2 preparation for SPanC-Lnc."""

from __future__ import annotations

import hashlib
import json
import logging
import platform
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import pandas as pd
import yaml

from lncspacemap.io.gtf import read_gene_gtf
from lncspacemap.preprocessing.eligibility import (
    EligibilityPolicy,
    annotate_gene_features,
    build_target_catalog,
    select_masked_gene_proxies,
)

LOG = logging.getLogger("lncspacemap.week1")


def fast_file_fingerprint(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    size = path.stat().st_size
    digest.update(str(size).encode())
    with path.open("rb") as handle:
        digest.update(handle.read(block_size))
        if size > block_size:
            handle.seek(max(0, size - block_size))
            digest.update(handle.read(block_size))
    return digest.hexdigest()


def _package_versions() -> dict[str, str]:
    result = {}
    for package in ("anndata", "scanpy", "pandas", "numpy", "scipy", "torch"):
        try:
            result[package] = version(package)
        except PackageNotFoundError:
            result[package] = "not-installed"
    return result


def run_week1_prepare(
    reference_path: Path,
    feature_qc_path: Path,
    gtf_path: Path,
    config_path: Path,
    output_dir: Path,
    review_dir: Path,
) -> None:
    import anndata as ad

    config = yaml.safe_load(config_path.read_text())
    output_dir.mkdir(parents=True, exist_ok=True)
    (review_dir / "metrics").mkdir(parents=True, exist_ok=True)
    (review_dir / "manifests").mkdir(parents=True, exist_ok=True)

    reference = ad.read_h5ad(reference_path, backed="r")
    try:
        reference_checks = {
            "has_counts_layer": "counts" in reference.layers,
            "unique_obs_names": reference.obs_names.is_unique,
            "unique_var_names": reference.var_names.is_unique,
            "has_feature_type": "feature_type" in reference.var,
            "has_sample_key": config["reference"]["sample_key"] in reference.obs,
            "has_quantification_mask": "quantified_sample_count" in reference.var,
        }
        reference_row = {
            "status": "PASS" if all(reference_checks.values()) else "FAIL",
            "cells": reference.n_obs,
            "features": reference.n_vars,
            **reference_checks,
        }
    finally:
        reference.file.close()
    pd.DataFrame([reference_row]).to_csv(
        review_dir / "metrics/week1_reference_contract.tsv", sep="\t", index=False
    )
    if reference_row["status"] != "PASS":
        raise ValueError("reference failed Week 1 data contract")

    feature_qc = pd.read_csv(feature_qc_path, sep="\t", index_col=0)
    gene_table = read_gene_gtf(gtf_path)
    annotated = annotate_gene_features(feature_qc, gene_table)
    gene_mask = annotated["feature_type"].eq("gene")
    gene_annotation_fraction = float(annotated.loc[gene_mask, "gene_type"].notna().mean())
    if gene_annotation_fraction < config["annotation"]["min_gene_id_coverage"]:
        raise ValueError(
            f"GENCODE gene ID coverage {gene_annotation_fraction:.3f} below threshold"
        )

    policy = EligibilityPolicy(**config["eligibility"])
    targets = build_target_catalog(annotated, policy)
    proxies = select_masked_gene_proxies(
        annotated,
        targets,
        n_genes=config["masking"]["proxy_genes"],
        n_folds=config["masking"]["folds"],
    )
    annotated.to_csv(output_dir / "annotated_feature_qc.tsv.gz", sep="\t")
    targets.to_csv(output_dir / "target_catalog.tsv.gz", sep="\t")
    proxies.to_csv(output_dir / "masked_gene_folds.tsv", sep="\t")

    target_summary = (
        targets.groupby("eligibility", observed=True)
        .agg(
            targets=("eligibility", "size"),
            median_total_counts=("total_counts", "median"),
            median_detected_cells=("detected_cells", "median"),
            median_supported_samples=("quantified_sample_count", "median"),
        )
        .reset_index()
    )
    target_summary.to_csv(
        review_dir / "metrics/week1_target_eligibility_summary.tsv",
        sep="\t",
        index=False,
    )
    fold_summary = (
        proxies.groupby("fold", observed=True)
        .agg(
            genes=("fold", "size"),
            median_total_counts=("total_counts", "median"),
            median_detected_cells=("detected_cells", "median"),
            median_difficulty_distance=("difficulty_distance", "median"),
        )
        .reset_index()
    )
    fold_summary.to_csv(
        review_dir / "metrics/week1_masked_gene_fold_summary.tsv",
        sep="\t",
        index=False,
    )
    annotation_summary = pd.DataFrame(
        [
            {
                "status": "PASS",
                "gencode_gene_records": len(gene_table),
                "reference_genes": int(gene_mask.sum()),
                "annotated_reference_genes": int(
                    annotated.loc[gene_mask, "gene_type"].notna().sum()
                ),
                "gene_id_coverage": gene_annotation_fraction,
                "protein_coding_proxy_candidates": int(
                    (
                        gene_mask
                        & annotated["gene_type"].eq("protein_coding")
                        & (annotated["detected_cells"] > 0)
                    ).sum()
                ),
            }
        ]
    )
    annotation_summary.to_csv(
        review_dir / "metrics/week1_annotation_summary.tsv", sep="\t", index=False
    )

    manifest = {
        "schema_version": "0.1",
        "stage": "week1_prepare",
        "status": "PASS",
        "seed": config["project"]["seed"],
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": _package_versions(),
        "inputs": {
            "reference": {
                "name": reference_path.name,
                "fingerprint": fast_file_fingerprint(reference_path),
            },
            "feature_qc": {
                "name": feature_qc_path.name,
                "fingerprint": fast_file_fingerprint(feature_qc_path),
            },
            "gtf": {
                "name": gtf_path.name,
                "fingerprint": fast_file_fingerprint(gtf_path),
            },
        },
        "outputs": {
            "target_catalog": "target_catalog.tsv.gz",
            "masked_gene_folds": "masked_gene_folds.tsv",
            "annotated_feature_qc": "annotated_feature_qc.tsv.gz",
        },
    }
    (review_dir / "manifests/week1_prepare_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    LOG.info(
        "PASS_WEEK1_P0_P2_READY_FOR_BASELINES eligible_targets=%d proxy_genes=%d",
        int(targets["eligibility"].eq("eligible").sum()),
        len(proxies),
    )
