"""Week 2A contract for observed MelD cuTAR counts and reference targets."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

from lncspacemap.io.spanc_lnc import (
    _load_cutar_bed,
    inspect_text_matrix,
    match_barcodes,
    read_cutar_matrix,
)

LOG = logging.getLogger("lncspacemap.week2_cutar")


def _sha1(path: Path) -> str:
    digest = hashlib.sha1()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _matrix_feature_stats(matrix) -> tuple[np.ndarray, np.ndarray]:
    matrix = matrix.tocsc() if sp.issparse(matrix) else sp.csc_matrix(matrix)
    detected = np.diff(matrix.indptr).astype(np.int64)
    totals = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    return detected, totals


def build_target_catalog(
    spatial_feature_ids: pd.Index,
    spatial_counts,
    reference_var: pd.DataFrame,
    reference_feature_qc: pd.DataFrame,
    bed: pd.DataFrame,
    *,
    min_reference_total_counts: int,
    min_reference_detected_cells: int,
    min_reference_detected_fraction: float,
    min_reference_samples: int,
    min_spatial_detected_spots: int,
) -> pd.DataFrame:
    """Build an exact-ID cuTAR catalog with reference and spatial evidence."""
    spatial_feature_ids = pd.Index(spatial_feature_ids.astype(str), name="target_id")
    for label, index in (
        ("spatial cuTAR", spatial_feature_ids),
        ("reference var", reference_var.index),
        ("reference feature QC", reference_feature_qc.index),
        ("BED", bed.index),
    ):
        if not pd.Index(index).is_unique:
            raise ValueError(f"{label} identifiers are not unique")
    if spatial_counts.shape[1] != len(spatial_feature_ids):
        raise ValueError("spatial count columns do not match spatial feature IDs")

    detected_spots, spatial_totals = _matrix_feature_stats(spatial_counts)
    catalog = pd.DataFrame(index=spatial_feature_ids)
    catalog["spatial_detected_spots"] = detected_spots
    catalog["spatial_total_counts"] = spatial_totals

    reference_types = reference_var.get(
        "feature_type",
        pd.Series("", index=reference_var.index),
    ).astype(str)
    catalog["in_reference"] = catalog.index.isin(reference_var.index)
    catalog["reference_feature_type"] = reference_types.reindex(catalog.index).fillna("")
    catalog["reference_is_cutar"] = catalog["reference_feature_type"].eq("cuTAR")
    catalog["in_reference_qc"] = catalog.index.isin(reference_feature_qc.index)
    catalog["in_bed"] = catalog.index.isin(bed.index)

    qc_columns = {
        "total_counts": "reference_total_counts",
        "detected_cells": "reference_detected_cells",
        "quantified_cells": "reference_quantified_cells",
        "quantified_sample_count": "reference_quantified_samples",
        "selection_tier": "reference_selection_tier",
    }
    for source, destination in qc_columns.items():
        values = (
            reference_feature_qc[source]
            if source in reference_feature_qc
            else pd.Series(np.nan, index=reference_feature_qc.index)
        )
        catalog[destination] = values.reindex(catalog.index).to_numpy()
    catalog["reference_selection_tier"] = (
        catalog["reference_selection_tier"].fillna("missing").astype(str)
    )

    for column in ("chrom", "start", "end", "strand"):
        catalog[column] = bed[column].reindex(catalog.index).to_numpy()
    catalog["chrom"] = catalog["chrom"].fillna("").astype(str)
    catalog["strand"] = catalog["strand"].fillna("").astype(str)

    numeric = (
        "reference_total_counts",
        "reference_detected_cells",
        "reference_quantified_cells",
        "reference_quantified_samples",
    )
    for column in numeric:
        catalog[column] = pd.to_numeric(catalog[column], errors="coerce").fillna(0)
    catalog["reference_detected_fraction"] = np.divide(
        catalog["reference_detected_cells"],
        catalog["reference_quantified_cells"],
        out=np.zeros(len(catalog), dtype=np.float64),
        where=catalog["reference_quantified_cells"].to_numpy() > 0,
    )
    required_detected = np.maximum(
        int(min_reference_detected_cells),
        np.ceil(
            float(min_reference_detected_fraction)
            * catalog["reference_quantified_cells"].to_numpy()
        ).astype(np.int64),
    )
    catalog["required_reference_detected_cells"] = required_detected
    catalog["reference_supported"] = (
        catalog["reference_is_cutar"]
        & catalog["in_reference_qc"]
        & catalog["in_bed"]
        & catalog["reference_total_counts"].ge(min_reference_total_counts)
        & catalog["reference_detected_cells"].ge(required_detected)
        & catalog["reference_quantified_samples"].ge(min_reference_samples)
    )
    catalog["spatial_evaluable"] = catalog["spatial_detected_spots"].ge(
        min_spatial_detected_spots
    )
    catalog["frozen_primary"] = (
        catalog["reference_supported"] & catalog["spatial_evaluable"]
    )
    catalog["spatial_detection_bin"] = pd.cut(
        catalog["spatial_detected_spots"],
        bins=[-1, 0, 2, 5, 10, np.inf],
        labels=["0 spots", "1-2 spots", "3-5 spots", "6-10 spots", ">10 spots"],
    ).astype(str)

    conditions = [
        catalog["frozen_primary"],
        catalog["reference_supported"] & ~catalog["spatial_evaluable"],
        ~catalog["reference_supported"] & catalog["spatial_evaluable"],
    ]
    labels = [
        "primary_evaluable",
        "reference_supported_insufficient_spatial_truth",
        "spatial_observed_insufficient_reference",
    ]
    catalog["target_status"] = np.select(
        conditions,
        labels,
        default="insufficient_reference_and_spatial_evidence",
    )
    return catalog


def _overlap_summary(
    spatial_ids: set[str],
    reference_ids: set[str],
    bed_ids: set[str],
    frozen_count: int,
) -> pd.DataFrame:
    triple = spatial_ids & reference_ids & bed_ids
    rows = [
        ("spatial_cutar_ids", len(spatial_ids)),
        ("reference_cutar_ids", len(reference_ids)),
        ("bed_cutar_ids", len(bed_ids)),
        ("spatial_reference_overlap", len(spatial_ids & reference_ids)),
        ("spatial_bed_overlap", len(spatial_ids & bed_ids)),
        ("reference_bed_overlap", len(reference_ids & bed_ids)),
        ("three_way_overlap", len(triple)),
        ("frozen_primary_targets", int(frozen_count)),
    ]
    return pd.DataFrame(rows, columns=["category", "targets"])


def audit_meld_cutar(
    reference_path: Path,
    reference_feature_qc_path: Path,
    spatial_gene_path: Path,
    spatial_cutar_path: Path,
    bed_path: Path,
    config_path: Path,
    output_dir: Path,
    review_dir: Path,
) -> None:
    """Audit, align, and freeze the observed MelD cuTAR target panel."""
    import anndata as ad

    for path in (
        reference_path,
        reference_feature_qc_path,
        spatial_gene_path,
        spatial_cutar_path,
        bed_path,
        config_path,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)
    config = yaml.safe_load(config_path.read_text())
    contract_cfg = config["contract"]
    eligibility = config["eligibility"]
    source_cfg = config["source"]
    source_sha1 = _sha1(spatial_cutar_path)
    if source_sha1 != str(source_cfg["sha1"]):
        raise ValueError(
            f"{spatial_cutar_path.name}: SHA1 {source_sha1} does not match the "
            f"released {source_cfg['expected_filename']} matrix "
            f"({source_cfg['sha1']}). MelD/MelD.txt is not the cuTAR count "
            "matrix; download data/cuTAR_counts/MelD_cuTAR_mat.txt."
        )

    layout = inspect_text_matrix(spatial_cutar_path)
    chunk_rows = int(config["io"]["chunk_rows"])
    if layout.orientation == "cells_by_features":
        chunk_rows = min(chunk_rows, 8)
    cutar = read_cutar_matrix(
        spatial_cutar_path,
        chunk_rows=chunk_rows,
    )
    reference = ad.read_h5ad(reference_path, backed="r")
    spatial = ad.read_h5ad(spatial_gene_path, backed="r")
    try:
        if "feature_type" not in reference.var:
            raise ValueError("reference var is missing feature_type")
        if "spatial" not in spatial.obsm:
            raise ValueError("spatial gene object is missing obsm['spatial']")
        reference_var = reference.var.copy()
        spatial_obs = spatial.obs.copy()
        spatial_names = spatial.obs_names.astype(str).copy()
        spatial_coords = np.asarray(spatial.obsm["spatial"]).copy()
        reference_cells = int(reference.n_obs)
    finally:
        reference.file.close()
        spatial.file.close()

    strategy, spatial_indices, cutar_indices = match_barcodes(
        spatial_names,
        cutar.obs_names.astype(str),
    )
    barcode_fraction = len(spatial_indices) / len(spatial_names)
    cutar_extra_barcodes = cutar.n_obs - len(cutar_indices)
    aligned = cutar[cutar_indices, :].copy()
    aligned.obs = spatial_obs.iloc[spatial_indices].copy()
    aligned.obs_names = pd.Index(
        spatial_names[spatial_indices],
        name="spot_id",
    )
    aligned.obsm["spatial"] = spatial_coords[spatial_indices]
    aligned.uns["barcode_match_strategy"] = strategy
    aligned.uns["source_matrix"] = str(spatial_cutar_path)
    aligned.uns["counts_source"] = "author_processed_raw_integer_counts"

    feature_qc = pd.read_csv(
        reference_feature_qc_path,
        sep="\t",
        index_col=0,
        low_memory=False,
    )
    feature_qc.index = feature_qc.index.astype(str)
    reference_var.index = reference_var.index.astype(str)
    bed = _load_cutar_bed(bed_path)
    bed.index = bed.index.astype(str)
    catalog = build_target_catalog(
        aligned.var_names,
        aligned.layers["counts"],
        reference_var,
        feature_qc,
        bed,
        min_reference_total_counts=int(eligibility["min_reference_total_counts"]),
        min_reference_detected_cells=int(
            eligibility["min_reference_detected_cells"]
        ),
        min_reference_detected_fraction=float(
            eligibility["min_reference_detected_fraction"]
        ),
        min_reference_samples=int(eligibility["min_reference_samples"]),
        min_spatial_detected_spots=int(
            eligibility["min_spatial_detected_spots"]
        ),
    )
    aligned.var = catalog.copy()

    spatial_ids = set(catalog.index)
    reference_ids = set(
        reference_var.index[
            reference_var["feature_type"].astype(str).eq("cuTAR")
        ]
    )
    bed_ids = set(bed.index)
    frozen = catalog.loc[catalog["frozen_primary"]].copy()
    frozen = frozen.sort_values(
        [
            "reference_selection_tier",
            "spatial_detected_spots",
            "reference_detected_cells",
            "spatial_total_counts",
        ],
        ascending=[True, False, False, False],
    )
    frozen.insert(0, "panel_rank", np.arange(1, len(frozen) + 1))

    spatial_reference_fraction = (
        len(spatial_ids & reference_ids) / len(spatial_ids) if spatial_ids else 0
    )
    spatial_bed_fraction = (
        len(spatial_ids & bed_ids) / len(spatial_ids) if spatial_ids else 0
    )
    coordinate_valid = bool(
        catalog.loc[catalog["in_bed"], "chrom"].str.len().gt(0).all()
        and catalog.loc[catalog["in_bed"], "start"].notna().all()
        and catalog.loc[catalog["in_bed"], "end"].notna().all()
        and (
            catalog.loc[catalog["in_bed"], "end"]
            > catalog.loc[catalog["in_bed"], "start"]
        ).all()
        and catalog.loc[catalog["in_bed"], "strand"].isin(["+", "-"]).all()
    )
    checks = {
        "barcode_overlap": barcode_fraction
        >= float(contract_cfg["min_spatial_barcode_overlap_fraction"]),
        "reference_overlap": spatial_reference_fraction
        >= float(contract_cfg["min_reference_feature_overlap_fraction"]),
        "bed_overlap": spatial_bed_fraction
        >= float(contract_cfg["min_bed_feature_overlap_fraction"]),
        "coordinates": coordinate_valid,
        "frozen_targets": len(frozen) >= int(contract_cfg["min_frozen_targets"]),
    }
    status = "PASS" if all(checks.values()) else "FAIL"

    output_dir.mkdir(parents=True, exist_ok=True)
    (review_dir / "metrics").mkdir(parents=True, exist_ok=True)
    (review_dir / "manifests").mkdir(parents=True, exist_ok=True)
    aligned.write_h5ad(
        output_dir / "meld_cutar_truth_aligned.h5ad",
        compression="gzip",
    )
    catalog.to_csv(output_dir / "w2a_meld_cutar_target_catalog.tsv.gz", sep="\t")
    frozen.to_csv(
        review_dir / "manifests/w2a_meld_frozen_targets.tsv",
        sep="\t",
    )

    overlap = _overlap_summary(
        spatial_ids,
        reference_ids,
        bed_ids,
        len(frozen),
    )
    overlap.to_csv(
        review_dir / "metrics/w2a_meld_target_overlap_summary.tsv",
        sep="\t",
        index=False,
    )
    tier_summary = (
        catalog.groupby(
            ["target_status", "spatial_detection_bin"],
            observed=True,
        )
        .agg(
            targets=("target_status", "size"),
            median_spatial_total_counts=("spatial_total_counts", "median"),
            median_reference_detected_cells=("reference_detected_cells", "median"),
            median_reference_samples=("reference_quantified_samples", "median"),
        )
        .reset_index()
    )
    tier_summary.to_csv(
        review_dir / "metrics/w2a_meld_target_tier_summary.tsv",
        sep="\t",
        index=False,
    )
    contract = pd.DataFrame(
        [
            {
                "status": status,
                "matrix_orientation": layout.orientation,
                "matrix_header_style": layout.header_style,
                "matrix_rows": layout.n_rows,
                "matrix_columns": layout.n_columns,
                "spatial_gene_spots": len(spatial_names),
                "spatial_cutar_barcodes": cutar.n_obs,
                "matched_spots": len(spatial_indices),
                "barcode_overlap_fraction": barcode_fraction,
                "extra_cutar_barcodes": cutar_extra_barcodes,
                "barcode_strategy": strategy,
                "spatial_cutars": len(spatial_ids),
                "spatial_reference_overlap_fraction": spatial_reference_fraction,
                "spatial_bed_overlap_fraction": spatial_bed_fraction,
                "coordinate_contract": "PASS" if coordinate_valid else "FAIL",
                "reference_cells": reference_cells,
                "frozen_primary_targets": len(frozen),
            }
        ]
    )
    contract.to_csv(
        review_dir / "metrics/w2a_meld_cutar_contract.tsv",
        sep="\t",
        index=False,
    )
    manifest = {
        "schema_version": "0.1",
        "stage": "w2a_meld_cutar_contract",
        "status": status,
        "decision": (
            "ADMIT_W2B_REAL_CUTAR_MAPPING"
            if status == "PASS"
            else "BLOCK_W2B_REPAIR_INPUT_CONTRACT"
        ),
        "checks": {key: "PASS" if value else "FAIL" for key, value in checks.items()},
        "matrix": {
            "path": str(spatial_cutar_path),
            "source_filename": str(source_cfg["expected_filename"]),
            "source_sha1": source_sha1,
            "orientation": layout.orientation,
            "header_style": layout.header_style,
            "chunk_rows": chunk_rows,
            "spots": int(aligned.n_obs),
            "cutars": int(aligned.n_vars),
        },
        "barcode_match": {
            "strategy": strategy,
            "matched_spots": len(spatial_indices),
            "spatial_gene_spots": len(spatial_names),
            "overlap_fraction": barcode_fraction,
            "extra_cutar_barcodes": cutar_extra_barcodes,
        },
        "feature_match": {
            "spatial_reference_overlap_fraction": spatial_reference_fraction,
            "spatial_bed_overlap_fraction": spatial_bed_fraction,
            "exact_id_matching": True,
            "genome_build_declared": str(config["annotation"]["genome_build"]),
        },
        "targets": {
            "frozen_primary": len(frozen),
            "min_spatial_detected_spots": int(
                eligibility["min_spatial_detected_spots"]
            ),
            "reference_thresholds": {
                "total_counts": int(eligibility["min_reference_total_counts"]),
                "detected_cells_floor": int(
                    eligibility["min_reference_detected_cells"]
                ),
                "detected_fraction": float(
                    eligibility["min_reference_detected_fraction"]
                ),
                "samples": int(eligibility["min_reference_samples"]),
            },
        },
        "large_outputs": {
            "aligned_truth_h5ad": str(
                output_dir / "meld_cutar_truth_aligned.h5ad"
            ),
            "full_target_catalog": str(
                output_dir / "w2a_meld_cutar_target_catalog.tsv.gz"
            ),
        },
    }
    (
        review_dir / "manifests/w2a_meld_cutar_manifest.json"
    ).write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    LOG.info(
        "%s_W2A_MELD_CUTAR_CONTRACT spots=%d cutars=%d frozen=%d",
        status,
        aligned.n_obs,
        aligned.n_vars,
        len(frozen),
    )
    if status != "PASS":
        failed = [key for key, value in checks.items() if not value]
        raise ValueError(f"W2A input contract failed: {failed}")
