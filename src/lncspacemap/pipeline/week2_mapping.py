"""Week 2B leakage-free Tangram projection of frozen MelD cuTAR targets."""

from __future__ import annotations

import json
import logging
import platform
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

from lncspacemap.preprocessing.anchors import select_shared_anchors

LOG = logging.getLogger("lncspacemap.week2_mapping")


def _normalized_subset(counts, indices, library_totals):
    matrix = sp.csr_matrix(counts[:, indices], dtype=np.float32)
    totals = np.asarray(library_totals).ravel()
    scale = np.divide(
        1e4,
        totals,
        out=np.zeros_like(totals, dtype=np.float32),
        where=totals > 0,
    )
    matrix = sp.diags(scale) @ matrix
    matrix.data = np.log1p(matrix.data)
    return matrix.tocsr()


def build_quantification_mask(reference, targets: list[str]) -> np.ndarray:
    """Return cells x targets mask from sample-specific quantification metadata."""
    if "sample_batch" not in reference.obs:
        raise ValueError("reference obs is missing sample_batch")
    samples = reference.obs["sample_batch"].astype(str).to_numpy()
    target_var = reference.var.loc[targets]
    mask = np.zeros((reference.n_obs, len(targets)), dtype=np.float32)
    for sample in pd.unique(samples):
        column = f"quantified_{sample}"
        if column not in target_var:
            raise ValueError(f"reference var is missing {column}")
        mask[samples == sample] = (
            target_var[column].fillna(False).astype(bool).to_numpy(dtype=np.float32)
        )
    if (mask.sum(axis=0) == 0).any():
        missing = np.asarray(targets)[mask.sum(axis=0) == 0].tolist()
        raise ValueError(f"targets have no quantified reference cells: {missing[:5]}")
    return mask


def project_mask_aware(
    mapping,
    target_expression,
    quantification_mask,
    *,
    min_support_fraction: float = 0.05,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Project targets and correct sample-specific structural zero filling.

    Returns raw Tangram projection, mapped-cell mean expression, and the
    fraction of mapped reference mass that genuinely quantified each target.
    """
    mapping = mapping.toarray() if sp.issparse(mapping) else np.asarray(mapping)
    expression = (
        target_expression.toarray()
        if sp.issparse(target_expression)
        else np.asarray(target_expression)
    )
    quantified = np.asarray(quantification_mask)
    mapping = mapping.astype(np.float32, copy=False)
    expression = expression.astype(np.float32, copy=False)
    quantified = quantified.astype(np.float32, copy=False)
    if mapping.ndim != 2 or expression.ndim != 2 or quantified.ndim != 2:
        raise ValueError("mapping, expression, and quantification mask must be 2-D")
    if expression.shape != quantified.shape:
        raise ValueError("target expression and quantification mask axes differ")
    if mapping.shape[0] != expression.shape[0]:
        raise ValueError("mapping and target expression cell axes differ")
    if not 0 <= float(min_support_fraction) <= 1:
        raise ValueError("min_support_fraction must be within [0, 1]")
    for name, values in (
        ("mapping", mapping),
        ("target expression", expression),
        ("quantification mask", quantified),
    ):
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f"{name} contains invalid values")

    raw = mapping.T @ expression
    support_mass = mapping.T @ quantified
    total_mass = mapping.sum(axis=0, dtype=np.float64).astype(np.float32)
    support_fraction = np.divide(
        support_mass,
        total_mass[:, None],
        out=np.zeros_like(support_mass, dtype=np.float32),
        where=total_mass[:, None] > 0,
    )
    relative = np.divide(
        raw,
        support_mass,
        out=np.full_like(raw, np.nan, dtype=np.float32),
        where=support_mass > 0,
    )
    relative[support_fraction < float(min_support_fraction)] = np.nan
    return (
        raw.astype(np.float32, copy=False),
        relative.astype(np.float32, copy=False),
        support_fraction.astype(np.float32, copy=False),
    )


def _fit_mapping(reference, spatial, anchors: list[str], config: dict):
    try:
        import anndata as ad
        import tangram as tg
    except ImportError as exc:
        raise ImportError("Tangram and anndata are required for Week 2B") from exc

    ref_idx = reference.var_names.get_indexer(anchors)
    spa_idx = spatial.var_names.get_indexer(anchors)
    if (ref_idx < 0).any() or (spa_idx < 0).any():
        raise ValueError("selected anchors are absent from a mapping input")
    ref_totals = np.asarray(reference.layers["counts"].sum(axis=1)).ravel()
    spa_totals = np.asarray(spatial.layers["counts"].sum(axis=1)).ravel()
    ad_sc = ad.AnnData(
        X=_normalized_subset(reference.layers["counts"], ref_idx, ref_totals),
        obs=reference.obs.copy(),
        var=reference.var.iloc[ref_idx].copy(),
    )
    ad_sp = ad.AnnData(
        X=_normalized_subset(spatial.layers["counts"], spa_idx, spa_totals),
        obs=spatial.obs.copy(),
        var=spatial.var.iloc[spa_idx].copy(),
    )
    tg.pp_adatas(ad_sc, ad_sp, genes=anchors, gene_to_lowercase=False)
    ad_map = tg.map_cells_to_space(
        ad_sc,
        ad_sp,
        mode="cells",
        density_prior=str(config["density_prior"]),
        learning_rate=float(config["learning_rate"]),
        num_epochs=int(config["num_epochs"]),
        device=str(config["device"]),
        random_state=int(config["random_state"]),
    )
    if not ad_map.obs_names.equals(reference.obs_names):
        raise ValueError("Tangram mapping cell order differs from reference")
    if not ad_map.var_names.equals(spatial.obs_names):
        raise ValueError("Tangram mapping spot order differs from spatial input")
    return ad_map, ref_totals


def _versions() -> dict[str, str]:
    result = {}
    for package in ("anndata", "tangram-sc", "numpy", "pandas", "scipy"):
        try:
            result[package] = version(package)
        except PackageNotFoundError:
            result[package] = "not-installed-or-unversioned"
    return result


def run_meld_cutar_mapping(
    reference_path: Path,
    annotated_feature_qc_path: Path,
    spatial_gene_path: Path,
    frozen_targets_path: Path,
    config_path: Path,
    output_dir: Path,
    review_dir: Path,
) -> None:
    """Fit one anchor-only map and project the frozen cuTAR panel in batches."""
    import anndata as ad

    for path in (
        reference_path,
        annotated_feature_qc_path,
        spatial_gene_path,
        frozen_targets_path,
        config_path,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)
    config = yaml.safe_load(config_path.read_text())
    map_cfg = config["mapping"]
    reference = ad.read_h5ad(reference_path)
    spatial = ad.read_h5ad(spatial_gene_path)
    if "counts" not in reference.layers or "counts" not in spatial.layers:
        raise ValueError("reference and spatial inputs require layers['counts']")
    if "spatial" not in spatial.obsm:
        raise ValueError("spatial gene input is missing obsm['spatial']")

    frozen = pd.read_csv(frozen_targets_path, sep="\t", index_col=0)
    frozen.index = frozen.index.astype(str)
    targets = frozen.index[frozen["frozen_primary"].astype(bool)].tolist()
    if not targets:
        raise ValueError("frozen target panel is empty")
    if len(targets) != len(set(targets)):
        raise ValueError("frozen target identifiers are not unique")
    missing = sorted(set(targets) - set(reference.var_names.astype(str)))
    if missing:
        raise ValueError(f"{len(missing)} frozen targets are absent from reference")

    annotated_qc = pd.read_csv(
        annotated_feature_qc_path, sep="\t", index_col=0, low_memory=False
    )
    annotated_qc.index = annotated_qc.index.astype(str)
    anchor_table = select_shared_anchors(
        reference,
        spatial,
        annotated_qc,
        targets,
        n_anchors=int(map_cfg["n_anchors"]),
        min_reference_fraction=float(map_cfg["min_reference_fraction"]),
        min_spatial_spots=int(map_cfg["min_spatial_spots"]),
    )
    anchors = anchor_table.index.astype(str).tolist()
    if set(anchors) & set(targets):
        raise AssertionError("FAIL_W2B_TARGET_LEAKAGE")

    LOG.info(
        "fitting one Tangram map cells=%d spots=%d anchors=%d targets=%d",
        reference.n_obs,
        spatial.n_obs,
        len(anchors),
        len(targets),
    )
    ad_map, ref_totals = _fit_mapping(reference, spatial, anchors, map_cfg)
    mapping = ad_map.X
    batch_size = int(map_cfg["target_batch_size"])
    min_support = float(map_cfg["min_reference_support_fraction"])
    raw_blocks, relative_blocks, support_blocks = [], [], []
    for start in range(0, len(targets), batch_size):
        batch = targets[start : start + batch_size]
        target_idx = reference.var_names.get_indexer(batch)
        expression = _normalized_subset(
            reference.layers["counts"], target_idx, ref_totals
        )
        quantified = build_quantification_mask(reference, batch)
        raw, relative, support = project_mask_aware(
            mapping,
            expression,
            quantified,
            min_support_fraction=min_support,
        )
        raw_blocks.append(raw)
        relative_blocks.append(relative)
        support_blocks.append(support)
        LOG.info(
            "projected targets %d-%d of %d",
            start + 1,
            start + len(batch),
            len(targets),
        )

    raw = np.concatenate(raw_blocks, axis=1)
    relative = np.concatenate(relative_blocks, axis=1)
    support = np.concatenate(support_blocks, axis=1)
    if raw.shape != (spatial.n_obs, len(targets)):
        raise ValueError("projected matrix has unexpected shape")
    if (
        not np.isfinite(raw).all()
        or (raw < 0).any()
        or not np.isfinite(support).all()
        or (support < 0).any()
        or (support > 1 + 1e-5).any()
    ):
        raise ValueError("projected values or reference support are invalid")
    finite_fraction = np.isfinite(relative).mean(axis=0)
    overall_finite_fraction = float(np.isfinite(relative).mean())
    status = (
        "PASS"
        if overall_finite_fraction
        >= float(map_cfg["min_finite_prediction_fraction"])
        else "FAIL"
    )
    target_summary = pd.DataFrame(
        {
            "target_id": targets,
            "finite_spot_fraction": finite_fraction,
            "min_reference_support": support.min(axis=0),
            "median_reference_support": np.median(support, axis=0),
            "max_reference_support": support.max(axis=0),
            "raw_projection_sum": raw.sum(axis=0),
        }
    ).set_index("target_id")

    output_dir.mkdir(parents=True, exist_ok=True)
    (review_dir / "metrics").mkdir(parents=True, exist_ok=True)
    (review_dir / "manifests").mkdir(parents=True, exist_ok=True)
    prediction = ad.AnnData(
        X=relative,
        obs=spatial.obs.copy(),
        var=frozen.loc[targets].copy(),
    )
    prediction.layers["raw_projection"] = raw
    prediction.layers["relative_expression"] = relative.copy()
    prediction.layers["reference_support"] = support
    prediction.obsm["spatial"] = np.asarray(spatial.obsm["spatial"]).copy()
    prediction.uns["lncspacemap_w2b"] = {
        "truth_used_for_training": False,
        "mapping_features": "protein_coding_anchor_genes_only",
        "target_projection": "sample_quantification_mask_aware",
        "minimum_reference_support_fraction": min_support,
    }
    prediction_path = output_dir / "meld_w2b_cutar_predictions.h5ad"
    mapping_path = output_dir / "meld_w2b_cell_to_spot_map.h5ad"
    prediction.write_h5ad(prediction_path, compression="gzip")
    ad_map.write_h5ad(mapping_path, compression="gzip")

    target_summary.to_csv(
        review_dir / "metrics/w2b_meld_target_support.tsv", sep="\t"
    )
    anchor_table.to_csv(
        review_dir / "manifests/w2b_meld_anchors.tsv", sep="\t"
    )
    summary = {
        "status": status,
        "spots": int(spatial.n_obs),
        "reference_cells": int(reference.n_obs),
        "anchors": len(anchors),
        "targets": len(targets),
        "target_batches": int(np.ceil(len(targets) / batch_size)),
        "finite_prediction_fraction": overall_finite_fraction,
        "median_reference_support": float(np.median(support)),
        "minimum_reference_support": float(np.min(support)),
        "truth_used_for_training": False,
        "target_leakage_check": "PASS",
        "cell_axis_check": "PASS",
        "spot_axis_check": "PASS",
    }
    pd.DataFrame([summary]).to_csv(
        review_dir / "metrics/w2b_meld_mapping_contract.tsv",
        sep="\t",
        index=False,
    )
    manifest = {
        "schema_version": "0.1",
        "stage": "w2b_meld_real_cutar_mapping",
        "decision": (
            "ADMIT_W2C_REAL_CUTAR_EVALUATION"
            if status == "PASS"
            else "BLOCK_W2C_INSUFFICIENT_REFERENCE_SUPPORT"
        ),
        **summary,
        "parameters": map_cfg,
        "large_outputs": {
            "prediction_h5ad": str(prediction_path),
            "cell_to_spot_map_h5ad": str(mapping_path),
        },
        "packages": _versions(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
    }
    (review_dir / "manifests/w2b_meld_mapping_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    LOG.info(
        "%s_W2B_MELD_REAL_CUTAR_MAPPING spots=%d anchors=%d targets=%d",
        status,
        spatial.n_obs,
        len(anchors),
        len(targets),
    )
    if status != "PASS":
        raise ValueError(
            "W2B mapping failed the finite prediction support contract: "
            f"{overall_finite_fraction:.4f}"
        )
