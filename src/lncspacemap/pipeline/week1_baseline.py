"""Week 1 MelD masked-gene baseline pipeline."""

from __future__ import annotations

import json
import logging
import platform
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import sparse

from lncspacemap.evaluation import evaluate_predictions
from lncspacemap.io.spatial import load_visium_counts, validate_spatial_counts
from lncspacemap.mapping.backends import run_spage, run_tangram
from lncspacemap.preprocessing.anchors import select_shared_anchors
from lncspacemap.preprocessing.proxies import build_spatial_proxy_folds

LOG = logging.getLogger("lncspacemap.week1_baseline")


def _versions() -> dict[str, str]:
    packages = ("anndata", "scanpy", "tangram-sc", "SpaGE", "numpy", "scipy")
    result = {}
    for package in packages:
        try:
            result[package] = version(package)
        except PackageNotFoundError:
            result[package] = "not-installed-or-unversioned"
    return result


def prepare_spatial(
    matrix_path: Path,
    positions_path: Path,
    output_path: Path,
    review_dir: Path,
) -> None:
    adata = load_visium_counts(matrix_path, positions_path)
    contract = validate_spatial_counts(adata)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")
    (review_dir / "metrics").mkdir(parents=True, exist_ok=True)
    pd.DataFrame([contract]).to_csv(
        review_dir / "metrics/week1_meld_spatial_contract.tsv",
        sep="\t",
        index=False,
    )
    LOG.info("PASS_MELD_SPATIAL_CONTRACT spots=%d genes=%d", adata.n_obs, adata.n_vars)


def _matrix_frame(adata, genes: list[str], *, log_normalize: bool = False) -> pd.DataFrame:
    idx = adata.var_names.get_indexer(genes)
    if (idx < 0).any():
        raise ValueError("truth genes absent from spatial object")
    matrix = adata.layers["counts"][:, idx]
    if log_normalize:
        matrix = sparse.csr_matrix(matrix, dtype=np.float32)
        totals = np.asarray(adata.layers["counts"].sum(axis=1)).ravel()
        scale = np.divide(
            1e4,
            totals,
            out=np.zeros_like(totals, dtype=np.float32),
            where=totals > 0,
        )
        matrix = sparse.diags(scale) @ matrix
        matrix.data = np.log1p(matrix.data)
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return pd.DataFrame(
        np.asarray(matrix, dtype=np.float32),
        index=adata.obs_names,
        columns=genes,
    )


def run_fold(
    reference_path: Path,
    spatial_path: Path,
    annotated_feature_qc_path: Path,
    folds_path: Path,
    config_path: Path,
    output_dir: Path,
    review_dir: Path,
    *,
    backend: str,
    fold: int,
) -> None:
    import anndata as ad

    config = yaml.safe_load(config_path.read_text())
    reference = ad.read_h5ad(reference_path)
    spatial = ad.read_h5ad(spatial_path)
    validate_spatial_counts(spatial)
    feature_qc = pd.read_csv(
        annotated_feature_qc_path, sep="\t", index_col=0, low_memory=False
    )
    folds = pd.read_csv(folds_path, sep="\t", index_col=0, low_memory=False)
    requested = folds.index[folds["fold"].eq(fold)].astype(str).tolist()
    present = [gene for gene in requested if gene in spatial.var_names]
    raw_truth = _matrix_frame(spatial, present)
    min_spots = int(config["targets"]["min_detected_spots"])
    targets = [
        gene for gene in present if int(np.count_nonzero(raw_truth[gene])) >= min_spots
    ]
    target_source = "reference_only_proxy_folds"
    if len(targets) < int(config["targets"]["min_evaluable_targets"]):
        fallback = build_spatial_proxy_folds(
            reference,
            spatial,
            feature_qc,
            folds,
            n_folds=int(config["spatial_proxy_fallback"]["folds"]),
            genes_per_fold=int(
                config["spatial_proxy_fallback"]["genes_per_fold"]
            ),
            min_detected_spots=min_spots,
            max_detected_fraction=float(
                config["spatial_proxy_fallback"]["max_detected_fraction"]
            ),
        )
        (review_dir / "manifests").mkdir(parents=True, exist_ok=True)
        fallback.to_csv(
            review_dir / "manifests/week1_meld_spatial_proxy_folds.tsv",
            sep="\t",
        )
        target_table = fallback.loc[fallback["fold"].eq(fold)].copy()
        requested = target_table.index.astype(str).tolist()
        targets = requested
        raw_truth = _matrix_frame(spatial, targets)
        target_source = "meld_spatial_evaluable_fallback"
        LOG.warning(
            "reference-only fold %d had insufficient spatial truth; "
            "using deterministic MelD proxy panel with %d targets",
            fold,
            len(targets),
        )
    else:
        target_table = folds.loc[targets].copy()
        target_table["proxy_selection"] = "reference_only"
    truth = _matrix_frame(spatial, targets, log_normalize=True)
    anchor_table = select_shared_anchors(
        reference,
        spatial,
        feature_qc,
        targets,
        n_anchors=int(config["anchors"]["n_genes"]),
        min_reference_fraction=float(config["anchors"]["min_reference_fraction"]),
        min_spatial_spots=int(config["anchors"]["min_spatial_spots"]),
    )
    anchors = anchor_table.index.tolist()
    if set(targets) & set(anchors):
        raise AssertionError("FAIL_LEAKAGE_GATE")

    (review_dir / "manifests").mkdir(parents=True, exist_ok=True)
    (review_dir / "metrics").mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    anchor_table.to_csv(
        review_dir / f"manifests/week1_meld_fold{fold}_anchors.tsv", sep="\t"
    )
    target_table["spatial_detected_spots"] = [
        int(np.count_nonzero(raw_truth[g])) for g in targets
    ]
    target_table.to_csv(
        review_dir / f"manifests/week1_meld_fold{fold}_targets.tsv", sep="\t"
    )

    if backend == "spage":
        predicted = run_spage(
            reference,
            spatial,
            anchors,
            targets,
            n_pv=int(config["spage"]["n_pv"]),
        )
    elif backend == "tangram":
        predicted = run_tangram(
            reference,
            spatial,
            anchors,
            targets,
            device=str(config["tangram"]["device"]),
            num_epochs=int(config["tangram"]["num_epochs"]),
            random_state=int(config["project"]["seed"]),
        )
    else:
        raise ValueError(f"unsupported backend: {backend}")

    predicted = predicted.loc[spatial.obs_names, targets]
    per_gene, summary = evaluate_predictions(
        predicted,
        truth,
        n_permutations=int(config["evaluation"]["permutations"]),
        seed=int(config["project"]["seed"]) + int(fold),
    )
    per_gene.insert(0, "backend", backend)
    per_gene.insert(1, "fold", fold)
    metrics_path = review_dir / f"metrics/week1_meld_fold{fold}_{backend}.tsv"
    per_gene.to_csv(metrics_path, sep="\t")
    pd.DataFrame([{**summary, "backend": backend, "fold": fold}]).to_csv(
        review_dir / f"metrics/week1_meld_fold{fold}_{backend}_summary.tsv",
        sep="\t",
        index=False,
    )

    result = ad.AnnData(
        X=predicted.to_numpy(dtype=np.float32),
        obs=spatial.obs.copy(),
        var=pd.DataFrame(index=pd.Index(targets, name="gene_id")),
    )
    result.layers["truth_log1p_1e4"] = truth.to_numpy(dtype=np.float32)
    result.obsm["spatial"] = spatial.obsm["spatial"].copy()
    result.uns["backend"] = backend
    result.uns["fold"] = int(fold)
    result.uns["anchors"] = anchors
    result.write_h5ad(
        output_dir / f"meld_fold{fold}_{backend}_predictions.h5ad", compression="gzip"
    )
    manifest = {
        "schema_version": "0.1",
        "stage": "week1_meld_masked_baseline",
        "status": "PASS",
        "backend": backend,
        "fold": int(fold),
        "spots": int(spatial.n_obs),
        "anchors": len(anchors),
        "targets_requested": len(requested),
        "targets_evaluated": len(targets),
        "target_source": target_source,
        "permutations": int(config["evaluation"]["permutations"]),
        "leakage_check": "PASS",
        "spot_order_check": "PASS",
        "packages": _versions(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "summary": summary,
    }
    (
        review_dir
        / f"manifests/week1_meld_fold{fold}_{backend}_manifest.json"
    ).write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    LOG.info(
        "PASS_WEEK1_BASELINE backend=%s fold=%d spots=%d anchors=%d targets=%d "
        "target_source=%s",
        backend,
        fold,
        spatial.n_obs,
        len(anchors),
        len(targets),
        target_source,
    )
