"""Tangram adapter with masked targets excluded from mapping genes."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse


def _normalized_subset(counts, indices, library_totals):
    matrix = sparse.csr_matrix(counts[:, indices], dtype=np.float32)
    totals = np.asarray(library_totals).ravel()
    scale = np.divide(
        1e4, totals, out=np.zeros_like(totals, dtype=np.float32), where=totals > 0
    )
    matrix = sparse.diags(scale) @ matrix
    matrix.data = np.log1p(matrix.data)
    return matrix.tocsr()


def _project_target_matrix(
    ad_map,
    target_expression,
    reference_cell_names: pd.Index,
    spatial_spot_names: pd.Index,
    targets: list[str],
) -> pd.DataFrame:
    """Project targets with Tangram's mapping equation under strict axes.

    This is the target-only equivalent of ``tg.project_genes``:
    ``cell_by_spot.T @ cell_by_target``. It avoids the upstream function's
    in-place lower-casing of gene identifiers.
    """
    reference_cell_names = pd.Index(reference_cell_names)
    spatial_spot_names = pd.Index(spatial_spot_names)
    if not ad_map.obs_names.equals(reference_cell_names):
        raise ValueError("Tangram mapping cell axis differs from reference cell order")
    if not ad_map.var_names.equals(spatial_spot_names):
        raise ValueError("Tangram mapping spot axis differs from spatial spot order")

    mapping = ad_map.X
    if sparse.issparse(mapping):
        mapping = mapping.toarray()
    mapping = np.asarray(mapping, dtype=np.float32)
    if sparse.issparse(target_expression):
        target_expression = target_expression.toarray()
    target_expression = np.asarray(target_expression, dtype=np.float32)
    expected_mapping_shape = (len(reference_cell_names), len(spatial_spot_names))
    expected_target_shape = (len(reference_cell_names), len(targets))
    if mapping.shape != expected_mapping_shape:
        raise ValueError(
            f"Tangram mapping shape {mapping.shape}; "
            f"expected {expected_mapping_shape}"
        )
    if target_expression.shape != expected_target_shape:
        raise ValueError(
            f"target expression shape {target_expression.shape}; "
            f"expected {expected_target_shape}"
        )
    if (
        not np.isfinite(mapping).all()
        or not np.isfinite(target_expression).all()
        or (mapping < 0).any()
        or (target_expression < 0).any()
    ):
        raise ValueError("Tangram mapping or target expression is invalid")

    projected = mapping.T @ target_expression
    if not np.isfinite(projected).all() or (projected < 0).any():
        raise ValueError("Tangram target projection is invalid")
    return pd.DataFrame(
        projected.astype(np.float32, copy=False),
        index=spatial_spot_names,
        columns=targets,
    )


def run_tangram(
    reference,
    spatial,
    anchors: list[str],
    targets: list[str],
    *,
    device: str = "cuda:0",
    num_epochs: int = 200,
    random_state: int = 0,
):
    """Map cells to spots on anchors, then project withheld target genes."""
    try:
        import tangram as tg
    except ImportError as exc:
        raise ImportError("Tangram is unavailable in the active environment") from exc

    if set(anchors) & set(targets):
        raise ValueError("masked targets leaked into Tangram anchors")
    ref_features = anchors + targets
    ad_sc = reference[:, ref_features].copy()
    ad_sp = spatial[:, anchors].copy()
    ad_sc.X = _normalized_subset(
        reference.layers["counts"],
        reference.var_names.get_indexer(ref_features),
        reference.layers["counts"].sum(axis=1),
    )
    ad_sp.X = _normalized_subset(
        spatial.layers["counts"],
        spatial.var_names.get_indexer(anchors),
        spatial.layers["counts"].sum(axis=1),
    )
    target_idx = ad_sc.var_names.get_indexer(targets)
    if (target_idx < 0).any():
        raise ValueError("Tangram target genes are absent before preprocessing")
    target_expression = ad_sc.X[:, target_idx].copy()
    tg.pp_adatas(
        ad_sc,
        ad_sp,
        genes=anchors,
        gene_to_lowercase=False,
    )
    if not ad_sc.var_names.is_unique or not ad_sp.var_names.is_unique:
        raise ValueError("Tangram preprocessing produced duplicate gene identifiers")
    ad_map = tg.map_cells_to_space(
        ad_sc,
        ad_sp,
        mode="cells",
        density_prior="rna_count_based",
        num_epochs=int(num_epochs),
        device=device,
        random_state=int(random_state),
    )
    return _project_target_matrix(
        ad_map,
        target_expression,
        ad_sc.obs_names,
        spatial.obs_names,
        targets,
    )
