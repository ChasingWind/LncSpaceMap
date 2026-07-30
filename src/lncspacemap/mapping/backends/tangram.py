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
    tg.pp_adatas(ad_sc, ad_sp, genes=anchors)
    ad_map = tg.map_cells_to_space(
        ad_sc,
        ad_sp,
        mode="cells",
        density_prior="rna_count_based",
        num_epochs=int(num_epochs),
        device=device,
        random_state=int(random_state),
    )
    projected = tg.project_genes(ad_map, ad_sc)
    if not projected.obs_names.equals(spatial.obs_names):
        projected = projected[spatial.obs_names].copy()
    matrix = projected[:, targets].X
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    matrix = np.maximum(np.asarray(matrix, dtype=np.float32), 0)
    return pd.DataFrame(matrix, index=spatial.obs_names, columns=targets)
