"""SpaGE adapter with a stable spots-by-targets return contract."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse


def _log_normalize_subset(matrix, library_totals, target_sum: float = 1e4) -> np.ndarray:
    matrix = sparse.csr_matrix(matrix, dtype=np.float32)
    totals = np.asarray(library_totals).ravel()
    scale = np.divide(
        target_sum, totals, out=np.zeros_like(totals, dtype=np.float32), where=totals > 0
    )
    matrix = sparse.diags(scale) @ matrix
    matrix.data = np.log1p(matrix.data)
    return matrix.toarray().astype(np.float32, copy=False)


def run_spage(reference, spatial, anchors: list[str], targets: list[str], *, n_pv=30):
    """Run the official SpaGE function on normalized reference/spatial data."""
    try:
        from SpaGE.main import SpaGE
    except ImportError as exc:
        raise ImportError(
            "SpaGE is unavailable; install/activate the configured SpaGE environment"
        ) from exc

    if set(anchors) & set(targets):
        raise ValueError("masked targets leaked into SpaGE anchors")
    ref_features = anchors + targets
    ref_idx = reference.var_names.get_indexer(ref_features)
    spa_idx = spatial.var_names.get_indexer(anchors)
    if (ref_idx < 0).any() or (spa_idx < 0).any():
        raise ValueError("SpaGE received absent anchors or targets")
    rna = pd.DataFrame(
        _log_normalize_subset(
            reference.layers["counts"][:, ref_idx],
            reference.layers["counts"].sum(axis=1),
        ),
        index=reference.obs_names,
        columns=ref_features,
    )
    spatial_data = pd.DataFrame(
        _log_normalize_subset(
            spatial.layers["counts"][:, spa_idx],
            spatial.layers["counts"].sum(axis=1),
        ),
        index=spatial.obs_names,
        columns=anchors,
    )
    predicted = SpaGE(
        Spatial_data=spatial_data,
        RNA_data=rna,
        n_pv=min(int(n_pv), len(anchors)),
        genes_to_predict=targets,
    )
    predicted = predicted.loc[spatial.obs_names, targets].astype(np.float32)
    predicted[predicted < 0] = 0
    return predicted
