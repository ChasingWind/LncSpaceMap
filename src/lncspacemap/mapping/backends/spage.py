"""SpaGE adapter with a stable spots-by-targets return contract."""

from __future__ import annotations

import warnings

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


def _bind_spage_axes(
    predicted,
    spot_names: pd.Index,
    targets: list[str],
) -> pd.DataFrame:
    """Bind SpaGE's positional rows to spatial barcodes after strict checks.

    The upstream SpaGE function constructs its return frame without an index,
    even when its spatial input has named rows. Its row loop nevertheless
    preserves the input order, so positional binding is the correct contract.
    """
    if not isinstance(predicted, pd.DataFrame):
        predicted = pd.DataFrame(predicted)
    expected_shape = (len(spot_names), len(targets))
    if predicted.shape != expected_shape:
        raise ValueError(
            f"SpaGE returned shape {predicted.shape}; expected {expected_shape}"
        )
    if predicted.columns.has_duplicates:
        raise ValueError("SpaGE returned duplicate target columns")
    if set(predicted.columns.astype(str)) != set(targets):
        raise ValueError(
            "SpaGE returned target columns that differ from requested targets"
        )
    predicted = predicted.loc[:, targets].copy()
    predicted.index = pd.Index(spot_names, copy=True)
    if not predicted.index.is_unique:
        raise ValueError("spatial barcodes are not unique")
    return predicted.astype(np.float32)


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
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"The behavior of DataFrame\.var with axis=None is deprecated.*",
            category=FutureWarning,
        )
        predicted = SpaGE(
            Spatial_data=spatial_data,
            RNA_data=rna,
            n_pv=min(int(n_pv), len(anchors)),
            genes_to_predict=targets,
        )
    predicted = _bind_spage_axes(predicted, spatial.obs_names, targets)
    predicted[predicted < 0] = 0
    return predicted
