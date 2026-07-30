"""Load Visium-style 10x matrices under a small, explicit data contract."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


POSITION_COLUMNS = [
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres",
]


def read_tissue_positions(path: Path) -> pd.DataFrame:
    """Read either legacy headerless or Space Ranger headered positions."""
    table = pd.read_csv(path)
    if "barcode" not in table.columns:
        table = pd.read_csv(path, header=None, names=POSITION_COLUMNS)
    missing = set(POSITION_COLUMNS) - set(table.columns)
    if missing:
        raise ValueError(f"position table is missing columns: {sorted(missing)}")
    if table["barcode"].duplicated().any():
        raise ValueError("position table contains duplicate barcodes")
    return table.set_index("barcode")


def load_visium_counts(
    matrix_path: Path,
    positions_path: Path,
    *,
    in_tissue_only: bool = True,
):
    """Load a 10x H5 matrix with stable gene IDs and aligned coordinates."""
    import scanpy as sc

    adata = sc.read_10x_h5(matrix_path, genome=None, gex_only=True)
    if "gene_ids" not in adata.var:
        raise ValueError("10x matrix does not contain gene_ids")
    adata.var["gene_symbol"] = adata.var_names.astype(str)
    adata.var_names = adata.var["gene_ids"].astype(str)
    adata.var_names = pd.Index(adata.var_names).str.replace(r"\.\d+$", "", regex=True)
    if not adata.var_names.is_unique:
        raise ValueError("stable gene IDs are not unique after version stripping")
    if not adata.obs_names.is_unique:
        raise ValueError("10x matrix contains duplicate barcodes")

    positions = read_tissue_positions(positions_path)
    missing = adata.obs_names.difference(positions.index)
    if len(missing):
        raise ValueError(f"{len(missing)} matrix barcodes are missing coordinates")
    positions = positions.loc[adata.obs_names]
    if in_tissue_only:
        keep = positions["in_tissue"].astype(int).eq(1).to_numpy()
        adata = adata[keep].copy()
        positions = positions.iloc[np.flatnonzero(keep)]
    else:
        adata = adata.copy()

    for column in POSITION_COLUMNS[1:]:
        adata.obs[column] = positions[column].to_numpy()
    adata.obsm["spatial"] = positions[
        ["pxl_col_in_fullres", "pxl_row_in_fullres"]
    ].to_numpy(dtype=np.float64)
    adata.X = sparse.csr_matrix(adata.X, dtype=np.float32)
    adata.layers["counts"] = adata.X.copy()
    adata.uns["counts_source"] = "10x_filtered_feature_bc_matrix"
    validate_spatial_counts(adata)
    return adata


def validate_spatial_counts(adata) -> dict[str, object]:
    """Validate invariants needed by all baseline backends."""
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("spatial observation and variable names must be unique")
    if "counts" not in adata.layers or "spatial" not in adata.obsm:
        raise ValueError("spatial object requires counts layer and spatial coordinates")
    counts = adata.layers["counts"]
    values = counts.data if sparse.issparse(counts) else np.asarray(counts).ravel()
    coords = np.asarray(adata.obsm["spatial"])
    if values.size and (not np.isfinite(values).all() or np.min(values) < 0):
        raise ValueError("spatial counts contain NaN, infinity, or negative values")
    if values.size and not np.allclose(values, np.rint(values), atol=1e-6):
        raise ValueError("spatial counts are not integer-like raw counts")
    if coords.shape != (adata.n_obs, 2) or not np.isfinite(coords).all():
        raise ValueError("spatial coordinates are missing, misaligned, or non-finite")
    return {
        "status": "PASS",
        "spots": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "unique_spots": bool(adata.obs_names.is_unique),
        "unique_genes": bool(adata.var_names.is_unique),
        "raw_nonnegative_integer_counts": True,
        "coordinates_aligned": True,
    }
