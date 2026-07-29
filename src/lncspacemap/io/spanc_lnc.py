"""SPanC-Lnc data audit and reference construction.

The raw source directory is treated as read-only. Large matrices are read
incrementally and all generated H5AD files are written outside git.
"""

from __future__ import annotations

import csv
import json
import logging
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import scipy.sparse as sp

LOG = logging.getLogger("lncspacemap.spanc_lnc")
_CUTAR_RE = re.compile(r"^(?:cu|u)?TAR[\w.-]*$", re.IGNORECASE)
_BARCODE_RE = re.compile(r"([ACGTN]{16})", re.IGNORECASE)


@dataclass(frozen=True)
class SamplePair:
    sample_id: str
    gene_file: str
    cutar_file: str
    subtype: str
    treatment: str


@dataclass
class TextMatrixInfo:
    path: str
    delimiter: str
    orientation: str
    n_rows: int
    n_columns: int
    row_id_name: str
    first_ids: list[str]


def load_sample_pairs(config_path: Path) -> list[SamplePair]:
    import yaml

    cfg = yaml.safe_load(config_path.read_text())
    return [SamplePair(**item) for item in cfg["sample_pairs"]]


def _delimiter(path: Path) -> str:
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        sample = handle.read(65536)
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,;").delimiter
    except csv.Error:
        return "\t"


def inspect_text_matrix(path: Path, count_rows: bool = True) -> TextMatrixInfo:
    sep = _delimiter(path)
    preview = pd.read_csv(path, sep=sep, nrows=8, dtype=str)
    if preview.shape[1] < 2:
        raise ValueError(f"{path.name}: expected at least two columns")
    column_ids = [str(x) for x in preview.columns[1:]]
    row_ids = preview.iloc[:, 0].astype(str).tolist()
    row_score = np.mean([bool(_CUTAR_RE.match(x)) for x in row_ids])
    col_score = np.mean([bool(_CUTAR_RE.match(x)) for x in column_ids])
    if row_score >= 0.5 and row_score > col_score:
        orientation = "features_by_cells"
    elif col_score >= 0.5 and col_score > row_score:
        orientation = "cells_by_features"
    else:
        raise ValueError(
            f"{path.name}: cannot determine matrix orientation "
            f"(row cuTAR fraction={row_score:.3f}, column fraction={col_score:.3f})"
        )
    n_rows = -1
    if count_rows:
        with path.open("rb") as handle:
            n_rows = max(sum(1 for _ in handle) - 1, 0)
    return TextMatrixInfo(
        path=path.name,
        delimiter="tab" if sep == "\t" else repr(sep),
        orientation=orientation,
        n_rows=n_rows,
        n_columns=preview.shape[1] - 1,
        row_id_name=str(preview.columns[0]),
        first_ids=row_ids[:3],
    )


def read_cutar_matrix(path: Path, chunk_rows: int = 64) -> ad.AnnData:
    """Read a dense text matrix incrementally and return sparse cells x cuTARs."""
    import anndata as ad

    info = inspect_text_matrix(path, count_rows=False)
    sep = "\t" if info.delimiter == "tab" else info.delimiter.strip("'")
    blocks: list[sp.csr_matrix] = []
    row_ids: list[str] = []
    column_ids: list[str] | None = None
    for chunk in pd.read_csv(
        path, sep=sep, chunksize=chunk_rows, index_col=0, dtype=np.float32
    ):
        ids = chunk.index.astype(str).tolist()
        values = chunk.to_numpy(dtype=np.float32, copy=False)
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f"{path.name}: counts must be finite and non-negative")
        if not np.allclose(values, np.rint(values), atol=1e-6):
            raise ValueError(f"{path.name}: non-integer values found; raw counts unresolved")
        blocks.append(sp.csr_matrix(values))
        row_ids.extend(ids)
        if column_ids is None:
            column_ids = [str(x) for x in chunk.columns]
    if not blocks or column_ids is None:
        raise ValueError(f"{path.name}: empty matrix")
    matrix = sp.vstack(blocks, format="csr")
    if info.orientation == "features_by_cells":
        matrix = matrix.T.tocsr()
        obs_names, var_names = column_ids, row_ids
    else:
        obs_names, var_names = row_ids, column_ids
    if len(set(obs_names)) != len(obs_names):
        raise ValueError(f"{path.name}: duplicate cell barcodes")
    if len(set(var_names)) != len(var_names):
        raise ValueError(f"{path.name}: duplicate cuTAR identifiers")
    out = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=pd.Index(obs_names, name="cell_id")),
        var=pd.DataFrame(index=pd.Index(var_names, name="feature_id")),
    )
    out.var["feature_type"] = "cuTAR"
    out.layers["counts"] = out.X.copy()
    return out


def read_gene_matrix(path: Path) -> ad.AnnData:
    import scanpy as sc

    out = sc.read_10x_h5(path, gex_only=True)
    out.X = out.X.tocsr().astype(np.float32)
    if out.X.nnz and (not np.isfinite(out.X.data).all() or (out.X.data < 0).any()):
        raise ValueError(f"{path.name}: invalid gene counts")
    symbols = out.var_names.astype(str)
    gene_ids = out.var.get("gene_ids", pd.Series(symbols, index=out.var_names)).astype(str)
    if gene_ids.duplicated().any():
        raise ValueError(f"{path.name}: duplicate gene IDs")
    out.var["gene_symbol"] = symbols.to_numpy()
    out.var_names = pd.Index(gene_ids.to_numpy(), name="feature_id")
    out.var["feature_type"] = "gene"
    out.layers["counts"] = out.X.copy()
    return out


def _barcode_keys(values: Iterable[str], strategy: str) -> list[str]:
    keys = []
    for value in values:
        text = str(value).strip()
        if strategy == "exact":
            key = text
        elif strategy == "strip_suffix":
            key = re.sub(r"-\d+$", "", text)
        elif strategy == "sequence16":
            matches = _BARCODE_RE.findall(text)
            key = matches[-1].upper() if matches else ""
        else:
            raise ValueError(strategy)
        keys.append(key)
    return keys


def match_barcodes(
    gene_barcodes: Iterable[str], cutar_barcodes: Iterable[str]
) -> tuple[str, list[int], list[int]]:
    gene_values = list(map(str, gene_barcodes))
    cutar_values = list(map(str, cutar_barcodes))
    candidates = []
    for strategy in ("exact", "strip_suffix", "sequence16"):
        left = _barcode_keys(gene_values, strategy)
        right = _barcode_keys(cutar_values, strategy)
        if "" in left or "" in right:
            continue
        if len(set(left)) != len(left) or len(set(right)) != len(right):
            continue
        right_index = {key: i for i, key in enumerate(right)}
        pairs = [(i, right_index[key]) for i, key in enumerate(left) if key in right_index]
        candidates.append((len(pairs), strategy, pairs))
    if not candidates:
        raise ValueError("no unambiguous barcode normalization strategy")
    overlap, strategy, pairs = max(candidates, key=lambda item: item[0])
    if overlap == 0:
        raise ValueError("gene and cuTAR matrices have zero barcode overlap")
    LOG.info("barcode strategy=%s overlap=%d", strategy, overlap)
    return strategy, [x[0] for x in pairs], [x[1] for x in pairs]


def combine_sample(pair: SamplePair, data_dir: Path) -> tuple[ad.AnnData, dict]:
    import anndata as ad

    genes = read_gene_matrix(data_dir / pair.gene_file)
    cutars = read_cutar_matrix(data_dir / pair.cutar_file)
    original_gene_cells = genes.n_obs
    original_cutar_cells = cutars.n_obs
    strategy, gene_idx, cutar_idx = match_barcodes(genes.obs_names, cutars.obs_names)
    overlap = len(gene_idx)
    min_fraction = min(overlap / genes.n_obs, overlap / cutars.n_obs)
    if min_fraction < 0.8:
        raise ValueError(
            f"{pair.sample_id}: barcode overlap fraction {min_fraction:.3f} < 0.8"
        )
    genes = genes[gene_idx].copy()
    cutars = cutars[cutar_idx].copy()
    cutars.obs_names = genes.obs_names.copy()
    combined = ad.concat(
        [genes, cutars],
        axis=1,
        join="outer",
        merge="first",
        label="feature_source",
        keys=["gene", "cuTAR"],
        index_unique=None,
    )
    combined.obs["sample"] = pair.sample_id
    combined.obs["disease"] = "melanoma"
    combined.obs["subtype"] = pair.subtype
    combined.obs["treatment"] = pair.treatment
    combined.obs["source_accession"] = pair.gene_file.split("_")[0]
    combined.uns["lncspacemap_preprocessing"] = {
        "schema_version": "0.1",
        "counts_source": "layers/counts",
        "barcode_strategy": strategy,
        "gene_file": pair.gene_file,
        "cutar_file": pair.cutar_file,
    }
    metrics = {
        "sample_id": pair.sample_id,
        "gene_cells": original_gene_cells,
        "cutar_cells": original_cutar_cells,
        "matched_cells": overlap,
        "barcode_strategy": strategy,
        "n_genes": genes.n_vars,
        "n_cutars": cutars.n_vars,
        "gene_nnz": int(genes.X.nnz),
        "cutar_nnz": int(cutars.X.nnz),
    }
    return combined, metrics


def audit_10x_h5(path: Path, display_name: str | None = None) -> dict:
    import h5py

    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        shape = tuple(int(x) for x in matrix["shape"][:])
        return {
            "file": display_name or path.name,
            "kind": "10x_h5",
            "observations": shape[1],
            "features": shape[0],
            "bytes": path.stat().st_size,
            "status": "PASS",
            "detail": "",
        }


def audit_h5ad(path: Path) -> dict:
    import anndata as ad

    try:
        obj = ad.read_h5ad(path, backed="r")
        detail = json.dumps(
            {
                "layers": list(obj.layers.keys()),
                "obsm": list(obj.obsm.keys()),
                "obs_columns": list(obj.obs.columns),
                "var_columns": list(obj.var.columns),
            },
            separators=(",", ":"),
        )
        result = {
            "file": path.name,
            "kind": "h5ad",
            "observations": obj.n_obs,
            "features": obj.n_vars,
            "bytes": path.stat().st_size,
            "status": "PASS",
            "detail": detail,
        }
        obj.file.close()
        return result
    except Exception as exc:  # audit must report all files
        return {
            "file": path.name,
            "kind": "h5ad",
            "observations": "",
            "features": "",
            "bytes": path.stat().st_size,
            "status": "FAIL",
            "detail": str(exc),
        }


def run_audit(data_dir: Path, review_dir: Path, pairs: list[SamplePair]) -> None:
    review_dir.mkdir(parents=True, exist_ok=True)
    (review_dir / "metrics").mkdir(exist_ok=True)
    (review_dir / "manifests").mkdir(exist_ok=True)
    rows = [audit_10x_h5(p) for p in sorted(data_dir.glob("GSM*.h5"))]
    for folder_name in ("MelD", "MelDN"):
        folder = data_dir / folder_name
        matrix = folder / "filtered_feature_bc_matrix.h5"
        if matrix.is_file():
            rows.append(audit_10x_h5(matrix, f"{folder_name}/{matrix.name}"))
        else:
            rows.append(
                {
                    "file": f"{folder_name}/filtered_feature_bc_matrix.h5",
                    "kind": "10x_spatial_h5",
                    "observations": "",
                    "features": "",
                    "bytes": "",
                    "status": "FAIL",
                    "detail": "missing file",
                }
            )
        position_candidates = [
            folder / "tissue_positions_list.csv",
            folder / "spatial/tissue_positions_list.csv",
            folder / "spatial/tissue_positions.csv",
        ]
        if not any(path.is_file() for path in position_candidates):
            rows.append(
                {
                    "file": f"{folder_name}/tissue_positions",
                    "kind": "spatial_coordinates",
                    "observations": "",
                    "features": "",
                    "bytes": "",
                    "status": "FAIL",
                    "detail": "no tissue position file found",
                }
            )
    rows.extend(audit_h5ad(p) for p in sorted(data_dir.glob("*.h5ad")))
    pair_rows = []
    for pair in pairs:
        try:
            info = inspect_text_matrix(data_dir / pair.cutar_file)
            pair_rows.append({"sample_id": pair.sample_id, **asdict(info), "status": "PASS"})
        except Exception as exc:
            pair_rows.append(
                {
                    "sample_id": pair.sample_id,
                    "path": pair.cutar_file,
                    "status": "FAIL",
                    "error": str(exc),
                }
            )
    pd.DataFrame(rows).to_csv(
        review_dir / "metrics/spanc_lnc_matrix_audit.tsv", sep="\t", index=False
    )
    pd.DataFrame(pair_rows).to_csv(
        review_dir / "metrics/spanc_lnc_cutar_audit.tsv", sep="\t", index=False
    )
    manifest = [
        {
            "file": p.name,
            "bytes": p.stat().st_size,
            "suffix": "".join(p.suffixes),
        }
        for p in sorted(data_dir.iterdir())
        if p.is_file()
    ]
    pd.DataFrame(manifest).to_csv(
        review_dir / "manifests/spanc_lnc_file_manifest.tsv", sep="\t", index=False
    )


def build_references(
    data_dir: Path, output_dir: Path, review_dir: Path, pairs: list[SamplePair]
) -> None:
    import anndata as ad

    output_dir.mkdir(parents=True, exist_ok=True)
    metrics = []
    objects = []
    for pair in pairs:
        LOG.info("building %s", pair.sample_id)
        obj, row = combine_sample(pair, data_dir)
        obj.write_h5ad(output_dir / f"{pair.sample_id}.gene_cutar.h5ad", compression="lzf")
        objects.append(obj)
        metrics.append(row)
    combined = ad.concat(
        objects,
        axis=0,
        join="outer",
        merge="same",
        label="sample_batch",
        keys=[p.sample_id for p in pairs],
        index_unique=":",
        fill_value=0,
    )
    combined.layers["counts"] = combined.X.copy()
    combined.write_h5ad(
        output_dir / "acral_melanoma_gene_cutar_combined.h5ad", compression="lzf"
    )
    review_dir.mkdir(parents=True, exist_ok=True)
    (review_dir / "metrics").mkdir(exist_ok=True)
    pd.DataFrame(metrics).to_csv(
        review_dir / "metrics/spanc_lnc_reference_build.tsv", sep="\t", index=False
    )
