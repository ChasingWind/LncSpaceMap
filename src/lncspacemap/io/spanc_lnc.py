"""SPanC-Lnc data audit and reference construction.

The raw source directory is treated as read-only. Large matrices are read
incrementally and all generated H5AD files are written outside git.
"""

from __future__ import annotations

import csv
import json
import logging
import re
from dataclasses import dataclass
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
    header_style: str
    n_rows: int
    n_columns: int
    row_id_name: str
    first_ids: list[str]
    row_cutar_fraction: float
    column_cutar_fraction: float
    row_barcode_fraction: float
    column_barcode_fraction: float


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


def _fraction(values: Iterable[str], pattern: re.Pattern[str]) -> float:
    values = [str(value).strip() for value in values]
    return float(np.mean([bool(pattern.match(value)) for value in values])) if values else 0.0


def _text_matrix_layout(path: Path) -> dict:
    """Resolve common TSV layouts without allowing pandas to infer an index."""
    sep = _delimiter(path)
    with path.open("rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter=sep)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"{path.name}: empty matrix") from exc
        preview_rows = []
        for _ in range(3):
            try:
                preview_rows.append(next(reader))
            except StopIteration:
                break
    if len(header) < 2 or not preview_rows:
        raise ValueError(f"{path.name}: expected at least two columns")

    widths = {len(row) for row in preview_rows}
    row_ids = [row[0].strip() for row in preview_rows if row]
    implicit_index = widths == {len(header) + 1}
    standard_width = widths == {len(header)}
    if not implicit_index and not standard_width:
        raise ValueError(
            f"{path.name}: inconsistent field counts "
            f"(header={len(header)}, preview={sorted(widths)})"
        )

    if implicit_index:
        column_ids = [str(value).strip() for value in header]
        header_style = "implicit_row_index"
    else:
        column_ids = [str(value).strip() for value in header[1:]]
        header_style = "explicit_row_index"

    row_cutar = _fraction(row_ids, _CUTAR_RE)
    column_cutar = _fraction(column_ids, _CUTAR_RE)
    row_barcode = _fraction(row_ids, _BARCODE_RE)
    column_barcode = _fraction(column_ids, _BARCODE_RE)

    if row_cutar >= 0.5 and column_barcode >= 0.5:
        orientation = "features_by_cells"
        cell_ids = column_ids
    elif row_barcode >= 0.5 and column_cutar >= 0.5 and standard_width:
        orientation = "cells_by_features"
        cell_ids = row_ids
    else:
        raise ValueError(
            f"{path.name}: cannot determine matrix orientation "
            f"(row cuTAR={row_cutar:.3f}, column cuTAR={column_cutar:.3f}, "
            f"row barcode={row_barcode:.3f}, column barcode={column_barcode:.3f}, "
            f"header_style={header_style})"
        )

    numeric_preview = [
        value
        for row in preview_rows
        for value in (row[1:] if standard_width else row[1:])
    ]
    try:
        numeric = np.asarray(numeric_preview, dtype=np.float64)
    except ValueError as exc:
        raise ValueError(f"{path.name}: non-numeric count found in preview") from exc
    if not np.isfinite(numeric).all() or (numeric < 0).any():
        raise ValueError(f"{path.name}: counts must be finite and non-negative")
    if not np.allclose(numeric, np.rint(numeric), atol=1e-6):
        raise ValueError(f"{path.name}: non-integer value found in count preview")

    return {
        "sep": sep,
        "orientation": orientation,
        "header_style": header_style,
        "row_ids": row_ids,
        "column_ids": column_ids,
        "cell_ids": cell_ids,
        "row_cutar_fraction": row_cutar,
        "column_cutar_fraction": column_cutar,
        "row_barcode_fraction": row_barcode,
        "column_barcode_fraction": column_barcode,
    }


def inspect_text_matrix(path: Path, count_rows: bool = True) -> TextMatrixInfo:
    layout = _text_matrix_layout(path)
    n_rows = -1
    if count_rows:
        with path.open("rb") as handle:
            n_rows = max(sum(1 for _ in handle) - 1, 0)
    return TextMatrixInfo(
        path=path.name,
        delimiter="tab" if layout["sep"] == "\t" else repr(layout["sep"]),
        orientation=layout["orientation"],
        header_style=layout["header_style"],
        n_rows=n_rows,
        n_columns=len(layout["column_ids"]),
        row_id_name=(
            "implicit_index"
            if layout["header_style"] == "implicit_row_index"
            else "first_column"
        ),
        first_ids=layout["row_ids"][:3],
        row_cutar_fraction=layout["row_cutar_fraction"],
        column_cutar_fraction=layout["column_cutar_fraction"],
        row_barcode_fraction=layout["row_barcode_fraction"],
        column_barcode_fraction=layout["column_barcode_fraction"],
    )


def read_cutar_matrix(path: Path, chunk_rows: int = 64):
    """Read a dense text matrix incrementally and return sparse cells x cuTARs."""
    import anndata as ad

    layout = _text_matrix_layout(path)
    sep = layout["sep"]
    blocks: list[sp.csr_matrix] = []
    row_ids: list[str] = []
    column_ids: list[str] | None = None
    read_options = {
        "sep": sep,
        "chunksize": chunk_rows,
        "index_col": 0,
        "dtype": np.float32,
    }
    if layout["header_style"] == "implicit_row_index":
        read_options.update(
            {
                "header": None,
                "skiprows": 1,
                "names": ["feature_id", *layout["column_ids"]],
            }
        )
    for chunk in pd.read_csv(path, **read_options):
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
    if layout["orientation"] == "features_by_cells":
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


def read_10x_barcodes(path: Path) -> list[str]:
    import h5py

    with h5py.File(path, "r") as handle:
        values = handle["matrix/barcodes"][:]
    return [
        value.decode("utf-8") if isinstance(value, (bytes, np.bytes_)) else str(value)
        for value in values
    ]


def load_bed_feature_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path.name}: BED line {line_number} has <4 columns")
            ids.add(fields[3])
    if not ids:
        raise ValueError(f"{path.name}: no feature IDs found in BED column 4")
    return ids


def scan_cutar_feature_ids(path: Path, bed_ids: set[str]) -> dict:
    sep = _delimiter(path)
    seen: set[str] = set()
    duplicates = 0
    bed_overlap = 0
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        next(handle, None)
        for line in handle:
            feature_id = line.split(sep, 1)[0].strip()
            if feature_id in seen:
                duplicates += 1
            else:
                seen.add(feature_id)
                bed_overlap += int(feature_id in bed_ids)
    return {
        "n_cutars": len(seen),
        "duplicate_cutar_ids": duplicates,
        "bed_overlap": bed_overlap,
        "bed_overlap_fraction": bed_overlap / len(seen) if seen else 0.0,
    }


def audit_cutar_pair(
    pair: SamplePair, data_dir: Path, bed_ids: set[str], min_overlap: float = 0.8
) -> dict:
    gene_path = data_dir / pair.gene_file
    cutar_path = data_dir / pair.cutar_file
    info = inspect_text_matrix(cutar_path)
    if info.orientation != "features_by_cells":
        raise ValueError(
            f"{pair.sample_id}: expected features_by_cells, got {info.orientation}"
        )
    layout = _text_matrix_layout(cutar_path)
    gene_barcodes = read_10x_barcodes(gene_path)
    cutar_barcodes = layout["cell_ids"]
    strategy, gene_idx, cutar_idx = match_barcodes(gene_barcodes, cutar_barcodes)
    matched = len(gene_idx)
    gene_fraction = matched / len(gene_barcodes)
    cutar_fraction = matched / len(cutar_barcodes)
    feature_stats = scan_cutar_feature_ids(cutar_path, bed_ids)
    failures = []
    if gene_fraction < min_overlap or cutar_fraction < min_overlap:
        failures.append(
            f"barcode overlap below {min_overlap:.2f} "
            f"(gene={gene_fraction:.3f}, cuTAR={cutar_fraction:.3f})"
        )
    if feature_stats["duplicate_cutar_ids"]:
        failures.append(f"{feature_stats['duplicate_cutar_ids']} duplicate cuTAR IDs")
    if feature_stats["bed_overlap_fraction"] < 0.95:
        failures.append(
            "BED overlap below 0.95 "
            f"({feature_stats['bed_overlap_fraction']:.3f})"
        )
    return {
        "sample_id": pair.sample_id,
        "gene_file": pair.gene_file,
        "cutar_file": pair.cutar_file,
        "status": "FAIL" if failures else "PASS",
        "orientation": info.orientation,
        "header_style": info.header_style,
        "n_cutars": feature_stats["n_cutars"],
        "n_cutar_cells": len(cutar_barcodes),
        "n_gene_cells": len(gene_barcodes),
        "barcode_strategy": strategy,
        "matched_cells": matched,
        "gene_match_fraction": gene_fraction,
        "cutar_match_fraction": cutar_fraction,
        **feature_stats,
        "row_cutar_fraction": info.row_cutar_fraction,
        "column_barcode_fraction": info.column_barcode_fraction,
        "error": "; ".join(failures),
    }


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
    LOG.info("starting SPanC-Lnc enhanced audit: %s", data_dir)
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
    bed_path = data_dir / "00_cuTARs.bed"
    if not bed_path.is_file():
        raise FileNotFoundError(bed_path)
    bed_ids = load_bed_feature_ids(bed_path)
    LOG.info("loaded %d unique cuTAR IDs from %s", len(bed_ids), bed_path.name)
    pair_rows = []
    for pair in pairs:
        try:
            result = audit_cutar_pair(pair, data_dir, bed_ids)
            pair_rows.append(result)
            LOG.info(
                "audit %s status=%s matched=%s gene_fraction=%.3f "
                "cutar_fraction=%.3f bed_fraction=%.3f",
                pair.sample_id,
                result["status"],
                result["matched_cells"],
                result["gene_match_fraction"],
                result["cutar_match_fraction"],
                result["bed_overlap_fraction"],
            )
        except Exception as exc:
            LOG.exception("audit %s failed", pair.sample_id)
            pair_rows.append(
                {
                    "sample_id": pair.sample_id,
                    "gene_file": pair.gene_file,
                    "cutar_file": pair.cutar_file,
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
    pair_pass = bool(pair_rows) and all(row["status"] == "PASS" for row in pair_rows)
    matrix_pass = bool(rows) and all(row["status"] == "PASS" for row in rows)
    marker = (
        "PASS_SPANCLNC_ENHANCED_AUDIT_READY_FOR_BUILD"
        if pair_pass and matrix_pass
        else "HOLD_BUILD_ENHANCED_AUDIT_FAILED"
    )
    LOG.info(
        "%s pair_pass=%s matrix_pass=%s pair_count=%d matrix_checks=%d",
        marker,
        pair_pass,
        matrix_pass,
        len(pair_rows),
        len(rows),
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
