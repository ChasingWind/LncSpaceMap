"""Minimal, deterministic GENCODE gene-table reader."""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd

_ATTRIBUTE_RE = re.compile(r'(\S+)\s+"([^"]*)";')


def _open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("rt")


def strip_ensembl_version(value: str) -> str:
    return str(value).split(".", 1)[0]


def read_gene_gtf(path: Path) -> pd.DataFrame:
    """Read gene records from a GTF without loading transcript/exon rows."""
    rows = []
    with _open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attributes = dict(_ATTRIBUTE_RE.findall(fields[8]))
            versioned_id = attributes.get("gene_id")
            if not versioned_id:
                continue
            rows.append(
                {
                    "gene_id": strip_ensembl_version(versioned_id),
                    "gene_id_versioned": versioned_id,
                    "gene_name": attributes.get("gene_name", ""),
                    "gene_type": attributes.get(
                        "gene_type", attributes.get("gene_biotype", "")
                    ),
                    "chrom": fields[0],
                    "start": int(fields[3]) - 1,
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "annotation_source": fields[1],
                }
            )
    if not rows:
        raise ValueError(f"{path}: no GTF gene records found")
    table = pd.DataFrame(rows)
    if table["gene_id"].duplicated().any():
        duplicates = int(table["gene_id"].duplicated().sum())
        raise ValueError(f"{path.name}: {duplicates} duplicate canonical gene IDs")
    return table.set_index("gene_id").sort_index()
