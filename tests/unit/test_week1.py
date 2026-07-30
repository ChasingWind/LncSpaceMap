from pathlib import Path

import numpy as np
import pandas as pd

from lncspacemap.io.gtf import read_gene_gtf, strip_ensembl_version
from lncspacemap.preprocessing.eligibility import (
    EligibilityPolicy,
    annotate_gene_features,
    build_target_catalog,
    select_masked_gene_proxies,
)


def test_read_gene_gtf_and_strip_version(tmp_path: Path):
    path = tmp_path / "genes.gtf"
    path.write_text(
        '##description: test\n'
        'chr1\tHAVANA\tgene\t11\t20\t.\t+\t.\t'
        'gene_id "ENSG000001.4"; gene_name "GENE1"; gene_type "protein_coding";\n'
    )
    table = read_gene_gtf(path)
    assert strip_ensembl_version("ENSG000001.4") == "ENSG000001"
    assert table.index.tolist() == ["ENSG000001"]
    assert table.loc["ENSG000001", "start"] == 10
    assert table.loc["ENSG000001", "gene_type"] == "protein_coding"


def test_target_eligibility_distinguishes_supported_and_exploratory():
    feature_qc = pd.DataFrame(
        {
            "feature_type": ["cuTAR", "cuTAR", "cuTAR"],
            "total_counts": [100, 30, 0],
            "detected_cells": [30, 12, 0],
            "quantified_cells": [1000, 1000, 1000],
            "quantified_sample_count": [3, 1, 1],
            "detected_cell_fraction": [0.03, 0.012, 0.0],
            "chrom": ["chr1", "chr1", "chr1"],
            "start": [1, 2, 3],
            "end": [10, 20, 30],
            "strand": ["+", "+", "-"],
        },
        index=["cuTAR1", "cuTAR2", "cuTAR3"],
    )
    catalog = build_target_catalog(feature_qc, EligibilityPolicy())
    assert catalog["eligibility"].tolist() == [
        "eligible",
        "exploratory",
        "insufficient_reference_signal",
    ]
    assert catalog.loc["cuTAR2", "eligibility_reason"] == "low_sample_support"


def test_masked_proxy_folds_use_only_annotated_protein_coding_genes():
    gene_ids = [f"ENSG{i:03d}" for i in range(20)]
    target_ids = [f"cuTAR{i}" for i in range(5)]
    index = gene_ids + target_ids
    feature_qc = pd.DataFrame(
        {
            "feature_type": ["gene"] * 20 + ["cuTAR"] * 5,
            "total_counts": np.arange(20, 45),
            "detected_cells": np.arange(10, 35),
            "quantified_cells": [1000] * 25,
            "quantified_sample_count": [6] * 20 + [3] * 5,
            "detected_cell_fraction": np.arange(10, 35) / 1000,
            "chrom": [np.nan] * 20 + ["chr1"] * 5,
            "start": [np.nan] * 20 + list(range(5)),
            "end": [np.nan] * 20 + list(range(10, 15)),
            "strand": [np.nan] * 20 + ["+"] * 5,
        },
        index=index,
    )
    gene_table = pd.DataFrame(
        {
            "gene_id_versioned": [f"{gene}.1" for gene in gene_ids],
            "gene_name": gene_ids,
            "gene_type": ["protein_coding"] * 18 + ["lncRNA"] * 2,
            "chrom": ["chr1"] * 20,
            "start": range(20),
            "end": range(10, 30),
            "strand": ["+"] * 20,
            "annotation_source": ["HAVANA"] * 20,
        },
        index=gene_ids,
    )
    annotated = annotate_gene_features(feature_qc, gene_table)
    targets = build_target_catalog(annotated, EligibilityPolicy())
    proxies = select_masked_gene_proxies(annotated, targets, n_genes=10, n_folds=5)
    assert len(proxies) == 10
    assert proxies["gene_type"].eq("protein_coding").all()
    assert set(proxies["fold"]) == {0, 1, 2, 3, 4}
    assert not set(proxies.index) & set(targets.index)
