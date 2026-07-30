from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from lncspacemap.io.spanc_lnc import (
    _assign_cutar_tiers,
    _feature_catalog,
    SamplePair,
    inspect_text_matrix,
    match_barcodes,
    read_cutar_matrix,
)


def test_detect_features_by_cells(tmp_path: Path):
    path = tmp_path / "matrix.tsv"
    pd.DataFrame(
        {
            "feature": ["cuTAR1", "cuTAR2"],
            "AAAAAAAAAAAAAAAA-1": [1, 0],
            "CCCCCCCCCCCCCCCC-1": [0, 2],
        }
    ).to_csv(path, sep="\t", index=False)
    info = inspect_text_matrix(path)
    assert info.orientation == "features_by_cells"
    assert info.header_style == "explicit_row_index"
    assert info.n_rows == 2


def test_detect_implicit_row_index_header(tmp_path: Path):
    path = tmp_path / "uq_cutar.tsv"
    path.write_text(
        "AAAAAAAAAAAAAAAA.1\tCCCCCCCCCCCCCCCC.1\n"
        "cuTAR1\t1\t0\n"
        "cuTAR2\t0\t2\n"
    )
    info = inspect_text_matrix(path)
    assert info.orientation == "features_by_cells"
    assert info.header_style == "implicit_row_index"
    assert info.n_rows == 2
    assert info.n_columns == 2
    assert info.first_ids == ["cuTAR1", "cuTAR2"]
    assert info.row_cutar_fraction == 1.0
    assert info.column_barcode_fraction == 1.0


def test_read_implicit_row_index_count_matrix(tmp_path: Path):
    path = tmp_path / "uq_cutar.tsv"
    path.write_text(
        "AAAAAAAAAAAAAAAA.1\tCCCCCCCCCCCCCCCC.1\n"
        "cuTAR1\t1\t0\n"
        "cuTAR2\t0\t2\n"
    )
    obj = read_cutar_matrix(path, chunk_rows=1)
    assert obj.shape == (2, 2)
    assert obj.obs_names.tolist() == [
        "AAAAAAAAAAAAAAAA.1",
        "CCCCCCCCCCCCCCCC.1",
    ]
    assert obj.var_names.tolist() == ["cuTAR1", "cuTAR2"]
    np.testing.assert_array_equal(obj.X.toarray(), [[1, 0], [0, 2]])
    np.testing.assert_array_equal(obj.layers["counts"].toarray(), [[1, 0], [0, 2]])


def test_barcode_suffix_matching():
    strategy, left, right = match_barcodes(
        ["AAAAAAAAAAAAAAAA-1", "CCCCCCCCCCCCCCCC-1"],
        ["AAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCC"],
    )
    assert strategy == "strip_suffix"
    assert left == [0, 1]
    assert right == [0, 1]


def test_barcode_dot_suffix_matches_10x_dash_suffix():
    strategy, left, right = match_barcodes(
        ["AAAAAAAAAAAAAAAA-1", "CCCCCCCCCCCCCCCC-1"],
        ["AAAAAAAAAAAAAAAA.1", "CCCCCCCCCCCCCCCC.1"],
    )
    assert strategy == "sequence16"
    assert left == [0, 1]
    assert right == [0, 1]


def test_cutar_tiers_require_cross_sample_support():
    tiers = _assign_cutar_tiers(
        detected_cells=np.array([0, 9, 10, 30]),
        quantified_cells=np.array([1000, 1000, 1000, 1000]),
        quantified_samples=np.array([6, 2, 2, 3]),
    )
    assert tiers.tolist() == ["excluded", "all_detected", "extended", "core"]


def test_feature_catalog_preserves_type_coordinates_and_quantification(tmp_path: Path):
    bed = tmp_path / "cutars.bed"
    bed.write_text(
        "chr1\t10\t20\tcuTAR1\t.\t+\n"
        "chr2\t30\t40\tcuTAR2\t.\t-\n"
    )
    first = SimpleNamespace(
        var_names=pd.Index(["ENSG1", "cuTAR1"]),
        var=pd.DataFrame(
            {
                "feature_type": ["gene", "cuTAR"],
                "gene_symbol": ["GENE1", np.nan],
            },
            index=["ENSG1", "cuTAR1"],
        ),
    )
    second = SimpleNamespace(
        var_names=pd.Index(["ENSG1", "cuTAR2"]),
        var=pd.DataFrame(
            {
                "feature_type": ["gene", "cuTAR"],
                "gene_symbol": ["GENE1", np.nan],
            },
            index=["ENSG1", "cuTAR2"],
        ),
    )
    pairs = [
        SamplePair("S1", "g1.h5", "u1.tsv", "acral", "untreated"),
        SamplePair("S2", "g2.h5", "u2.tsv", "acral", "untreated"),
    ]
    catalog = _feature_catalog([first, second], pairs, bed)
    assert catalog.loc["ENSG1", "feature_type"] == "gene"
    assert catalog.loc["ENSG1", "quantified_sample_count"] == 2
    assert catalog.loc["cuTAR1", "quantified_sample_count"] == 1
    assert bool(catalog.loc["cuTAR1", "quantified_S1"])
    assert not bool(catalog.loc["cuTAR1", "quantified_S2"])
    assert catalog.loc["cuTAR1", "chrom"] == "chr1"
    assert catalog.loc["cuTAR2", "strand"] == "-"
