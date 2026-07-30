from pathlib import Path

import pandas as pd

from lncspacemap.io.spanc_lnc import inspect_text_matrix, match_barcodes


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
