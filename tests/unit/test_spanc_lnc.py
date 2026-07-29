from pathlib import Path

import pandas as pd

from lncspacemap.io.spanc_lnc import inspect_text_matrix, match_barcodes


def test_detect_features_by_cells(tmp_path: Path):
    path = tmp_path / "matrix.tsv"
    pd.DataFrame(
        {"feature": ["cuTAR1", "cuTAR2"], "AAAC-1": [1, 0], "CCGT-1": [0, 2]}
    ).to_csv(path, sep="\t", index=False)
    info = inspect_text_matrix(path)
    assert info.orientation == "features_by_cells"
    assert info.n_rows == 2


def test_barcode_suffix_matching():
    strategy, left, right = match_barcodes(
        ["AAAAAAAAAAAAAAAA-1", "CCCCCCCCCCCCCCCC-1"],
        ["AAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCC"],
    )
    assert strategy == "strip_suffix"
    assert left == [0, 1]
    assert right == [0, 1]
