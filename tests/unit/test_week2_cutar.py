import numpy as np
import pandas as pd
from scipy import sparse

from lncspacemap.pipeline.week2_cutar import (
    _admission_checks,
    _sha1,
    build_target_catalog,
)


def test_sha1_source_contract(tmp_path):
    path = tmp_path / "matrix.txt"
    path.write_bytes(b"released-matrix\n")
    assert _sha1(path) == "a34764384d224c765f3886570e89c578dc402af0"


def test_cross_sample_reference_overlap_uses_target_count_not_fraction():
    checks = _admission_checks(
        barcode_fraction=1.0,
        spatial_reference_overlap_count=362,
        spatial_bed_fraction=1.0,
        coordinate_valid=True,
        frozen_count=362,
        contract_cfg={
            "min_spatial_barcode_overlap_fraction": 0.99,
            "min_reference_feature_overlap_targets": 100,
            "min_bed_feature_overlap_fraction": 0.95,
            "min_frozen_targets": 30,
        },
    )
    assert all(checks.values())


def test_build_target_catalog_freezes_only_four_way_supported_targets():
    spatial_ids = pd.Index(["cuTAR1", "cuTAR2", "cuTAR3"])
    counts = sparse.csr_matrix(
        [
            [1, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 0, 1],
        ],
        dtype=np.float32,
    )
    reference_var = pd.DataFrame(
        {"feature_type": ["cuTAR", "cuTAR", "cuTAR"]},
        index=spatial_ids,
    )
    feature_qc = pd.DataFrame(
        {
            "total_counts": [100, 100, 5],
            "detected_cells": [20, 20, 3],
            "quantified_cells": [1000, 1000, 1000],
            "quantified_sample_count": [3, 3, 3],
            "selection_tier": ["extended", "extended", "all_detected"],
        },
        index=spatial_ids,
    )
    bed = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2", "chr3"],
            "start": [10, 20, 30],
            "end": [15, 25, 35],
            "strand": ["+", "-", "+"],
        },
        index=spatial_ids,
    )
    catalog = build_target_catalog(
        spatial_ids,
        counts,
        reference_var,
        feature_qc,
        bed,
        min_reference_total_counts=20,
        min_reference_detected_cells=10,
        min_reference_detected_fraction=0.002,
        min_reference_samples=2,
        min_spatial_detected_spots=3,
    )
    assert bool(catalog.loc["cuTAR1", "frozen_primary"])
    assert not bool(catalog.loc["cuTAR2", "frozen_primary"])
    assert not bool(catalog.loc["cuTAR3", "frozen_primary"])
    assert (
        catalog.loc["cuTAR2", "target_status"]
        == "reference_supported_insufficient_spatial_truth"
    )
    assert (
        catalog.loc["cuTAR3", "target_status"]
        == "spatial_observed_insufficient_reference"
    )


def test_reference_detected_fraction_sets_dynamic_floor():
    spatial_ids = pd.Index(["cuTAR1"])
    counts = sparse.csr_matrix([[1], [1], [1]], dtype=np.float32)
    reference_var = pd.DataFrame(
        {"feature_type": ["cuTAR"]},
        index=spatial_ids,
    )
    feature_qc = pd.DataFrame(
        {
            "total_counts": [100],
            "detected_cells": [15],
            "quantified_cells": [10000],
            "quantified_sample_count": [3],
            "selection_tier": ["extended"],
        },
        index=spatial_ids,
    )
    bed = pd.DataFrame(
        {"chrom": ["chr1"], "start": [10], "end": [20], "strand": ["+"]},
        index=spatial_ids,
    )
    catalog = build_target_catalog(
        spatial_ids,
        counts,
        reference_var,
        feature_qc,
        bed,
        min_reference_total_counts=20,
        min_reference_detected_cells=10,
        min_reference_detected_fraction=0.002,
        min_reference_samples=2,
        min_spatial_detected_spots=3,
    )
    assert catalog.loc["cuTAR1", "required_reference_detected_cells"] == 20
    assert not bool(catalog.loc["cuTAR1", "reference_supported"])
