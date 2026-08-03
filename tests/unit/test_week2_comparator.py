import numpy as np
import pandas as pd

from lncspacemap.pipeline.week2_comparator import paired_method_comparison


def _metrics(value):
    return pd.DataFrame(
        {
            "evaluation_panels": ["all_frozen;primary_ge6"] * 3,
            "spatial_detection_bin": ["6-10 spots"] * 3,
            "reference_quantified_samples": [3] * 3,
            "pearson": [value, value + 0.1, value + 0.2],
            "spearman": [value, value + 0.1, value + 0.2],
            "z_nrmse": [1.4 - value, 1.3 - value, 1.2 - value],
            "detection_auroc": [0.5 + value] * 3,
            "detection_auprc": [0.1 + value] * 3,
            "detection_auprc_lift": [1.0 + value] * 3,
            "topk_recall": [0.1 + value] * 3,
            "topk_recall_lift": [1.0 + value] * 3,
        },
        index=["a", "b", "c"],
    )


def test_paired_comparison_uses_metric_direction():
    relative = _metrics(0.0)
    raw = _metrics(0.1)
    paired, summary = paired_method_comparison(
        relative,
        raw,
        panel="primary_ge6",
        material_delta={
            "spearman": 0.05,
            "z_nrmse": 0.05,
            "detection_auroc": 0.05,
        },
    )
    np.testing.assert_allclose(paired["improvement_spearman"], 0.1)
    np.testing.assert_allclose(paired["improvement_z_nrmse"], 0.1)
    assert bool(
        summary.set_index("metric").loc["spearman", "material_improvement"]
    )
    assert bool(
        summary.set_index("metric").loc["z_nrmse", "material_improvement"]
    )
