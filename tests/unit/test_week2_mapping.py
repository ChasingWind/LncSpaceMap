import numpy as np
import pandas as pd
from scipy import sparse

from lncspacemap.pipeline.week2_mapping import (
    build_quantification_mask,
    project_mask_aware,
)


class ReferenceStub:
    def __init__(self):
        self.n_obs = 3
        self.obs = pd.DataFrame(
            {"sample_batch": ["AM1", "AM1", "AM2"]},
            index=["c1", "c2", "c3"],
        )
        self.var = pd.DataFrame(
            {
                "quantified_AM1": [True, False],
                "quantified_AM2": [False, True],
            },
            index=["cuTAR1", "cuTAR2"],
        )


def test_build_quantification_mask_uses_cell_sample():
    mask = build_quantification_mask(ReferenceStub(), ["cuTAR1", "cuTAR2"])
    np.testing.assert_array_equal(mask, [[1, 0], [1, 0], [0, 1]])


def test_mask_aware_projection_excludes_structural_zero_cells():
    mapping = np.asarray([[1.0], [1.0], [1.0]], dtype=np.float32)
    expression = sparse.csr_matrix([[2.0], [4.0], [0.0]], dtype=np.float32)
    quantified = np.asarray([[1.0], [1.0], [0.0]], dtype=np.float32)
    raw, relative, support = project_mask_aware(
        mapping, expression, quantified, min_support_fraction=0.5
    )
    np.testing.assert_allclose(raw, [[6.0]])
    np.testing.assert_allclose(relative, [[3.0]])
    np.testing.assert_allclose(support, [[2.0 / 3.0]])


def test_low_reference_support_abstains():
    mapping = np.asarray([[0.9], [0.1]], dtype=np.float32)
    expression = np.asarray([[0.0], [5.0]], dtype=np.float32)
    quantified = np.asarray([[0.0], [1.0]], dtype=np.float32)
    _, relative, support = project_mask_aware(
        mapping, expression, quantified, min_support_fraction=0.2
    )
    np.testing.assert_allclose(support, [[0.1]])
    assert np.isnan(relative[0, 0])
