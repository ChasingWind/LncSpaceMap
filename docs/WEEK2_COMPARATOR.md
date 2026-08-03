# Week 2C2: bounded Tangram projection comparator

This is the only post-W2C Tangram comparison. It does not refit the mapping or
change anchors, epochs, density prior, target eligibility, truth, metrics, or
permutation count. It evaluates the existing `raw_projection` layer and
compares it target-by-target with the frozen `relative_expression` results.

The primary panel remains the 163 targets detected in at least six MelD spots.
Paired directional differences and one-sided Wilcoxon tests are reported for
all six metrics, AUPRC lift, and top-k lift. z-NRMSE improvements are oriented
so that positive values always mean raw projection is better.

Raw projection is admitted only if it passes the original W2C scientific gate:
at least three of median Spearman >=0.05, AUROC >=0.55, AUPRC lift >=1.25, and
top-k lift >=1.25. A merely relative improvement is insufficient.

The terminal decisions are:

- `ADMIT_W2D_CALIBRATION_RAW_PROJECTION`; or
- `FREEZE_TANGRAM_NEGATIVE_BASELINE_ENTER_W3_RELATION_MODEL`.

No further Tangram parameter tuning is permitted after this comparison.
