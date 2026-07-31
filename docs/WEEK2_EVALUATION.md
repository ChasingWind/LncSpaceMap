# Week 2C: real cuTAR evaluation

Week 2C compares the frozen Week 2B prediction object with the independently
stored Week 2A MelD cuTAR truth. No mapping model is refit in this stage.

Continuous truth is normalized as
`log1p(cuTAR count / (gene counts + all cuTAR counts) * 10000)`. Detection
truth remains the raw condition `cuTAR count > 0`.

The six primary metrics are Pearson, Spearman, z-NRMSE, detection AUROC,
detection AUPRC, and top-k spatial recall. AUPRC lift, top-k lift, 200
permutation nulls, prediction coverage, and truth-positive coverage are
reported as calibration metadata.

Abstained spot-target pairs remain missing. Five metrics use valid spots only;
top-k recall retains all observed-positive spots in its denominator so that an
abstention on a true-positive spot cannot improve performance.

Summaries are frozen for all 362 targets, the 163 targets detected in at least
six spots, the 51 targets detected in more than ten spots, the 199 targets
detected in three to five spots, and reference quantification in two versus at
least three samples.

The primary `>=6 spots` panel passes the bounded scientific gate when at least
three of four criteria pass: median Spearman >=0.05, AUROC >=0.55, AUPRC lift
>=1.25, and top-k lift >=1.25. A failed scientific gate does not invalidate the
run; it routes the project to one bounded comparator rather than calibration.
