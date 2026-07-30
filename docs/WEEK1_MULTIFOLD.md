# Week 1C: five-fold baseline evaluation

This stage reruns all five MelD spatially evaluable proxy folds with SpaGE and
Tangram. Each backend is evaluated on 125 fully withheld protein-coding genes
using:

1. Spearman correlation.
2. Pearson correlation.
3. z-normalized RMSE.
4. Detection AUROC.
5. Detection AUPRC.
6. Top-k spatial recovery, where k equals the number of truth-positive spots.

Every per-target metric is calibrated against 200 spatial permutations of the
prediction vector. Results are also stratified into 3-5, 6-10, and more than
10 truth-detected spots.

The stage chooses only a provisional backbone. It does not declare production
lncRNA mapping accuracy. A backend wins a metric using the five-fold median
and becomes the provisional winner only if it wins at least four of six
metrics. Leakage, spot order, and target-panel identity must pass before a
winner is reported.
