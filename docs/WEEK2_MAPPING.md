# Week 2B: real MelD cuTAR mapping

Week 2B fits one Tangram cell-to-spot map from 2,000 shared protein-coding
anchor genes and projects the frozen Week 2A cuTAR panel. Observed spatial
cuTAR counts are not opened or supplied to this stage.

## Leakage and missingness controls

- cuTAR targets are excluded from mapping features;
- the cell-to-spot map is trained once and reused for every target batch;
- reference zeros from samples that did not quantify a cuTAR are structural
  missing values, not biological zero expression;
- sample-specific `quantified_<sample>` metadata defines a cells-by-targets
  quantification mask;
- predictions below the configured mapped-reference support fraction are
  represented as missing and must be treated as abstentions.

For mapping matrix \(M\), reference target expression \(X\), and quantification
mask \(Q\), the adjusted prediction is:

\[
\hat X = \frac{M^T X}{M^T Q}.
\]

The output retains the unadjusted Tangram projection, adjusted relative
expression, and reference-support fraction as separate layers.

## Outputs

Server-only:

- `meld_w2b_cutar_predictions.h5ad`;
- `meld_w2b_cell_to_spot_map.h5ad`.

Lightweight review outputs:

- `git_eval/metrics/w2b_meld_mapping_contract.tsv`;
- `git_eval/metrics/w2b_meld_target_support.tsv`;
- `git_eval/manifests/w2b_meld_anchors.tsv`;
- `git_eval/manifests/w2b_meld_mapping_manifest.json`;
- `git_eval/logs/week2b_meld_mapping.log`.

`ADMIT_W2C_REAL_CUTAR_EVALUATION` permits the frozen predictions to be compared
with the separately stored aligned MelD cuTAR truth in Week 2C.
