# Week 2A: observed MelD cuTAR contract

Week 2A is a data-admission stage. It does not fit Tangram or use observed
spatial cuTAR counts as model input.

## Inputs

- finalized acral-melanoma gene-plus-cuTAR single-cell reference;
- finalized reference feature-QC table;
- aligned MelD gene spatial H5AD from Week 1;
- author-provided `MelD.txt` cuTAR count matrix;
- `00_cuTARs.bed`;
- `configs/week2_cutar.yaml`.

The text matrix is parsed incrementally. Orientation, delimiter, header style,
raw integer counts, duplicate identifiers, and barcode normalization are
resolved before any target is admitted.

## Four-way target contract

A primary target must simultaneously have:

1. an exact cuTAR ID in the observed MelD matrix;
2. the same exact cuTAR ID in the single-cell reference;
3. a valid BED coordinate and strand;
4. sufficient reference and observed spatial evidence.

Reference evidence requires at least 20 total counts, detection in at least
`max(10 cells, 0.2% of quantified cells)`, and quantification in at least two
reference samples. Spatial evaluability requires detection in at least three
matched MelD spots.

No target aliases, coordinate-overlap substitutions, or arbitrary duplicate
resolution are allowed in this stage.

## Outputs

Large server-only outputs:

- `meld_cutar_truth_aligned.h5ad`;
- `w2a_meld_cutar_target_catalog.tsv.gz`.

Lightweight review outputs:

- `git_eval/metrics/w2a_meld_cutar_contract.tsv`;
- `git_eval/metrics/w2a_meld_target_overlap_summary.tsv`;
- `git_eval/metrics/w2a_meld_target_tier_summary.tsv`;
- `git_eval/manifests/w2a_meld_frozen_targets.tsv`;
- `git_eval/manifests/w2a_meld_cutar_manifest.json`;
- `git_eval/logs/week2a_meld_cutar.log`.

`ADMIT_W2B_REAL_CUTAR_MAPPING` means the aligned truth object and frozen target
panel can be used for leakage-free projection and evaluation. Spatial cuTAR
counts remain truth-only throughout W2-B.
