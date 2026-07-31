# Week 1D: bounded Tangram tuning

This stage permits one selection round and, only when warranted, one
confirmation round. It is intentionally small because the Week 1 proxy truth
is extremely sparse and is not the final lncRNA benchmark.

## Why the baseline scores are low

The median target is detected in only five of 617 spatial spots. The reference
and MelD spatial data also come from different specimens and technologies, so
single-cell and spatial dropout patterns differ. Tangram's documentation notes
that spatial sparsity and cross-dataset sparsity mismatch reduce mapping
quality. The completed run is therefore technically valid, while its
continuous-expression accuracy remains modest.

## Round 1: fold 0 selection

The frozen baseline is compared with three low-risk candidates:

| Candidate | Anchors | Epochs | Density prior |
|---|---:|---:|---|
| `anchors_1000` | 1,000 | 300 | RNA-count based |
| `anchors_4000` | 4,000 | 300 | RNA-count based |
| `uniform_2000` | 2,000 | 300 | Uniform |

The target panel, normalization, model mode, seed, and six metrics are
unchanged. A candidate must improve at least four metric directions, clear the
predefined material threshold for at least two metrics, and avoid a material
regression in both Spearman and detection AUROC.

## Round 2: folds 1-4 confirmation

Only the best eligible fold-0 candidate is run on folds 1-4. The same gate is
applied to pooled untouched-fold results. If fold 0 has no eligible candidate,
or the selected candidate does not reproduce its gain, the stage returns
`STOP_TUNING_ENTER_W2`. No third tuning round is allowed.

An accepted candidate returns `ACCEPT_TUNED_TANGRAM_ENTER_W2`. Either decision
advances to W2; the only difference is which Tangram configuration is frozen.
