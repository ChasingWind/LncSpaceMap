# Week 1 baseline: MelD masked-gene mapping

This stage answers one question quickly: can a method reconstruct spatial
expression for low-expression protein-coding genes that were completely hidden
from model fitting? It is a software and leakage-control milestone, not yet the
final lncRNA benchmark.

## Required gates

1. The MelD 10x matrix contains non-negative integer raw counts.
2. Every retained barcode has one aligned tissue coordinate.
3. Reference and spatial genes use version-free Ensembl gene IDs.
4. At least 100 expressed shared protein-coding anchors remain.
5. Every evaluated target is absent from the anchor and fitting feature list.
6. Predictions preserve the exact spatial spot and target order.
7. Predictions contain no NaN, infinity, or negative values.
8. Each backend reports per-target Spearman correlation and z-NRMSE.

Both prediction and truth are evaluated on library-size-normalized
`log1p(counts per 10,000)` expression. Raw counts remain the source used for
target detection and input-contract validation.

Targets from the requested fold must be observed in at least three MelD spots.
If a reference-only fold has fewer than five evaluable genes, the pipeline
deterministically builds a MelD-specific panel of five folds by 25 genes from
shared protein-coding genes. Candidates must be detected in at least three but
no more than 25% of spots and are ranked for similarity to the original
reference low-expression proxies. Spatial truth is used only to establish
evaluation eligibility, never as model input for a masked target. The complete
fallback panel is exported for review and is shared unchanged by SpaGE and
Tangram. Tight distribution matching and the full metric panel remain deferred
until the final benchmark is frozen.

## Outputs

Large prediction H5AD files are written only to the server-side output
directory. `git_eval` receives the spatial contract, anchor/target manifests,
per-target metrics, summaries, package versions, and logs.

SpaGE uses its official normalized cell-by-gene DataFrame API. Tangram maps
reference cells to spots using anchors only and projects the withheld genes
after fitting. A 24 GB RTX 3090 is sufficient for the 56,557-cell by 623-spot
MelD mapping; A6000 is optional rather than required.
