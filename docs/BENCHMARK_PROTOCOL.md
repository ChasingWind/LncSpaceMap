# LncSpaceMap Benchmark Protocol

## 1. Benchmark questions

1. Does the method recover the spatial rank and shape of fully hidden genes?
2. Does it work for genes with lncRNA-like abundance and sparsity?
3. Does it recover actual observed lncRNAs, not only proxy coding genes?
4. Does refinement improve prediction without manufacturing spatial smoothness?
5. Are confidence scores and intervals calibrated?
6. Does performance transfer to an untouched dataset or platform?

## 2. Evaluation levels

### Level A: matched low-expression coding genes

Use measured coding genes matched to target lncRNAs by:

- mean reference expression;
- reference detection rate;
- dispersion;
- cell-type specificity;
- sample reproducibility.

This level supplies enough genes for parameter selection and sensitivity
analysis but is not sufficient by itself.

### Level B: observed lncRNA masking

Select lncRNAs measured in both reference and ST. Remove each held-out fold from
all anchors and fitting inputs, then predict it. This is the main internal
lncRNA evidence.

### Level C: independent cross-platform validation

Freeze all parameters, predict targets using a reference and one spatial
support, and compare with another platform or sample not used for fitting.
Xenium, Curio/Seeker, and STOmics data are candidate validation supports.

## 3. No-leakage gene cross-validation

- Use five gene-wise folds.
- Stratify folds by feature class and expression difficulty.
- Remove the held-out fold from anchors, embeddings, backend fitting, gene
  graphs, and parameter selection.
- Do not perform random spot masking as a substitute for unmeasured-gene
  prediction.
- Keep a final external dataset frozen until the method and thresholds are
  selected.

## 4. Baselines and ablations

Minimum comparison:

1. Tangram.
2. SpaGE.
3. selected or weighted dual-backend prediction.
4. plus metacell and low-expression calibration.
5. plus gene-network correction.
6. plus edge-aware spatial refinement.
7. final method plus uncertainty and abstention.

Optional external comparison:

- SpaIM;
- SPRITE post-processing;
- gimVI when its data contract is suitable.

## 5. Six primary metrics

### Gene-wise Spearman correlation

Correlation between predicted and observed expression across spatial
observations. This is the primary rank/distribution metric.

### z-score NRMSE

RMSE after applying the same per-gene z-score transformation to prediction and
truth. Lower is better.

### Structural Similarity Index

SSIM compares overall expression structure. Coordinates must be converted to a
common raster or neighborhood representation by one frozen policy.

### Jensen-Shannon divergence

Prediction and truth are converted to spatial probability distributions.
Lower divergence is better.

### Domain Enrichment Concordance

For domain-wise expression proportions `p_gk`:

`DEC_g = 1 - 0.5 * sum_k(abs(p_hat_gk - p_gk))`.

The score ranges from 0 to 1. Domains are external annotations or are frozen
from measured anchor genes before target masking.

### Low-Expression Hotspot F1

Apply one frozen hotspot definition to prediction and truth, then compute the
F1 score between hotspot sets. The hotspot fraction is informed by observed
detection and constrained to 5%-30%.

## 6. Secondary diagnostics

- absolute Moran's I error;
- predicted-to-observed mean and variance ratios;
- zero/detection calibration;
- bootstrap map stability;
- backend disagreement;
- empirical prediction-interval coverage;
- runtime, peak host RAM, and peak GPU memory.

Metrics that directly reward smoothness are never interpreted without Moran's
I error and hotspot recovery.

## 7. Statistical reporting

- Report gene-level distributions, not only pooled means.
- Report median and bootstrap 95% confidence interval.
- Use paired gene-level comparisons for ablations.
- Separate datasets, actual lncRNAs, and proxy genes.
- Report performance by abundance, detection, specificity, and confidence bin.
- Include failure and abstention rates.

## 8. Model-selection rule

The selected model must improve median Spearman and DEC over the best baseline
and must not materially worsen NRMSE or SSIM. If network correction or spatial
refinement fails this rule, the component is disabled.

No component is retained solely because it improves a composite rank.

## 9. External-test rule

After parameter freeze:

- no threshold or weight changes based on external results;
- report every eligible target, not only favorable examples;
- define primary targets before inspecting prediction maps;
- treat different tissue conditions as stress tests, not matched validation;
- use cross-platform concordance as evidence of localization, not proof of
  exact transcript counts.
