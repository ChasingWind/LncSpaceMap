# LncSpaceMap Project Design

## 1. Scientific objective

LncSpaceMap predicts the spatial distribution of low-abundance lncRNAs that
are observed in a dissociated single-cell or single-nucleus reference but are
unmeasured, weakly measured, or deliberately masked in spatial transcriptomics
(ST).

The primary product is a calibrated relative spatial expression map, not a
claim of reconstructed molecule counts. The method must preserve the spatial
resolution of the input assay and provide an explicit abstention state.

## 2. Scope

### In scope for the minimum viable product

- Human and mouse count matrices represented as AnnData.
- Sequencing-based spot ST and single-cell-resolution ST.
- Annotated lncRNAs and pre-quantified custom transcripts such as cuTARs.
- Tangram and SpaGE as interchangeable baseline backends.
- Low-expression-aware gene masking and calibration.
- Metacell stabilization, gene-network correction, spatial refinement,
  uncertainty, and six primary benchmark metrics.
- CPU execution for preprocessing and SpaGE; single-GPU execution for Tangram
  and optional SpaIM comparison.

### Out of scope for the minimum viable product

- De novo transcript discovery from FASTQ or BAM.
- Transcript-boundary reconstruction or isoform assignment.
- Fabrication of predictions for targets absent from the reference.
- Histology-only expression prediction.
- Automatic super-resolution of Visium spots into cells.
- A universal pretrained lncRNA foundation model.

Custom transcripts can be mapped only after a separate quantification workflow
has produced compatible reference and validation matrices.

## 3. Design principles

1. **No target leakage:** a held-out target is removed from every anchor,
   embedding, model-fitting, and calibration input.
2. **Matched difficulty:** calibration genes match target lncRNAs by abundance,
   detection rate, dispersion, specificity, and sample reproducibility.
3. **Relative output:** continuous predictions are relative expression scores;
   they are not reported as observed UMIs.
4. **Boundary awareness:** spatial smoothing is allowed only across physical
   neighbors with compatible measured transcriptomes.
5. **Uncertainty first:** every target receives stability, interval, and
   confidence information.
6. **Abstention:** insufficient reference evidence produces no confident map.
7. **Frozen external test:** at least one dataset is never used for parameter
   selection.
8. **Backend independence:** LncSpaceMap adds lncRNA-specific calibration and
   governance around prediction backends instead of hiding backend behavior.

## 4. End-to-end workflow

### Stage 0: run registration

- Load a versioned YAML configuration.
- Record input checksums or stable file fingerprints.
- Record package versions, device information, random seeds, and git commit.
- Assign a deterministic run identifier.

**Gate P0:** configuration parses, required files exist, and no output path
overlaps an input path.

### Stage 1: data audit and harmonization

- Resolve raw count sources without silently treating log-normalized data as
  counts.
- Enforce unique observation and feature identifiers.
- normalize Ensembl IDs by an explicit policy.
- Join feature metadata to GENCODE/Ensembl or a user-provided target catalog.
- Validate spatial coordinates and sample labels.
- Produce reference, spatial, and target audit tables.

**Gate P1:** the datasets share enough reliable protein-coding anchors, the
target catalog is traceable, and all target IDs have an unambiguous definition.

### Stage 2: target eligibility

Default eligibility evidence for a target:

- total reference count at least 20;
- detected in at least `max(10 cells, 0.2% of cells)`;
- detectable in at least two biological samples when multiple samples exist,
  unless the target is explicitly declared sample-specific;
- not dominated by one low-quality cell or one technical batch;
- coordinate and strand information available for custom targets.

Targets are assigned `eligible`, `exploratory`, or
`insufficient_reference_signal`.

**Gate P2:** only eligible or explicitly requested exploratory targets proceed.

### Stage 3: metacell reference

- Split the reference by biological sample and cell type.
- Build local k-nearest-neighbor groups without mixing donors or broad cell
  types.
- Aggregate approximately 20-50 cells per metacell.
- retain the mapping from metacells to source cells.
- Compute target abundance, detection, dispersion, and specificity features.

**Gate P3:** metacells improve target detection without collapsing required
cell states or mixing biological samples.

### Stage 4: anchor construction

- Start from protein-coding genes measured in both reference and ST.
- Exclude mitochondrial, ribosomal, ambiguous, and target-derived features.
- Balance variable genes with cell-type and tissue-region markers.
- Use a default of approximately 2,000 anchors, subject to data availability.

**Gate P4:** anchor-only cross-dataset correspondence is adequate. A failed
alignment produces dataset-level abstention rather than forced prediction.

### Stage 5: baseline prediction

- Tangram maps reference cells or metacells onto the spatial support.
- SpaGE predicts held-out genes from a joint reference-spatial representation.
- Backends produce the same canonical `PredictionResult` contract.
- Optional backends can be evaluated without changing downstream modules.

### Stage 6: low-expression calibration and ensemble

- Select shared calibration genes matched to target lncRNA difficulty.
- Fully mask calibration genes in five gene-wise folds.
- Predict all masked calibration genes with each backend.
- Learn constrained backend weights by difficulty bin, not per target.
- Calibrate relative expression and presence probability.

The base ensemble is

`Y_base = w_bin * Y_tangram + (1 - w_bin) * Y_spage`.

If an ensemble is not demonstrably better, the best validated single backend
is retained.

**Gate P5:** ensemble or selected backend improves the primary metrics without
trading them for excessive smoothing.

### Stage 7: gene-network correction

- Construct a sparse robust lncRNA-mRNA network from metacells.
- Estimate masked-gene prediction errors.
- Propagate only bounded, well-supported corrections from calibration genes to
  targets.
- cap correction magnitude to prevent low-expression inflation.

This stage adapts the error-propagation concept used by SPRITE but applies
lncRNA-specific calibration and eligibility controls.

### Stage 8: edge-aware spatial refinement

Construct a graph whose edges require both physical adjacency and similarity
of measured anchor-gene PCs. Refinement uses

`Y_final = (1 - lambda) * Y_corrected + lambda * S * Y_corrected`,

where `lambda` is selected on masked low-expression genes and is conservatively
bounded.

**Gate P6:** refinement improves structural metrics and hotspot recovery
without erasing domain boundaries.

### Stage 9: uncertainty and abstention

- Bootstrap reference cells, metacells, or anchor subsets.
- measure disagreement among backends and bootstrap predictions.
- estimate calibrated intervals from matched masked genes.
- assign target confidence: `high`, `medium`, `low`, or `abstain`.
- distinguish target-level abstention from dataset-level alignment failure.

### Stage 10: benchmark and report

- Run leakage-free masked-gene evaluation.
- Evaluate actual observed lncRNAs separately from matched low-expression
  protein-coding genes.
- Run an untouched external or cross-platform validation.
- Produce metrics, plots, run manifests, resource summaries, and an audit
  report.

## 5. Four-week delivery plan

### Week 1: contract and baseline closure

- Implement AnnData audit and feature harmonization.
- Implement target catalog and eligibility.
- Add Tangram and SpaGE adapters.
- Complete one masked-gene end-to-end run.

**Exit:** one command produces baseline predictions and a valid audit report.

### Week 2: lncRNA-specific method

- Add metacells and matched calibration genes.
- Add constrained backend selection or ensemble.
- Add gene-network correction and edge-aware refinement.
- Complete core ablations.

**Exit:** final candidate beats or matches the best baseline under frozen
low-expression evaluation.

### Week 3: confidence and external evidence

- Add bootstrap stability and prediction intervals.
- Add confidence and abstention rules.
- Evaluate two development/validation datasets and one frozen external set.
- Run SPanC-Lnc cross-platform validation.

**Exit:** confidence is calibrated and actual lncRNA results agree with the
low-expression proxy benchmark.

### Week 4: package and reproducibility

- Freeze CLI, configuration, output schemas, and package versions.
- Add tests, example data, reports, and documentation.
- Run the full regression suite on RTX 3090.
- Run one large-memory stress test on A6000 if available.

**Exit:** reproducible release candidate with no undocumented manual steps.

## 6. Acceptance criteria

- At least three dataset contexts: development, validation, and frozen external.
- Final method improves median gene-wise Spearman and domain enrichment
  concordance over the best single backend in at least two contexts.
- NRMSE and SSIM do not materially regress.
- Spatial refinement does not create excessive Moran's I or diffuse hotspots
  across measured boundaries.
- Nominal 90% intervals achieve acceptable empirical coverage.
- Low-evidence targets are rejected rather than assigned confident maps.
- Every reported result is linked to a configuration and run manifest.

No universal absolute correlation threshold is declared before observing
dataset difficulty. Relative improvement, confidence intervals, and
cross-platform replication are the primary evidence.
