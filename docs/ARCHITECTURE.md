# LncSpaceMap Software Architecture

## 1. Architectural style

LncSpaceMap is a stage-based Python package with explicit data contracts.
Stages are deterministic under a fixed configuration and seed. Heavy backends
are isolated behind adapters so that method-specific data conversions cannot
silently alter the canonical inputs.

## 2. Planned package tree

```text
src/lncspacemap/
├── cli/                 command-line entrypoints
├── contracts/           canonical data and result schemas
├── io/                  AnnData, annotation, and manifest I/O
├── preprocessing/       audit, harmonization, normalization, anchors
├── metacells/           reference stabilization
├── mapping/
│   └── backends/        Tangram, SpaGE, optional comparison adapters
├── calibration/         matched genes, masking, ensemble, calibration
├── refinement/          gene-network and spatial-graph refinement
├── uncertainty/         bootstrap, intervals, confidence, abstention
├── evaluation/          metrics, splits, ablations, statistical summaries
├── reporting/           tables, plots, HTML/Markdown reports
├── pipeline/            stage orchestration and checkpoint registry
└── utils/               seeds, sparse operations, validation helpers
```

## 3. Canonical objects

### ReferenceBundle

- reference AnnData or backed view;
- raw count source;
- sample and cell-type fields;
- canonical gene table;
- metacell membership and profiles;
- reference QC summary.

### SpatialBundle

- spatial AnnData or backed view;
- raw count source;
- coordinates and optional image metadata;
- sample and spatial-domain fields;
- canonical gene table;
- spatial-neighbor graph;
- spatial QC summary.

### TargetCatalog

One row per target with:

- canonical target ID and optional display symbol;
- source annotation and genome build;
- feature class and biotype;
- coordinate and strand for custom targets;
- detection, abundance, dispersion, and specificity;
- eligibility and exclusion reason;
- confidence annotations supplied by the source resource.

### MappingResult

- backend name and version;
- predicted relative-expression matrix;
- mapping diagnostics;
- backend-specific artifacts referenced by path, not embedded blindly;
- anchor list and no-leakage confirmation.

### CalibrationResult

- matched calibration-gene table;
- gene-fold assignments;
- backend weights or selected backend;
- presence calibration model;
- selected spatial regularization parameters;
- fold-level metrics.

### PredictionResult

- final relative-expression matrix;
- presence-probability matrix;
- lower and upper interval matrices;
- target confidence table;
- dataset and target abstention flags;
- provenance manifest.

## 4. Component responsibilities

### `io`

- Read AnnData in memory or backed mode.
- Detect and resolve count layers by declared policy.
- Read GTF-derived feature tables, BED target catalogs, and target lists.
- Write prediction AnnData, tabular summaries, and manifests.
- Never guess whether `X` contains counts or normalized values.

### `preprocessing`

- Validate indices, sparse formats, non-negativity, and finite values.
- Harmonize identifiers and record every dropped or merged feature.
- Select anchors without target leakage.
- Compute reusable normalized views while preserving raw counts.

### `metacells`

- Build sample- and cell-type-constrained local aggregates.
- Expose diagnostics for cell-state preservation and target recovery.
- Preserve source-cell membership for audit and bootstrapping.

### `mapping.backends`

Each adapter implements:

- `prepare(reference, spatial, anchors, config)`;
- `fit(...)`;
- `predict(targets)`;
- `diagnostics()`;
- `save_manifest()`.

Backend outputs are converted to the canonical spatial observation order and
target order before leaving the adapter.

### `calibration`

- Match calibration genes to target difficulty.
- Create leakage-free gene folds.
- Select a backend or learn constrained bin-level weights.
- Calibrate presence scores and prediction scale.

### `refinement`

- Build robust target-anchor gene networks.
- Apply bounded prediction-error correction.
- Build physical-plus-transcriptomic spatial graphs.
- Select and apply conservative graph refinement.

### `uncertainty`

- Run reproducible bootstraps.
- summarize backend and bootstrap disagreement.
- fit calibrated intervals from matched held-out genes.
- implement high/medium/low/abstain governance.

### `evaluation`

- Implement six primary metrics and secondary diagnostics.
- Aggregate gene-, fold-, dataset-, and expression-bin results.
- Compare ablations with paired bootstrap intervals.
- keep proxy low-expression genes separate from actual lncRNA results.

### `reporting`

- Generate machine-readable CSV/TSV/JSON outputs.
- Generate lightweight plots and a human-readable report.
- Redact raw paths and sensitive sample identifiers from review artifacts.

## 5. Pipeline state and checkpointing

Each stage writes:

- `stage_manifest.json`;
- a compact audit table;
- references to large local artifacts;
- success, failure, or abstention status;
- configuration and input fingerprints.

Large stage outputs live outside git. Only selected lightweight summaries are
copied into `git_eval`.

## 6. Prediction AnnData contract

The final AnnData has observations equal to spatial spots/cells and variables
equal to target lncRNAs.

- `X`: final relative-expression prediction.
- `layers["presence_probability"]`: calibrated probability.
- `layers["prediction_lower"]`: lower interval.
- `layers["prediction_upper"]`: upper interval.
- `obs`: copied spatial identifiers and safe metadata.
- `obsm["spatial"]`: original coordinates.
- `var`: target identity, eligibility, confidence, stability, and provenance.
- `uns["lncspacemap"]`: schema version, configuration hash, backend versions,
  metrics summary, and resolution statement.

This object is a large local artifact and is ignored by git.

## 7. Failure model

Failures are typed:

- `DATA_CONTRACT_ERROR`;
- `ANNOTATION_MISMATCH`;
- `INSUFFICIENT_ANCHORS`;
- `DATASET_NOT_ALIGNABLE`;
- `BACKEND_FAILURE`;
- `CALIBRATION_FAILURE`;
- `INSUFFICIENT_REFERENCE_SIGNAL`;
- `UNCERTAINTY_NOT_CALIBRATED`.

Scientific abstention is not treated as a software crash.

## 8. Reproducibility

- One root seed expanded deterministically by stage and fold.
- Stable observation and feature ordering.
- Configuration and dependency lock files.
- GPU, CUDA, PyTorch, backend, and package versions in manifests.
- No implicit downloads during a production run.
- No absolute server paths in committed reports.
