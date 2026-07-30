# LncSpaceMap

LncSpaceMap is a Python framework for mapping low-abundance long non-coding RNA
(lncRNA) signals from single-cell or single-nucleus RNA sequencing references
into spatial transcriptomics data.

The project follows a conservative rule: every spatial prediction must be
evaluated by fully masking observed targets, calibrated against genes with
comparable detection properties, accompanied by uncertainty, and allowed to
abstain when the reference evidence is insufficient.

## Project status

SPanC-Lnc reference preparation is complete. The finalized development
reference contains 56,557 matched cells, 33,538 genes, and 43,382 cuTAR
targets with explicit sample-quantification masks and genomic coordinates.
Week 1 contract, GENCODE annotation, target-eligibility, and leakage-free
masked-gene split implementation is active.

## Planned workflow

1. Harmonize gene identifiers and lncRNA annotations.
2. Determine whether each target lncRNA has sufficient reference evidence.
3. Construct within-sample and within-cell-type metacells.
4. Predict spatial expression with pluggable Tangram and SpaGE backends.
5. Calibrate or combine predictions using fully masked low-expression genes.
6. Correct predictions with a metacell-derived lncRNA-mRNA graph.
7. Apply edge-aware spatial refinement without crossing tissue boundaries.
8. Estimate uncertainty and assign high, medium, low, or abstain confidence.
9. Benchmark with leakage-free masking and external cross-platform validation.

## Repository map

- `docs/`: scientific design, architecture, data contract, benchmark protocol,
  compute profile, and dataset assessments.
- `configs/`: versioned configuration templates.
- `src/lncspacemap/`: future Python package modules.
- `tests/`: future unit, integration, and regression tests.
- `git_eval/`: lightweight results copied from external compute servers for
  review.

Raw sequencing data, large matrices, genome annotations, model checkpoints,
and complete result directories must not be committed. See `.gitignore` and
`git_eval/README.md`. Review artifacts must remain below 50 MB per file and
200 MB per reviewed run.

## Primary deliverables

- A spot- or cell-by-lncRNA AnnData prediction object.
- Relative-expression and presence-probability layers.
- Uncertainty intervals and per-target confidence metadata.
- Leakage-free benchmark tables covering six primary metrics.
- A reproducible run manifest and lightweight review report.

## Documentation

- [`docs/PROJECT_DESIGN.md`](docs/PROJECT_DESIGN.md)
- [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md)
- [`docs/DATA_CONTRACT.md`](docs/DATA_CONTRACT.md)
- [`docs/BENCHMARK_PROTOCOL.md`](docs/BENCHMARK_PROTOCOL.md)
- [`docs/COMPUTE_PROFILE.md`](docs/COMPUTE_PROFILE.md)
- [`docs/SPANCLNC_DATA_ASSESSMENT.md`](docs/SPANCLNC_DATA_ASSESSMENT.md)
- [`docs/SPANCLNC_PREPROCESSING.md`](docs/SPANCLNC_PREPROCESSING.md)
- [`docs/WEEK1_EXECUTION.md`](docs/WEEK1_EXECUTION.md)

## Resolution statement

LncSpaceMap preserves the resolution of the spatial input. A Visium input
produces spot-level predictions. Single-cell predictions are only claimed when
the spatial assay or an accepted segmentation provides genuine single-cell
coordinates.
