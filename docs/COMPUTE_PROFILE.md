# LncSpaceMap Compute Profile

## 1. Supported primary profile

- GPU: one NVIDIA RTX 3090 with 24 GB VRAM.
- Host RAM: 64 GB minimum, 128 GB recommended.
- CPU: 16-32 threads.
- Storage: at least 200 GB fast local workspace for development; more when raw
  sequencing data are retained outside the repository.

RTX 3090 is the release target. A6000 must not become a runtime requirement.

## 2. A6000 profile

An RTX A6000 with 48 GB VRAM is useful for:

- large cell-level Tangram mapping;
- many spatial observations;
- larger SpaIM comparisons;
- stress tests without aggressive metacell reduction;
- memory profiling and upper-bound experiments.

It is a capacity resource, not a scientific requirement.

## 3. Default resource controls

- metacell mapping by default;
- no more than approximately 10,000 metacells in the initial release profile;
- approximately 2,000 anchor genes;
- target prediction in chunks;
- one CV fold per GPU at a time;
- sparse matrices and backed H5AD access;
- float32 for the first validated release;
- deterministic seeds where supported.

Mixed precision is optional only after demonstrating that low-expression
targets are numerically stable.

## 4. Tangram memory planning

The learned mapping scales approximately with
`n_reference_units * n_spatial_units`. A rough planning lower bound for a
float32 parameter with gradients and Adam state is:

`16 * n_reference_units * n_spatial_units bytes`,

before temporary tensors and framework overhead.

Practical guidance:

- 10,000 metacells x 5,000 spots: comfortable on RTX 3090.
- 50,000 cells x 5,000 spots: generally feasible but requires profiling.
- 100,000 cells x 10,000 spots: high risk on 24 GB; use metacells or A6000.

## 5. CPU and RAM-heavy components

The following can be limited by host RAM rather than GPU:

- H5AD conversion and accidental dense copies;
- Curio uTAR matrix inspection;
- robust gene-correlation construction;
- cross-validation metric aggregation;
- bootstrap summaries;
- report generation from large gene-level tables.

The 11 GB Curio breast uTAR H5AD should be read in backed mode and subset before
materialization.

## 6. Profiling requirements

Every production run records:

- GPU model and VRAM;
- CUDA, driver, and PyTorch versions;
- peak allocated and reserved GPU memory;
- peak host resident memory;
- stage wall time;
- number of reference units, spots/cells, anchors, and targets;
- fallback or OOM recovery actions.

Resource results are copied into `git_eval/metrics/resource_usage.tsv`.
