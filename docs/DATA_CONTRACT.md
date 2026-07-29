# LncSpaceMap Data Contract

## 1. Required inputs

### Reference AnnData

Required:

- observations are cells, nuclei, or declared metacells;
- variables are genes/custom transcripts;
- a declared raw non-negative count source;
- unique observation and variable identifiers;
- `obs` sample identifier;
- `obs` cell type or cluster identifier, strongly recommended;
- feature ID and symbol metadata;
- species and genome build.

Accepted count locations are configured explicitly, for example:

- `layers["counts"]`;
- `raw.X`;
- `X`.

The program must not choose between these sources by silent heuristic.

### Spatial AnnData

Required:

- observations are spatial spots or segmented cells;
- variables are measured genes/custom transcripts;
- a declared raw count source;
- unique observation and variable identifiers;
- `obsm["spatial"]` or configured coordinate columns;
- sample identifier;
- assay and resolution metadata.

Recommended:

- tissue image metadata;
- region or domain labels;
- cell-type annotations for single-cell ST;
- batch or slide identifier.

### Target catalog

Targets can originate from:

- GENCODE/Ensembl lncRNA annotations;
- a user-provided feature table;
- a BED-derived custom catalog such as SPanC-Lnc cuTARs.

Minimum fields:

- `target_id`;
- `feature_class`;
- `annotation_source`;
- `genome_build`.

Custom genomic targets also require chromosome, start, end, and strand.

## 2. Identifier policy

- Ensembl version suffix removal is optional and recorded.
- Symbols are display fields, not primary join keys when stable IDs exist.
- Duplicate IDs are never resolved by keeping an arbitrary first occurrence.
- One-to-many aliases require explicit aggregation or exclusion.
- Custom targets must use the same coordinate system and genome build across
  reference and spatial data.

## 3. Expression policy

- Raw counts are preserved.
- Normalized views are generated as separate matrices or temporary objects.
- Baseline backends receive the transformation required by their interface.
- Evaluation compares predictions and truth after one documented common
  transformation.
- Predicted continuous values are called relative expression, not UMIs.

## 4. Anchor policy

Anchor genes must:

- exist in reference and ST;
- be reliably detected;
- be independent of target definitions;
- exclude every masked calibration target;
- exclude genes with ambiguous mapping;
- balance broad variability and biological marker coverage.

The exact anchor list is exported for every run.

## 5. Target eligibility defaults

Default target evidence:

- reference total counts >= 20;
- detected cells >= `max(10, 0.002 * n_reference_cells)`;
- sample support >= 2 when multiple reference samples exist;
- no single-cell or single-sample domination;
- valid annotation.

Thresholds are configurable and evaluated in sensitivity analyses.

## 6. Spatial coordinate checks

- coordinates must be finite;
- duplicate coordinates are allowed only if the assay semantics permit them;
- sample-specific coordinate systems must not be connected across slides;
- coordinate units and origin are recorded when available;
- spatial graphs are constructed within each sample.

## 7. Dataset compatibility gates

The pipeline rejects or abstains when:

- shared reliable coding anchors are insufficient;
- reference and ST represent incompatible species or genome builds;
- custom target IDs cannot be reconciled;
- raw versus normalized expression is unresolved;
- spatial coordinates are missing;
- the reference contains no evidence for the requested targets;
- cross-dataset anchor validation indicates a failed alignment.

## 8. SPanC-Lnc-specific checks

Before using supplied files:

- verify all H5AD matrices are sparse or can be accessed in backed mode;
- report `shape`, `X` type, layers, raw, obs, var, obsm, and uns keys;
- verify Curio gene and uTAR observation IDs overlap exactly or document loss;
- verify melanoma uTAR text orientation and barcode naming;
- split STOmics objects by sample/tissue before constructing graphs;
- verify MelD/MelDN uTAR count barcodes match tissue positions;
- confirm `00_cuTARs.bed` uses GRCh38/hg38 coordinates;
- retain source confidence flags such as database overlap and less-confident
  downstream-of-coding-gene annotations.

## 9. Lightweight audit outputs

The external server should copy only these summaries into `git_eval`:

- dataset shape and matrix source;
- obs/var schema;
- feature-class counts;
- barcode and gene-overlap tables;
- count-distribution summaries;
- coordinate diagnostics;
- target eligibility counts;
- errors and warnings.

Raw matrices and complete metadata containing sensitive identifiers stay on
the compute server.
