# SPanC-Lnc preprocessing

The preprocessing stage never modifies the downloaded files.

## Stage 1: audit

```bash
cd ~/Spatial/LncSpaceMap
conda activate gstnet
python scripts/process_spanc_lnc.py audit \
  --data-dir ~/Spatial/data/Spanc-Lnc
```

Review:

- `git_eval/metrics/spanc_lnc_matrix_audit.tsv`
- `git_eval/metrics/spanc_lnc_cutar_audit.tsv`
- `git_eval/manifests/spanc_lnc_file_manifest.tsv`
- `git_eval/logs/spanc_lnc_preprocess.log`

The enhanced audit verifies all of the following before build:

1. the UQ matrices use a supported orientation and header layout;
2. count previews are finite, non-negative integers;
3. cuTAR identifiers are unique and at least 95% occur in `00_cuTARs.bed`;
4. 10x H5 and cuTAR barcodes match unambiguously;
5. barcode overlap is at least 80% on both sides;
6. every matrix and spatial-coordinate input check passes.

Do not run the build stage unless the log ends with
`PASS_SPANCLNC_ENHANCED_AUDIT_READY_FOR_BUILD`. A
`HOLD_BUILD_ENHANCED_AUDIT_FAILED` marker means the failing audit rows must be
reviewed first.

The UQ `Acral_Mel*_cuTAR_mat.txt` files use a valid but unusual header: the
first line contains only cell barcodes, while every data line contains one
additional leading cuTAR identifier. The reader handles this explicitly and
does not rely on pandas' implicit-index inference.

## Stage 2: build

```bash
python scripts/process_spanc_lnc.py build \
  --data-dir ~/Spatial/data/Spanc-Lnc \
  --output-dir ~/Spatial/data/Spanc-Lnc/processed/reference
```

The build:

1. reads each dense cuTAR text matrix in small row chunks;
2. converts it to a sparse cells-by-cuTAR matrix;
3. validates finite, non-negative integer counts;
4. matches cells by the least destructive unambiguous barcode strategy;
5. requires at least 80% overlap on both sides;
6. combines genes and cuTARs while preserving raw counts;
7. writes six sample H5AD files and one combined H5AD.

AM4 and the four cutaneous melanoma samples are intentionally excluded because
the downloaded UQ manifest does not identify matching cuTAR matrices for them.

## Expected large outputs

Large outputs remain outside git:

```text
~/Spatial/data/Spanc-Lnc/processed/reference/
├── AM1.gene_cutar.h5ad
├── AM2.gene_cutar.h5ad
├── AM3_pre.gene_cutar.h5ad
├── AM3_post.gene_cutar.h5ad
├── AM5.gene_cutar.h5ad
├── AM6.gene_cutar.h5ad
└── acral_melanoma_gene_cutar_combined.h5ad
```

Only the lightweight audit tables and logs under `git_eval` should be
committed for review.

## Stage 3: finalize and reference QC

Use this mode after the six sample H5AD files have already been built. It does
not reread the large raw cuTAR text matrices.

```bash
python scripts/process_spanc_lnc.py finalize \
  --data-dir ~/Spatial/data/Spanc-Lnc \
  --output-dir ~/Spatial/data/Spanc-Lnc/processed/reference
```

Finalize performs the following:

1. reconstructs `feature_type` and available gene metadata in the combined
   H5AD;
2. attaches cuTAR chromosome, start, end, and strand from `00_cuTARs.bed`;
3. records one `quantified_<sample>` mask per sample plus
   `quantified_sample_count`, preventing structural absence from being treated
   as a biological zero;
4. calculates cell-level gene/cuTAR counts and detection rates;
5. calculates feature-level counts, prevalence, and sample support;
6. assigns provisional `all_detected`, `extended`, and `core` cuTAR tiers;
7. atomically replaces the earlier combined H5AD only after validation and
   writes `PASS_SPANCLNC_REFERENCE_QC_READY_FOR_ANNOTATION` to the log.

Detailed QC remains outside git:

```text
processed/reference/qc/
├── spanc_lnc_cell_qc.tsv.gz
└── spanc_lnc_feature_qc.tsv.gz
```

Lightweight review outputs are written to:

- `git_eval/metrics/spanc_lnc_combined_contract.tsv`
- `git_eval/metrics/spanc_lnc_reference_cell_qc_summary.tsv`
- `git_eval/metrics/spanc_lnc_cutar_tier_summary.tsv`

The tier rules are deliberately provisional. `extended` requires detection in
at least 10 cells and quantification in at least two samples; `core` requires
detection in at least 30 cells and quantification in at least three samples.
Both require at least 0.1% prevalence among cells in which the feature was
quantified. No feature is removed at this stage.
