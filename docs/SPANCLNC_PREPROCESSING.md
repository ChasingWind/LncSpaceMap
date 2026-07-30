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
