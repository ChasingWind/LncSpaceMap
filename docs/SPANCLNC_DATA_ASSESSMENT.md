# SPanC-Lnc Dataset Assessment for LncSpaceMap

## 1. Overall decision

Keep all downloaded files. They are appropriate and valuable, but they are
primarily spatial ground-truth and cross-platform validation resources. As
downloaded, they do not yet form a complete LncSpaceMap development input
because a compatible scRNA reference containing the same annotated/custom
lncRNA targets is still required.

The SPanC-Lnc study used Visium plus scRNA-seq for discovery and used
Curio/Seeker and STOmics mainly for validation. This distinction should be
preserved in LncSpaceMap.

## 2. File-by-file role

| File | Keep | Proposed role | Main condition |
|---|---:|---|---|
| `00_cuTARs.bed` | yes | cuTAR target catalog and coordinate bridge | confirm GRCh38/hg38, strand, stable cuTAR IDs, and coordinate convention |
| `Curio_breast_cancer_genes.h5ad` | yes | measured coding-gene spatial support at single-cell resolution | coordinates, counts, barcodes, and annotations must be present |
| `Curio_breast_cancer_uTARs.h5ad` | yes | high-value observed lncRNA/cuTAR ground truth | merge exactly to the gene object by observation ID; use backed mode |
| `Curio_melanoma_genes.h5ad` | yes | measured coding-gene spatial support | coordinates and observation IDs must align with the uTAR table |
| `Curio_melanoma_uTARs_counts.txt.gz` | yes | melanoma cuTAR ground truth | verify orientation and merge by barcode; an H5AD replacement is not required |
| `MelD/` | yes | melanoma Visium target/support and internal masking data | require tissue positions, gene matrix, image metadata, and matching cuTAR counts |
| `MelDN/` | yes | dysplastic-nevus Visium stress-test/support | it is a different biological state, not an interchangeable replicate of MelD |
| `STomics_gene_cuTARs_bin50.h5ad` | yes | independent total-RNA, high-resolution external validation | split melanoma and CRC samples; validate bin coordinates and count source |

## 3. Specific cautions

### Curio breast cancer

The gene and uTAR objects are complementary parts of one spatial assay. Their
value comes from joining them by the same spatial observations. The 11 GB uTAR
object must not be loaded densely. First read metadata in backed mode, quantify
barcode overlap, then materialize only eligible targets.

This is an excellent Level B/Level C ground-truth dataset after a breast cancer
scRNA reference has been added.

### Curio melanoma

The absence of a melanoma uTAR H5AD is not itself a reason to replace the data.
The compressed count table can be attached to the gene H5AD if:

- rows/columns and delimiter are known;
- barcode naming can be reconciled without ambiguous truncation;
- every retained barcode has coordinates from the gene object;
- cuTAR IDs map to `00_cuTARs.bed`;
- counts are raw, non-negative, and sparse after conversion.

### MelD and MelDN

These Visium samples are central to the melanoma use case described in the
paper. They are suitable when the folder contains the standard spatial
positions and gene matrix plus a spot-by-cuTAR count table or an equivalent
custom-feature matrix.

MelDN is dysplastic nevus and should be labeled as a biological stress test.
Do not merge MelD and MelDN as technical replicates and do not connect their
spatial graphs.

If the folders contain only a standard gene matrix, they are not sufficient for
observed-cuTAR validation. In that case, retain them as spatial supports and
obtain or rebuild the cuTAR count matrices from the authors' outputs/BAMs.

### STOmics

The study used STOmics random-hexamer total-RNA data to complement poly(A)
capture and reported it as validation rather than the main discovery set.
This makes the object particularly useful for testing transfer across capture
chemistry, but also creates a real domain shift.

Use it as frozen external validation after splitting samples and verifying that
gene and cuTAR features share one count scale. Do not use it to tune every
hyperparameter and then report it as external evidence.

### cuTAR annotation limitations

SPanC-Lnc constructed uTARs using coverage bins and merging nearby active
regions. The paper notes that neighboring transcripts can be merged and that
higher-resolution technologies can produce different boundaries. Therefore:

- map by stable cuTAR ID and documented overlap, not symbol alone;
- retain annotation confidence;
- flag targets close to coding genes or with ambiguous boundaries;
- evaluate annotated GENCODE lncRNAs and novel cuTARs separately.

## 4. Required additions

### Priority 1: compatible scRNA references

For melanoma:

- obtain the public acral/cutaneous melanoma scRNA data used by the study
  (`PRJNA862451`) or the 10x primary melanoma dataset;
- obtain raw or author-processed cuTAR counts for the same cells;
- if only standard gene counts are available, re-quantify BAM/FASTQ using one
  combined coding-gene plus cuTAR annotation.

For breast cancer:

- obtain the 10x invasive ductal carcinoma 3' reference used by the study,
  the Wu et al. breast cancer atlas, or another well-annotated compatible
  breast tumor scRNA reference;
- ensure the same cuTAR/lncRNA targets are quantified in that reference.

Without this addition, Curio and Visium data can validate spatial expression
but cannot complete the intended reference-to-space prediction workflow.

### Priority 2: Xenium validation data

Download the processed Xenium object/count matrix, cell coordinates,
segmentation/cell metadata, and the 76-cuTAR custom panel definition released
with the study. Xenium is the strongest targeted independent validation
resource because the paper reports direct single-cell in situ detection of the
selected cuTARs.

### Priority 3: metadata and source tables

Obtain:

- sample-to-cancer/platform metadata;
- Curio cell annotations and spatial coordinates;
- STOmics sample/tissue labels;
- MelD/MelDN tissue-region labels if available;
- cuTAR confidence/database-overlap metadata;
- Supplementary Table 1 dataset registry;
- Supplementary Table 5 Xenium target panel;
- source-data tables linking exemplar cuTARs across platforms.

Region annotations are required for Domain Enrichment Concordance.

### Priority 4: reference genome resources

Keep outside git:

- GRCh38 primary assembly FASTA matching the original analysis;
- GENCODE v43 for reproduction and v47 for current annotation comparison;
- a combined gene-plus-cuTAR GTF generated from a documented conversion;
- chromosome naming and coordinate-convention manifest.

### Priority 5: at least one non-SPanC development pair

Add at least one public paired scRNA-ST dataset with many shared genes for
general masked-gene method development. SPanC-Lnc should supply lncRNA-specific
and cross-platform evidence, while a conventional pair supplies a larger,
less-selected calibration benchmark.

## 5. Recommended dataset roles

| Role | Dataset |
|---|---|
| Pipeline smoke test | small subset of MelD or Curio observations and targets |
| Method development | conventional paired scRNA-ST dataset with masked low-expression genes |
| Melanoma internal lncRNA validation | melanoma scRNA reference -> MelD |
| Biological stress test | frozen model -> MelDN |
| Single-cell cross-platform test | breast or melanoma scRNA reference -> Curio/Seeker |
| Total-RNA external test | frozen model -> STOmics |
| Targeted highest-confidence test | frozen model -> Xenium 76-cuTAR panel |

## 6. Mandatory audit before model training

For every file, export only lightweight audit results:

- file and object dimensions;
- matrix type and count source;
- obs/var/layer/raw/obsm/uns schema;
- total counts and detected features;
- gene and cuTAR feature counts;
- duplicate IDs;
- barcode overlap between complementary objects;
- coordinate completeness;
- target overlap with `00_cuTARs.bed`;
- sample labels and biological state;
- memory required after subsetting.

The final keep/replace decision should only change if the audit reveals corrupt
matrices, missing coordinates, irreconcilable barcodes, or non-count values
without provenance.

## 7. Evidence sources

- Nature Methods SPanC-Lnc article:
  https://www.nature.com/articles/s41592-026-03071-4
- UQ eSpace processed-data record:
  https://doi.org/10.48610/4570a11
- SPanC-Lnc analysis repository:
  https://github.com/GenomicsMachineLearning/SPanc_Lnc_PanCancer_LncRNA_Atlas
