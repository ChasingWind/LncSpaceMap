# Week 1 execution

Week 1 is split into two governed runs.

## Run W1-A: P0-P2 preparation

Inputs:

- finalized SPanC-Lnc combined reference H5AD;
- detailed finalized feature QC table;
- pinned GENCODE release 50 comprehensive CHR GTF;
- `configs/week1_spanc_lnc.yaml`.

W1-A validates the reference contract, joins gene biotypes by canonical
Ensembl ID, classifies cuTAR eligibility, selects low-expression
protein-coding proxy genes, and creates five deterministic gene-wise folds.

Large local outputs:

```text
processed/week1/
├── annotated_feature_qc.tsv.gz
├── target_catalog.tsv.gz
└── masked_gene_folds.tsv
```

Lightweight review outputs:

```text
git_eval/
├── logs/week1_prepare.log
├── metrics/week1_reference_contract.tsv
├── metrics/week1_annotation_summary.tsv
├── metrics/week1_target_eligibility_summary.tsv
├── metrics/week1_masked_gene_fold_summary.tsv
├── manifests/week1_prepare_manifest.json
└── manifests/week1_masked_gene_folds.tsv
```

Required terminal marker:

```text
PASS_WEEK1_P0_P2_READY_FOR_BASELINES
```

## Eligibility policy

A default eligible cuTAR must have:

- at least 20 total reference counts;
- detection in at least `max(10 cells, 0.2% of quantified cells)`;
- quantification support in at least two biological samples;
- chromosome, start, end, and strand.

Targets with some signal but insufficient support are `exploratory`. Targets
without reference signal are `insufficient_reference_signal`. No target is
deleted by W1-A.

## Masked proxy policy

Only GENCODE `protein_coding` genes can become proxy genes. Candidate
difficulty is measured from log total counts and detected-cell fraction
relative to eligible cuTARs. The closest 250 candidates are divided
deterministically into five gene-wise folds.

The held-out fold must later be excluded from anchors, embeddings, backend
fitting, and parameter selection. W1-A creates the split but does not yet run a
mapping backend.

## Run W1-B: baseline closure

W1-B begins only after W1-A review. It will:

1. audit and harmonize one spatial target;
2. construct leakage-free anchors for each fold;
3. run Tangram and SpaGE through canonical adapters;
4. write baseline predictions and diagnostics;
5. evaluate the first masked-gene run.
