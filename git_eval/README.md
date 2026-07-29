# git_eval review exchange

This directory is the only place for lightweight outputs copied from an
external compute server into git for review.

## Allowed

- `metrics/`: CSV or TSV summaries.
- `logs/`: concise run or error logs.
- `figures/`: PNG, JPEG, SVG, or PDF figures.
- `reports/`: Markdown or HTML review reports.
- `manifests/`: JSON, YAML, or TXT run/configuration summaries.

Recommended maximums:

- less than 50 MB per file;
- less than 200 MB per reviewed run;
- plots at review resolution, not full microscopy resolution.

## Never commit

- H5AD, H5, Loom, RDS, Zarr, MTX, or Parquet matrices;
- FASTQ, BAM, CRAM, BED, GTF, GFF, or genome FASTA files;
- model checkpoints or serialized Python/R objects;
- complete Space Ranger, Cell Ranger, Seeker, or SAW output directories;
- tissue images at raw resolution;
- patient-identifying metadata or unredacted absolute server paths.

## Minimum review bundle

```text
git_eval/
├── manifests/
│   ├── run_manifest.json
│   └── resolved_config.yaml
├── metrics/
│   ├── dataset_audit.tsv
│   ├── metrics_summary.tsv
│   ├── metrics_by_gene.tsv
│   └── resource_usage.tsv
├── logs/
│   └── run.log
├── figures/
│   ├── benchmark_overview.png
│   ├── confidence_calibration.png
│   └── selected_spatial_maps.pdf
└── reports/
    └── report.md
```

Use names such as `YYYYMMDD_dataset_stage_seed.ext`. Redact local paths and
sensitive sample identifiers before committing.
