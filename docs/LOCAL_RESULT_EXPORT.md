# Local-only model comparison export

Detailed SpaGE/Tangram comparison artifacts are generated under
`import_result/`. The complete directory is ignored by git and must not be
copied into `git_eval`.

Run from the repository root:

```bash
python scripts/plot_week1_model_comparison.py
```

Outputs:

- `model_comparison_per_gene.tsv`: paired per-target metrics and winners.
- `model_comparison_summary.tsv`: model-level summary.
- `model_comparison_by_detection.tsv`: results stratified by target support.
- `model_comparison_decision.md`: reusable interpretation and model choice.
- `model_comparison_4panel.png`: 300 dpi raster figure.
- `model_comparison_4panel.pdf`: editable vector figure.

The script automatically combines every fold for which both backend metric
files exist. Running it now summarizes fold 0; rerunning after the remaining
folds are available creates the full five-fold comparison without changing the
command.
