# Configuration

`default.yaml` records the planned initial defaults. Dataset-specific configs
will declare input paths, count sources, field mappings, genome build, and
dataset role. Paths remain local and must never point to committed raw data.

Every run must archive the resolved configuration outside git and copy a
redacted version or configuration hash into `git_eval/manifests`.
