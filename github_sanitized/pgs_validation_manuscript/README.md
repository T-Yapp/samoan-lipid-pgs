# PGS Validation Manuscript Pipeline (Sanitized)

This folder contains a path-free, shareable version of the figure-generation code used for the lipid PGS validation manuscript.

## What is included

- `config_template.yml`  
  Configuration file with all input/output paths and labels.
- `scripts/01_generate_forest_plots.R`  
  Generates manuscript forest plots (partial and incremental R2).
- `scripts/02_generate_supplementary_scatter_density.R`  
  Generates supplementary S1-S8 scatter+density plots (raw and residual versions).

## Inputs expected

### Forest plot script

- A CSV like `all_r2_results_PGS000888.csv` with columns:
  - `Trait`
  - `cohort_label`
  - `Partial_R2`, `Partial_R2_LCI`, `Partial_R2_UCI`
  - `Incremental_R2`, `Incremental_R2_LCI`, `Incremental_R2_UCI`

### Supplementary scatter script

Trait-specific PRS/phenotype tables with at least:

- `ID`
- `prs`
- phenotype column (configured per trait)
- covariates: `age`, `age2`, `sex`, `sex_age`, `sex_age2`, `V1`, `V2`, `V3`

For LDL replication (if score-by-chromosome files are used), each chromosome file should contain an object named `samoa.scores` with columns:

- `ID`
- one score column (per chromosome)

## Usage

1. Copy and edit `config_template.yml` (for example, `config.yml`).
2. Run:

```bash
Rscript scripts/01_generate_forest_plots.R config.yml
Rscript scripts/02_generate_supplementary_scatter_density.R config.yml
```

## Notes

- No absolute paths are hardcoded.
- No user-specific identifiers are embedded.
- TG is log-transformed before modeling in the supplementary script (to match analysis models).
