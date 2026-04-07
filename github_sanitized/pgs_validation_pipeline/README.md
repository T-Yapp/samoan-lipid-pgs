# PGS Validation Pipeline (GitHub Ready)

This folder contains a sanitized, configurable, end-to-end pipeline for the Samoan lipid PGS validation project intended for public release.

## Pipeline stages

1. `scripts/01_build_scoring_file_per_chr.R`  
   Harmonize PGS weights to cohort GDS variants and build per-chromosome scoring objects.
2. `scripts/02_score_per_chr.R`  
   Compute per-chromosome PRS values for all participants.
3. `scripts/03_aggregate_prs.R`  
   Merge chromosome-level PRS into one per-sample score.
4. `scripts/04_fit_models_and_r2.R`  
   Fit covariate and covariate+PRS models; compute partial and incremental R2 with bootstrap CIs.
5. `scripts/05_extract_imputation_quality_gds.R`  
   Extract quality metrics (R2) for variants used in scoring.
6. `scripts/06_format_tables.R`  
   Build manuscript-ready summary tables (Table 2 and Table 3 style outputs).

## Required inputs

Configured in `config_template.yml`:

- PGS weights file
- Cohort GDS templates
- Per-cohort phenotype/covariate/PRS-ready tables (or merged tables with `prs`)
- Variant summary files (optional for Table 2)

## Minimal run example

```bash
cp config_template.yml config.yml
# edit config.yml paths and column names

Rscript scripts/01_build_scoring_file_per_chr.R config.yml discovery 1
Rscript scripts/02_score_per_chr.R config.yml discovery 1
Rscript scripts/03_aggregate_prs.R config.yml discovery

Rscript scripts/04_fit_models_and_r2.R config.yml
Rscript scripts/05_extract_imputation_quality_gds.R config.yml
Rscript scripts/06_format_tables.R config.yml
```

## Notes for publication

- All paths are config-based (no absolute local paths in code).
- No personal identifiers are embedded.
- Outputs are written under `io.output_root`.
