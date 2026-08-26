# Polygenic Score Transferability for Lipid Traits in Samoan Cohorts

**Yapp TJ, Krishnan M, Liu S, Manna SL, Cheng H, Naseri T, Reupena MS, Viali S, Tuitele J, Deka R, Hawley NL, McGarvey ST, Weeks DE, Minster RL, Carlson JC**

Department of Human Genetics, University of Pittsburgh — and collaborating institutions

> **Manuscript:** "Variant Harmonization Critically Determines Polygenic Score Transferability for Lipid Traits in Samoan Populations" *(HGG Advances)*
> Correspondence: tjy19@pitt.edu

---

## Overview

This repository contains the analysis code for a systematic evaluation of multi-ancestry polygenic scores (PGS) for four lipid traits — LDL cholesterol (LDL-C), HDL cholesterol (HDL-C), triglycerides (TG), and total cholesterol (TC) — across 4,342 Samoan adults in five cohorts spanning 1990–2010. This represents the first lipid PGS benchmarking study in any Pacific Islander population.

**Key finding:** Apparent PGS transferability is critically dependent on variant harmonization. A curated pruning-and-thresholding LDL-C score (PGS000889; 9,009 variants) achieved only ~9% variant matching in Samoan imputed data and near-zero predictive performance, while a genome-wide PRS-CS score (PGS000888; ~1.24 million variants) achieved 99.6–99.7% matching and meaningful performance (incremental R²: 5.7–8.6%).

---

## Study Design

| Feature | Details |
|---|---|
| Traits | LDL-C, HDL-C, TG, TC |
| Cohorts | 5 cohorts (1990–2010); Samoa and American Samoa |
| Total N | 4,342 participants |
| PGS sources | Graham et al. 2021 (PGS000888); Kanoni et al. 2022 (PGS002781, PGS002783, PGS002784) |
| Genotyping | Affymetrix 6.0 (2010 cohort); Illumina GSA (earlier cohorts) |
| Imputation | Samoan-specific reference panel (Soifua Manuia; Carlson et al. 2025) |
| Performance metric | Incremental R² and partial R² with 95% bootstrapped CIs (n = 1,000) |
| Statistical model | Linear mixed models with kinship random effect (lmekin, coxme R package) |
| Analysis software | R v4.3.1 |

---

## Repository Structure

```
.
├── github_sanitized/
│   ├── pgs_validation_pipeline/       # End-to-end PGS scoring and evaluation pipeline
│   │   ├── scripts/
│   │   │   ├── 01_build_scoring_file_per_chr.R   # Harmonize PGS weights to GDS variants
│   │   │   ├── 02_score_per_chr.R                # Compute per-chromosome PRS
│   │   │   ├── 03_aggregate_prs.R                # Aggregate chromosome scores
│   │   │   ├── 04_fit_models_and_r2.R            # Fit models; compute R² with bootstrap CIs
│   │   │   ├── 05_extract_imputation_quality_gds.R  # Extract imputation quality metrics
│   │   │   └── 06_format_tables.R                # Format manuscript Tables 2 and 3
│   │   ├── config_template.yml                   # Configuration template (copy and edit)
│   │   └── README.md
│   │
│   └── pgs_validation_manuscript/     # Manuscript figure generation
│       ├── scripts/
│       │   ├── 01_generate_forest_plots.R             # Figure 1: incremental/partial R² forest plots
│       │   └── 02_generate_supplementary_scatter_density.R  # Figures S5–S8: PGS vs. phenotype scatter plots
│       ├── config_template.yml
│       └── README.md
```

---

## PGS Sources

Variant weights were downloaded from the [PGS Catalog](https://www.pgscatalog.org/) and lifted over from GRCh37 to GRCh38 using UCSC liftOver prior to harmonization.

| Trait | PGS ID | Score type | Variants | Reference |
|---|---|---|---|---|
| LDL-C | [PGS000888](https://www.pgscatalog.org/score/PGS000888/) | PRS-CS (genome-wide) | ~1.24M | Graham et al. 2021, *Nature* |
| HDL-C | [PGS002781](https://www.pgscatalog.org/score/PGS002781/) | Multi-ancestry meta-analysis | ~1.24M | Kanoni et al. 2022, *Genome Biol* |
| TC | [PGS002783](https://www.pgscatalog.org/score/PGS002783/) | Multi-ancestry meta-analysis | 10,699 | Kanoni et al. 2022, *Genome Biol* |
| TG | [PGS002784](https://www.pgscatalog.org/score/PGS002784/) | Multi-ancestry meta-analysis | 30,071 | Kanoni et al. 2022, *Genome Biol* |

> **Note on LDL-C score selection:** An alternative curated score (PGS000889; 9,009 variants) was initially evaluated but achieved only ~9% variant matching in Samoan imputed data, yielding near-zero predictive performance. All reported LDL-C results use PGS000888.

---

## Software Requirements

**R v4.3.1** with the following packages:

| Package | Use |
|---|---|
| `tidyverse` | Data manipulation and plotting throughout pipeline |
| `SeqArray` | Reading per-chromosome GDS genotype files (scripts 01, 05) |
| `SeqVarTools` | Variant-level utilities for GDS files (script 01) |
| `yaml` | Reading configuration files (all scripts) |
| `coxme` (v2.2-18) | `lmekin()` for linear mixed models with kinship (script 04; see note below) |
| `ggExtra` | Marginal density plots for supplementary figures |
| `scales`, `gridExtra`, `cowplot` | Forest plot layout and formatting |

> **Note on statistical modeling:** The published pipeline scripts (script 04) implement plain `lm()` as a portable reference implementation. The manuscript analyses used `lmekin()` from the `coxme` package to account for relatedness via an empirical kinship matrix. Applying this pipeline to related samples (as in the Samoan cohorts) requires replacing `lm()` with `lmekin()` and supplying a kinship matrix via the configuration.

Install all packages in R:
```r
install.packages(c("tidyverse", "yaml", "coxme", "ggExtra", "scales", "gridExtra", "cowplot"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SeqArray", "SeqVarTools"))
```

---

## Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/T-Yapp/samoan-lipid-pgs.git
cd samoan-lipid-pgs/github_sanitized/pgs_validation_pipeline

# 2. Copy and edit the configuration file
cp config_template.yml config.yml
# Edit config.yml: set paths to your GDS files, PGS weights, phenotype tables, and output directory

# 3. Run the pipeline (example: discovery cohort, chromosome 1)
Rscript scripts/01_build_scoring_file_per_chr.R config.yml discovery 1
Rscript scripts/02_score_per_chr.R config.yml discovery 1
# ... repeat for chromosomes 2–22 ...
Rscript scripts/03_aggregate_prs.R config.yml discovery

# 4. Evaluate performance
Rscript scripts/04_fit_models_and_r2.R config.yml
Rscript scripts/05_extract_imputation_quality_gds.R config.yml
Rscript scripts/06_format_tables.R config.yml

# 5. Generate manuscript figures
cd ../pgs_validation_manuscript
cp config_template.yml config.yml
# Edit config.yml to point to outputs from step 4
Rscript scripts/01_generate_forest_plots.R config.yml
Rscript scripts/02_generate_supplementary_scatter_density.R config.yml
```

See `github_sanitized/pgs_validation_pipeline/README.md` and `github_sanitized/pgs_validation_manuscript/README.md` for full input specifications.

---

## Data Availability

Genotype and phenotype data cannot be distributed directly due to participant privacy protections.

- **2010 Samoa cohort:** Available through dbGaP ([phs000914.v1.p1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000914.v1.p1) and [phs000972.v5.p1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000972.v5.p1))
- **Earlier cohorts (1990, 1994, 1995, 2002, 2003):** Available from the corresponding author upon reasonable request (tjy19@pitt.edu)
- **PGS weights:** Publicly available from the [PGS Catalog](https://www.pgscatalog.org/) (IDs listed above)

---

## Citation

> Yapp T-AJ, Krishnan M, Liu S, Manna SL, Cheng H, Naseri T, Reupena MS, Viali S, Tuitele J, Deka R, Hawley NL, McGarvey ST, Weeks DE, Minster RL, Carlson JC, Variant Harmonization Critically Determines Polygenic Score Transferability for Lipid Traits in Samoan Populations, Human Genetics and Genomics Advances (2026), doi: https://doi.org/10.1016/ j.xhgg.2026.100663.

---

## Funding

This work was supported by NHLBI grants R01HL093093, R01HL133040, and R01HL52611; NIA grant R01AG09375; and NIDDK grants R01DK59642 and R01DK55406. TOPMed molecular data were supported by the NHLBI.

---

## License

Code in this repository is provided for research reproducibility. Please cite the manuscript above if you use or adapt this pipeline.
