#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/06_format_tables.R <config.yml>")
}

cfg <- yaml::read_yaml(args[1])
out_dir <- cfg$io$output_root
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

r2_file <- file.path(out_dir, "all_r2_results.csv")
if (!file.exists(r2_file)) stop("Missing: ", r2_file)
r2 <- readr::read_csv(r2_file, show_col_types = FALSE)

table3 <- r2 %>%
  mutate(
    Partial_R2 = round(Partial_R2, 4),
    Incremental_R2 = round(Incremental_R2, 4),
    Partial_R2_CI = paste0("[", round(Partial_R2_LCI, 4), ", ", round(Partial_R2_UCI, 4), "]"),
    Incremental_R2_CI = paste0("[", round(Incremental_R2_LCI, 4), ", ", round(Incremental_R2_UCI, 4), "]")
  ) %>%
  select(Trait, cohort_label, N, Partial_R2, Partial_R2_CI, Incremental_R2, Incremental_R2_CI) %>%
  arrange(Trait, cohort_label)

readr::write_csv(table3, file.path(out_dir, "Table3_r2_summary.csv"))

imp_file <- file.path(out_dir, "imputation_quality_summary.csv")
variant_files <- Sys.glob(cfg$tables$variant_summary_glob)

if (file.exists(imp_file) && length(variant_files) > 0) {
  imp <- readr::read_csv(imp_file, show_col_types = FALSE)
  vars <- purrr::map_dfr(variant_files, ~ readr::read_csv(.x, show_col_types = FALSE), .id = "source_file")

  # Expects a cohort-like key in variant summary files
  cohort_col <- intersect(c("Cohort", "cohort_label", "cohort"), names(vars))
  if (length(cohort_col) > 0) {
    cc <- cohort_col[[1]]
    vars2 <- vars %>% rename(Cohort = all_of(cc))
    table2 <- vars2 %>%
      left_join(imp %>% select(Cohort, Mean_R2), by = "Cohort") %>%
      mutate(`Mean imputation quality (R2)` = Mean_R2) %>%
      select(any_of(c(
        "Trait", "Cohort", "Total_PGS_Variants", "Successfully_Matched",
        "Proportion_Matched", "Monomorphic_Variants", "Proportion_Monomorphic",
        "Mean_MAF", "Mean imputation quality (R2)"
      ))) %>%
      distinct()
    readr::write_csv(table2, file.path(out_dir, "Table2_variant_summary.csv"))
  }
}

message("Saved table outputs under: ", out_dir)
