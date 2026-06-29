#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/04_fit_models_and_r2.R <config.yml>")
}

cfg <- yaml::read_yaml(args[1])
out_dir <- cfg$io$output_root
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

covars <- unlist(cfg$modeling$covariates)
b_n <- as.integer(cfg$modeling$bootstrap$n)
set.seed(as.integer(cfg$modeling$bootstrap$seed))

boot_ci <- function(df, y_col, prs_col, covars, b_n = 1000) {
  base_f <- as.formula(paste(y_col, "~", paste(covars, collapse = " + ")))
  full_f <- as.formula(paste(y_col, "~", paste(c(covars, prs_col), collapse = " + ")))

  partial_v <- numeric(b_n)
  incr_v <- numeric(b_n)

  n <- nrow(df)
  for (i in seq_len(b_n)) {
    idx <- sample.int(n, n, replace = TRUE)
    d <- df[idx, , drop = FALSE]
    m0 <- lm(base_f, data = d)
    m1 <- lm(full_f, data = d)
    r0 <- summary(m0)$r.squared
    r1 <- summary(m1)$r.squared
    partial_v[i] <- r1
    incr_v[i] <- r1 - r0
  }

  list(
    partial_lci = unname(quantile(partial_v, 0.025, na.rm = TRUE)),
    partial_uci = unname(quantile(partial_v, 0.975, na.rm = TRUE)),
    incr_lci = unname(quantile(incr_v, 0.025, na.rm = TRUE)),
    incr_uci = unname(quantile(incr_v, 0.975, na.rm = TRUE))
  )
}

results <- list()
trait_cfgs <- cfg$modeling$trait_configs

for (trait_name in names(trait_cfgs)) {
  tc <- trait_cfgs[[trait_name]]
  d <- readr::read_csv(tc$input_csv, show_col_types = FALSE)

  y_col <- tc$phenotype_col
  prs_col <- tc$prs_col
  cohort_col <- tc$cohort_col

  req <- c(y_col, prs_col, cohort_col, covars)
  d <- d %>% filter(complete.cases(across(all_of(req))))
  if (isTRUE(tc$log_transform)) d[[y_col]] <- log(d[[y_col]])

  cohorts <- sort(unique(as.character(d[[cohort_col]])))
  for (cc in cohorts) {
    dc <- d %>% filter(.data[[cohort_col]] == cc)
    if (nrow(dc) < 50) next

    base_f <- as.formula(paste(y_col, "~", paste(covars, collapse = " + ")))
    full_f <- as.formula(paste(y_col, "~", paste(c(covars, prs_col), collapse = " + ")))
    m0 <- lm(base_f, data = dc)
    m1 <- lm(full_f, data = dc)
    r0 <- summary(m0)$r.squared
    r1 <- summary(m1)$r.squared
    ci <- boot_ci(dc, y_col, prs_col, covars, b_n = b_n)

    results[[length(results) + 1]] <- tibble(
      Trait = trait_name,
      cohort_label = cc,
      N = nrow(dc),
      Partial_R2 = r1,
      Partial_R2_LCI = ci$partial_lci,
      Partial_R2_UCI = ci$partial_uci,
      Incremental_R2 = r1 - r0,
      Incremental_R2_LCI = ci$incr_lci,
      Incremental_R2_UCI = ci$incr_uci
    )
  }
}

out <- bind_rows(results) %>% arrange(Trait, cohort_label)
out_file <- file.path(out_dir, "all_r2_results.csv")
readr::write_csv(out, out_file)
message("Saved R2 summary: ", out_file)
