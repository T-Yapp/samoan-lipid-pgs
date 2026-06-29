#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript scripts/03_aggregate_prs.R <config.yml> <cohort_key>")
}

cfg <- yaml::read_yaml(args[1])
cohort_key <- args[2]

if (!cohort_key %in% names(cfg$scoring$cohorts)) {
  stop("Unknown cohort key: ", cohort_key)
}

out_root <- cfg$scoring$cohorts[[cohort_key]]$out_dir
prs_dir <- file.path(out_root, "prs-scores")
chr_vec <- as.integer(unlist(cfg$scoring$chromosomes))

all_scores <- NULL
for (k in chr_vec) {
  f <- file.path(prs_dir, paste0("samoa.prs.chr", k, ".RData"))
  if (!file.exists(f)) next
  load(f)
  names(samoa.scores)[2] <- paste0("s", k)
  if (is.null(all_scores)) {
    all_scores <- samoa.scores
  } else {
    all_scores <- left_join(all_scores, samoa.scores, by = "ID")
  }
}

if (is.null(all_scores)) stop("No per-chromosome PRS files found in ", prs_dir)

all_scores <- all_scores %>%
  mutate(prs = rowSums(select(., starts_with("s")), na.rm = TRUE)) %>%
  select(ID, prs, starts_with("s"))

out_file <- file.path(out_root, "prs_total.tsv")
readr::write_tsv(all_scores, out_file)
message("Saved aggregated PRS: ", out_file)
