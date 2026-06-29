#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript scripts/02_score_per_chr.R <config.yml> <cohort_key> <chr>")
}

cfg <- yaml::read_yaml(args[1])
cohort_key <- args[2]
chr_i <- as.integer(args[3])

if (!cohort_key %in% names(cfg$scoring$cohorts)) {
  stop("Unknown cohort key: ", cohort_key)
}

out_root <- cfg$scoring$cohorts[[cohort_key]]$out_dir
scoring_file <- file.path(out_root, "scoring-files", paste0("final.scoring.file.chr", chr_i, ".RData"))
dos_file <- file.path(out_root, "dosages", paste0("final.dos.chr", chr_i, ".RData"))
prs_dir <- file.path(out_root, "prs-scores")
dir.create(prs_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(scoring_file) || !file.exists(dos_file)) {
  stop("Missing required inputs for chr ", chr_i, " in ", cohort_key)
}

load(scoring_file)
load(dos_file)

weights <- scoring_file$effect_weight
dos_mat <- as.matrix(final_dos[, -1, drop = FALSE])
score <- as.numeric(dos_mat %*% weights)

samoa.scores <- tibble(ID = final_dos$ID, score = score)
save(samoa.scores, file = file.path(prs_dir, paste0("samoa.prs.chr", chr_i, ".RData")))

message("Saved chr ", chr_i, " PRS for ", cohort_key)
