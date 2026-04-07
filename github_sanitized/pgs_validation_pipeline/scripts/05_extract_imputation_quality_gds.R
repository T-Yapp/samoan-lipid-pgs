#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(SeqArray)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/05_extract_imputation_quality_gds.R <config.yml>")
}

cfg <- yaml::read_yaml(args[1])
out_dir <- cfg$io$output_root
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

replace_chr <- function(template, chr) gsub("\\{chr\\}", as.character(chr), template)
chr_vec <- as.integer(unlist(cfg$scoring$chromosomes))
rsq_paths <- unlist(cfg$imputation_quality$rsq_paths)

extract_for_one <- function(scoring_dir, gds_template, label) {
  all_rsq <- c()
  n_total <- 0
  n_found <- 0

  for (chr_i in chr_vec) {
    sf <- file.path(scoring_dir, paste0("final.scoring.file.chr", chr_i, ".RData"))
    gf <- replace_chr(gds_template, chr_i)
    if (!file.exists(sf) || !file.exists(gf)) next

    load(sf)
    if (!"scoring_file" %in% ls()) next
    pgs_pos <- scoring_file$pos
    n_total <- n_total + length(pgs_pos)

    g <- seqOpen(gf)
    on.exit(seqClose(g), add = TRUE)
    pos <- seqGetData(g, "position")

    rsq <- NULL
    for (rp in rsq_paths) {
      rsq <- tryCatch(seqGetData(g, rp), error = function(e) NULL)
      if (!is.null(rsq)) break
    }
    if (is.null(rsq)) {
      seqClose(g)
      next
    }

    idx <- match(pgs_pos, pos)
    keep <- idx[!is.na(idx)]
    n_found <- n_found + length(keep)
    v <- rsq[keep]
    v <- v[!is.na(v)]
    all_rsq <- c(all_rsq, v)
    seqClose(g)
  }

  tibble(
    Cohort = label,
    Total_PGS_Variants = n_total,
    Variants_With_R2 = n_found,
    Proportion_With_R2 = ifelse(n_total > 0, paste0(round(100 * n_found / n_total, 1), "%"), "0%"),
    Mean_R2 = ifelse(length(all_rsq) > 0, round(mean(all_rsq), 3), NA_real_)
  )
}

rows <- list()
for (k in names(cfg$imputation_quality$cohorts)) {
  cc <- cfg$imputation_quality$cohorts[[k]]
  rows[[length(rows) + 1]] <- extract_for_one(cc$scoring_dir, cc$gds_template, cc$label)
}

out <- bind_rows(rows)
out_file <- file.path(out_dir, "imputation_quality_summary.csv")
readr::write_csv(out, out_file)
message("Saved imputation quality summary: ", out_file)
