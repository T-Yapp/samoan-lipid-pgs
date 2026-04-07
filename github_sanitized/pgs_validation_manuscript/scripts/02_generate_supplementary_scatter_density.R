#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggExtra)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/02_generate_supplementary_scatter_density.R <config.yml>")
}
cfg <- yaml::read_yaml(args[1])

out_dir <- cfg$io$output_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cohort_levels <- unlist(cfg$supplementary$cohort_levels)
cohort_cols <- unlist(cfg$supplementary$colors)
covars <- unlist(cfg$supplementary$covariates)
trait_names <- names(cfg$supplementary$traits)

add_partial_resid <- function(df, pheno_col, covariates, log_transform = FALSE) {
  req <- c(pheno_col, covariates, "prs")
  df <- df[complete.cases(df[, req]), ]
  if (log_transform) {
    df[[pheno_col]] <- log(df[[pheno_col]])
  }
  form <- as.formula(paste(pheno_col, "~", paste(covariates, collapse = " + ")))
  df$partial_resid <- residuals(lm(form, data = df))
  df$raw_phenotype <- df[[pheno_col]]
  df$prs_z <- as.numeric(scale(df$prs))
  df
}

load_chr_sum_prs <- function(score_dir) {
  all_s <- NULL
  for (k in 1:22) {
    rdata <- file.path(score_dir, paste0("samoa.prs.PGS000888.chr", k, ".RData"))
    if (!file.exists(rdata)) next
    load(rdata)
    if (!exists("samoa.scores")) {
      stop("Object 'samoa.scores' not found in: ", rdata)
    }
    names(samoa.scores)[2] <- paste0("s", k)
    if (is.null(all_s)) {
      all_s <- samoa.scores
    } else {
      all_s <- dplyr::left_join(all_s, samoa.scores, by = "ID")
    }
  }
  if (is.null(all_s)) {
    stop("No chromosome score files found in: ", score_dir)
  }
  all_s$prs <- rowSums(all_s[, setdiff(names(all_s), "ID")], na.rm = TRUE)
  all_s[, c("ID", "prs")] %>% mutate(ID = as.character(ID))
}

load_one_cohort <- function(file_path, pheno_col, covariates, log_transform, cohort_label) {
  read_tsv(file_path, show_col_types = FALSE) %>%
    rename(phenotype = all_of(pheno_col)) %>%
    select(ID, prs, all_of(covariates), phenotype) %>%
    mutate(ID = as.character(ID)) %>%
    add_partial_resid("phenotype", covariates, log_transform = log_transform) %>%
    mutate(cohort = cohort_label)
}

load_ldl_chr_cohort <- function(score_dir, covar_file, pheno_col, covariates, cohort_label, log_transform = FALSE) {
  prs_df <- load_chr_sum_prs(score_dir)
  covar_df <- read_tsv(covar_file, show_col_types = FALSE) %>%
    select(ID, all_of(covariates), all_of(pheno_col)) %>%
    mutate(ID = as.character(ID))
  inner_join(prs_df, covar_df, by = "ID") %>%
    rename(phenotype = all_of(pheno_col)) %>%
    add_partial_resid("phenotype", covariates, log_transform = log_transform) %>%
    mutate(cohort = cohort_label)
}

build_trait_df <- function(trait_name, trait_cfg) {
  discovery <- load_one_cohort(
    file_path = trait_cfg$discovery_file,
    pheno_col = trait_cfg$discovery_pheno_col,
    covariates = covars,
    log_transform = isTRUE(trait_cfg$log_transform),
    cohort_label = "2010 Samoa"
  )

  replication <- list()
  if (isTRUE(trait_cfg$use_chr_scores)) {
    for (cl in names(trait_cfg$replication_chr_score_dirs)) {
      replication[[cl]] <- load_ldl_chr_cohort(
        score_dir = trait_cfg$replication_chr_score_dirs[[cl]],
        covar_file = trait_cfg$replication_covar_files[[cl]],
        pheno_col = trait_cfg$replication_pheno_col,
        covariates = covars,
        cohort_label = cl,
        log_transform = isTRUE(trait_cfg$log_transform)
      )
    }
  } else {
    for (cl in names(trait_cfg$replication_files)) {
      replication[[cl]] <- load_one_cohort(
        file_path = trait_cfg$replication_files[[cl]],
        pheno_col = trait_cfg$replication_pheno_col,
        covariates = covars,
        log_transform = isTRUE(trait_cfg$log_transform),
        cohort_label = cl
      )
    }
  }

  bind_rows(c(list(discovery), replication)) %>%
    mutate(cohort = factor(cohort, levels = cohort_levels))
}

make_fig <- function(df, trait_name, y_label, out_stem, mode = "residual", log_transform = FALSE) {
  ns <- df %>% count(cohort) %>% mutate(lab = paste0(cohort, "\n(n=", n, ")"))
  df <- left_join(df, ns[, c("cohort", "lab")], by = "cohort") %>%
    mutate(lab = factor(lab, levels = ns$lab))
  col_map <- setNames(cohort_cols[as.character(ns$cohort)], ns$lab)

  y_units <- if (log_transform) "log(mg/dL)" else "mg/dL"
  y_raw_lbl <- if (log_transform) paste0("log(", y_label, ") (", y_units, ")") else paste0(y_label, " (", y_units, ")")

  if (mode == "residual") {
    y_var <- "partial_resid"
    y_title <- if (log_transform) paste0("log(", y_label, ") partial residuals") else paste0(y_label, " partial residuals (", y_units, ")")
    caption <- paste0(
      if (log_transform) "Trait was log-transformed prior to analysis. " else "",
      "Partial residuals after adjusting for age, age2, sex, age x sex, age2 x sex, and PC1-PC3.\n",
      "Regression line: lm(partial residuals ~ standardized PGS), per cohort."
    )
  } else {
    y_var <- "raw_phenotype"
    y_title <- y_raw_lbl
    caption <- paste0(
      if (log_transform) "Trait was log-transformed prior to analysis; y-axis shows log-transformed values. " else "",
      "Observed values vs. polygenic score. Regression lines are unadjusted."
    )
  }

  p <- ggplot(df, aes(x = prs_z, y = .data[[y_var]], colour = lab)) +
    geom_point(alpha = 0.18, size = 0.7) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, aes(fill = lab), alpha = 0.15) +
    scale_colour_manual(values = col_map, name = "Cohort") +
    scale_fill_manual(values = col_map, name = "Cohort") +
    labs(
      title = paste0(trait_name, " PGS vs Lipid Levels Across Samoan Cohorts"),
      x = "Polygenic Score (standardized)",
      y = y_title,
      caption = caption
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 13, face = "bold"),
      axis.title = element_text(size = 11),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.caption = element_text(size = 8, colour = "grey40", hjust = 0),
      legend.position = "right"
    )

  p_final <- ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE, alpha = 0.3, size = 5)

  out_png <- file.path(out_dir, paste0(out_stem, "_", mode, ".png"))
  out_pdf <- file.path(out_dir, paste0(out_stem, "_", mode, ".pdf"))
  ggsave(out_png, p_final, width = 10, height = 7, dpi = 300)
  ggsave(out_pdf, p_final, width = 10, height = 7)
  message("Saved: ", out_png)
  message("Saved: ", out_pdf)
}

for (trait_name in trait_names) {
  trait_cfg <- cfg$supplementary$traits[[trait_name]]
  y_label <- trait_cfg$y_label
  out_stem <- cfg$supplementary$out_prefix[[trait_name]]
  log_transform <- isTRUE(trait_cfg$log_transform)

  message("Building trait: ", trait_name)
  df <- build_trait_df(trait_name, trait_cfg)

  make_fig(df, trait_name, y_label, out_stem, mode = "raw", log_transform = log_transform)
  make_fig(df, trait_name, y_label, out_stem, mode = "residual", log_transform = log_transform)
}
