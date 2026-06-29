#!/usr/bin/env Rscript

# Convenience wrapper to run manuscript figure scripts from the sanitized module.
# Usage:
#   Rscript scripts/07_make_manuscript_figures.R <manuscript_config.yml>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/07_make_manuscript_figures.R <manuscript_config.yml>")
}

cfg <- args[1]

cmd1 <- sprintf("Rscript ../pgs_validation_manuscript/scripts/01_generate_forest_plots.R %s", shQuote(cfg))
cmd2 <- sprintf("Rscript ../pgs_validation_manuscript/scripts/02_generate_supplementary_scatter_density.R %s", shQuote(cfg))

status1 <- system(cmd1)
status2 <- system(cmd2)

if (status1 != 0 || status2 != 0) {
  stop("One or more manuscript figure scripts failed.")
}

message("Manuscript figures generated successfully.")
