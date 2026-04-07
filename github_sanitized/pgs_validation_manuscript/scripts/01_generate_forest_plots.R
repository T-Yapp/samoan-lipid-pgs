#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(gridExtra)
  library(cowplot)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/01_generate_forest_plots.R <config.yml>")
}
cfg <- yaml::read_yaml(args[1])

out_dir <- cfg$io$output_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

r2 <- readr::read_csv(cfg$forest$input_r2_csv, show_col_types = FALSE)

label_map <- cfg$forest$cohort_label_map
cohort_levels <- unlist(cfg$forest$cohort_levels)
trait_order <- unlist(cfg$forest$trait_order)
cohort_cols <- unlist(cfg$forest$colors)

if ("cohort_label" %in% names(r2)) {
  r2 <- r2 %>%
    mutate(cohort_label = recode(cohort_label, !!!label_map))
}

r2 <- r2 %>%
  mutate(
    cohort_label = factor(cohort_label, levels = cohort_levels),
    Trait = factor(Trait, levels = trait_order)
  )

make_forest <- function(dat, metric = c("partial", "incremental"), x_max = 0.27) {
  metric <- match.arg(metric)
  r2_col <- if (metric == "partial") "Partial_R2" else "Incremental_R2"
  lci_col <- if (metric == "partial") "Partial_R2_LCI" else "Incremental_R2_LCI"
  uci_col <- if (metric == "partial") "Partial_R2_UCI" else "Incremental_R2_UCI"
  x_title <- if (metric == "partial") expression("Partial" ~ R^2) else expression("Incremental" ~ R^2)
  subtitle <- if (metric == "partial") {
    "Partial R2 (variance explained by PGS in adjusted residuals)"
  } else {
    "Incremental R2 (additional variance beyond covariates)"
  }

  x_breaks <- seq(0, x_max, by = 0.05)

  make_panel <- function(d, trait_name) {
    pd <- d %>%
      filter(Trait == trait_name, !is.na(.data[[r2_col]])) %>%
      mutate(cohort_label = factor(as.character(cohort_label), levels = rev(cohort_levels)))

    if (nrow(pd) == 0) return(NULL)

    ggplot(pd, aes(y = cohort_label, x = .data[[r2_col]], colour = cohort_label)) +
      geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
      geom_errorbar(aes(xmin = .data[[lci_col]], xmax = .data[[uci_col]]), width = 0.25, linewidth = 0.9) +
      geom_point(size = 4, shape = 18) +
      scale_x_continuous(
        limits = c(0, x_max),
        breaks = x_breaks,
        labels = number_format(accuracy = 0.01),
        expand = expansion(mult = c(0, 0.02))
      ) +
      scale_colour_manual(values = cohort_cols, drop = FALSE) +
      labs(title = trait_name, x = x_title, y = NULL) +
      theme_classic(base_size = 15) +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 13, colour = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.4),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 5, r = 10, b = 5, l = 15, unit = "pt")
      )
  }

  panels <- lapply(trait_order, function(tr) make_panel(dat, tr))
  panels <- Filter(Negate(is.null), panels)

  legend_df <- data.frame(cohort_label = factor(cohort_levels, levels = cohort_levels), x = 0, y = 0)
  legend_p <- ggplot(legend_df, aes(x = x, y = y, colour = cohort_label)) +
    geom_point(size = 4, shape = 18) +
    scale_colour_manual(values = cohort_cols, name = "Cohort") +
    theme_void(base_size = 14) +
    theme(
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 12),
      legend.key.size = unit(1.2, "lines")
    )
  legend_grob <- cowplot::get_legend(legend_p)

  grid.arrange(
    grobs = c(panels, list(legend_grob)),
    ncol = 2,
    top = grid::textGrob(subtitle, gp = grid::gpar(fontsize = 13, col = "grey30"))
  )
}

x_max <- min(0.27, max(r2$Partial_R2_UCI, na.rm = TRUE) * 1.1)

fp_partial <- make_forest(r2, "partial", x_max = x_max)
fp_partial_png <- file.path(out_dir, paste0(cfg$forest$out_prefix, "_partial_R2.png"))
fp_partial_pdf <- file.path(out_dir, paste0(cfg$forest$out_prefix, "_partial_R2.pdf"))
ggsave(fp_partial_png, fp_partial, width = cfg$forest$width_in, height = cfg$forest$height_in, dpi = cfg$forest$dpi)
ggsave(fp_partial_pdf, fp_partial, width = cfg$forest$width_in, height = cfg$forest$height_in)

fp_incr <- make_forest(r2, "incremental", x_max = x_max)
fp_incr_png <- file.path(out_dir, paste0(cfg$forest$out_prefix, "_incremental_R2.png"))
fp_incr_pdf <- file.path(out_dir, paste0(cfg$forest$out_prefix, "_incremental_R2.pdf"))
ggsave(fp_incr_png, fp_incr, width = cfg$forest$width_in, height = cfg$forest$height_in, dpi = cfg$forest$dpi)
ggsave(fp_incr_pdf, fp_incr, width = cfg$forest$width_in, height = cfg$forest$height_in)

message("Saved: ", fp_partial_png)
message("Saved: ", fp_partial_pdf)
message("Saved: ", fp_incr_png)
message("Saved: ", fp_incr_pdf)
