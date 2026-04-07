#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(SeqArray)
  library(SeqVarTools)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript scripts/01_build_scoring_file_per_chr.R <config.yml> <cohort_key> <chr>")
}

cfg <- yaml::read_yaml(args[1])
cohort_key <- args[2]
chr_i <- as.integer(args[3])

if (!cohort_key %in% names(cfg$scoring$cohorts)) {
  stop("Unknown cohort key: ", cohort_key)
}

replace_chr <- function(template, chr) gsub("\\{chr\\}", as.character(chr), template)

strand_flip <- function(allele) {
  recode(allele, A = "T", T = "A", C = "G", G = "C", .default = NA_character_)
}

wcfg <- cfg$pgs$weights_columns
weights <- readr::read_tsv(cfg$pgs$weights_file, show_col_types = FALSE) %>%
  transmute(
    chr = as.character(.data[[wcfg$chr]]),
    pos = as.integer(.data[[wcfg$pos]]),
    effect_allele = as.character(.data[[wcfg$effect_allele]]),
    other_allele = as.character(.data[[wcfg$other_allele]]),
    effect_weight = as.numeric(.data[[wcfg$effect_weight]])
  )

tmp_scores <- weights %>% filter(chr == as.character(chr_i))
if (nrow(tmp_scores) == 0) stop("No PGS variants for chr ", chr_i)

cohort_cfg <- cfg$scoring$cohorts[[cohort_key]]
gds_file <- replace_chr(cohort_cfg$gds_template, chr_i)
if (!file.exists(gds_file)) stop("Missing GDS: ", gds_file)

out_root <- cohort_cfg$out_dir
scoring_dir <- file.path(out_root, "scoring-files")
dosage_dir <- file.path(out_root, "dosages")
dir.create(scoring_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(dosage_dir, showWarnings = FALSE, recursive = TRUE)

gds <- seqOpen(gds_file)
on.exit(seqClose(gds), add = TRUE)

sample_id <- seqGetData(gds, "sample.id")
keep_id <- sample_id[grepl(cfg$scoring$sample_id_regex, sample_id)]
seqSetFilter(gds, sample.id = keep_id)

pos_all <- seqGetData(gds, "position")
seqSetFilter(gds, variant.sel = pos_all %in% tmp_scores$pos)
pos <- seqGetData(gds, "position")

alleles <- seqGetData(gds, "allele")
ref <- sapply(strsplit(alleles, ","), `[`, 1)
alt <- sapply(strsplit(alleles, ","), `[`, 2)

maf <- tryCatch(seqGetData(gds, "annotation/info/MAF"), error = function(e) rep(NA_real_, length(pos)))
samoa_vars <- tibble(
  var.id = seqGetData(gds, "variant.id"),
  chr = chr_i,
  pos = pos,
  ref = ref,
  alt = alt,
  maf = maf
) %>%
  filter(!is.na(ref), !is.na(alt))

dos <- imputedDosage(gds, dosage.field = "DS", use.names = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("ID")

na_cols <- which(colSums(is.na(dos)) > 0)
if (length(na_cols) > 0) {
  for (j in na_cols) dos[[j]][is.na(dos[[j]])] <- mean(dos[[j]], na.rm = TRUE)
}

scoring_file <- inner_join(
  samoa_vars,
  tmp_scores %>% select(pos, effect_allele, other_allele, effect_weight),
  by = "pos"
)

scoring_file <- scoring_file %>%
  mutate(
    alleles.match = alt == effect_allele & ref == other_allele,
    alleles.flip = alt == other_allele & ref == effect_allele,
    flip.strand = (!alleles.match & !alleles.flip) &
      (alt == strand_flip(effect_allele) & ref == strand_flip(other_allele)),
    flip.strand.swap.alleles = (!alleles.match & !alleles.flip & !flip.strand) &
      (alt == strand_flip(other_allele) & ref == strand_flip(effect_allele)),
    drop = !(alleles.match | alleles.flip | flip.strand | flip.strand.swap.alleles)
  ) %>%
  filter(!drop)

final_dos <- dos %>% select(ID, all_of(as.character(scoring_file$var.id)))
flip_cols <- scoring_file %>%
  filter(alleles.flip | flip.strand.swap.alleles) %>%
  pull(var.id) %>%
  as.character()
if (length(flip_cols) > 0) final_dos[, flip_cols] <- 2 - final_dos[, flip_cols]

save(scoring_file, file = file.path(scoring_dir, paste0("final.scoring.file.chr", chr_i, ".RData")))
save(final_dos, file = file.path(dosage_dir, paste0("final.dos.chr", chr_i, ".RData")))

message("Saved chr ", chr_i, " scoring objects for ", cohort_key)
