
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Chromosome number must be supplied", call.=FALSE)
} else {
  chr.i <- args[1]
}

library(SeqArray)
library(tidyverse)
library(SeqVarTools)

#Step 0 - create scores object for PRS 
load("/home/tjy19/PRS/LDL_scores.RData")

strand_flip <- function(allele){
  if (allele == "A") {
    flip <- "T"
  } else if (allele == "C") {
    flip <- "G"
  } else if (allele == "G") {
    flip <- "C"
  } else if (allele == "T") {
    flip <- "A"
  }
  return(flip)
}


cat("\n\n")
cat(paste0("extracting variants for chromosome ", chr.i))
cat("\n")

# extract scoring info for that chr
tmp.scores <- scores %>% filter(chr_name == chr.i)

# extract just the PRS variants from this chr
#gds <- seqOpen(paste0("gds/samoa-chr", chr.i, ".dose_10a.gds"))
gds <- seqOpen( paste0("~/PRS/gds/discovery-9b-samoa-chr", chr.i, ".dose.gds"))  # freeze 9 

# drop non-SGIDs
sample.id <- seqGetData(gds, "sample.id")
keep.id <- sample.id[grep("^SG", sample.id)]
seqSetFilter(gds, sample.id=keep.id)

# extract variants by position in scoring file
pos <- seqGetData(gds, "position")
seqSetFilter(gds, variant.sel = (pos %in% tmp.scores$chromEnd.hg38))
pos <- seqGetData(gds, "position")

# extract alleles from PRS variants (format: "ref,alt")
alleles <- seqGetData(gds, "allele")
ref <- sapply(1:length(alleles), function(i) strsplit(alleles[i], ",")[[1]][1])
alt <- sapply(1:length(alleles), function(i) strsplit(alleles[i], ",")[[1]][2])

maf <- seqGetData(gds, "annotation/info/MAF")

samoa.vars <- tibble(var.id = seqGetData(gds, "variant.id"), # matches column name in dos (below)
                     chr = rep(chr.i, length(pos)), 
                     pos = pos,
                     ref = ref, 
                     alt = alt,
                     maf = maf)

#filter samora.vars. if missing ref or alt it needs to be dropped. then have it cat and say # rows dropped bc of missing alleles 
samoa.vars <- samoa.vars[!is.na(samoa.vars$ref) & !is.na(samoa.vars$alt),]
n_dropped <- nrow(samoa.vars[is.na(samoa.vars$ref) | is.na(samoa.vars$alt),])

cat(paste0("# variants dropped due to missing alleles on chr ", chr.i, ": ", n_dropped, "\n"))

# get genotypes (imputed dosages)
dos <- imputedDosage(gds, dosage.field="DS", use.names=TRUE) %>% as.data.frame() %>% rownames_to_column(., "ID") 
seqClose(gds)

# mean impute missingness
cols.to.impute <- which(colSums(is.na(dos))>0) %>% unname()
for(col.i in cols.to.impute) {
  dos[which(is.na(dos[, col.i])), col.i] <- mean(dos[,col.i], na.rm = TRUE)
}


# count how many variants found
tmp.scores$inGDS <- FALSE
tmp.scores$inGDS[tmp.scores$chromEnd.hg38 %in% pos] <- TRUE

cat(paste0("variants located on chr ", chr.i))
print(table(tmp.scores$inGDS))
cat("\n\n")

#variants dropped due tot missing alleles. put nrow. copy code above (maybe dont do this now)

### Harmonization between Samoan variants and PRS variants

# check allele agreement
scoring.file <- inner_join(samoa.vars, tmp.scores[, c("effect_allele", "other_allele", "effect_weight", "chromEnd.hg38")], by = c("pos" = "chromEnd.hg38"))

scoring.file$alleles.match <- FALSE
scoring.file$alleles.match[scoring.file$alt == scoring.file$effect_allele & scoring.file$ref == scoring.file$other_allele] <- TRUE
table(scoring.file$alleles.match)

# create column to indicate whether to flip dosage (2 -> 0, etc) when generating PRS
scoring.file$alleles.flip <- FALSE
scoring.file$alleles.flip[scoring.file$alt == scoring.file$other_allele & scoring.file$ref == scoring.file$effect_allele] <- TRUE

cat(paste0("allele flips on chr ", chr.i))
print(table(scoring.file$alleles.flip))
cat("\n\n")

# create column to indicate strand flips
scoring.file$flip.strand <- FALSE
for(i in 1:nrow(scoring.file)){
  if(!scoring.file$alleles.match[i] & !scoring.file$alleles.flip[i]){
    if (scoring.file$alt[i] == strand_flip(scoring.file$effect_allele[i]) & scoring.file$ref[i] == strand_flip(scoring.file$other_allele[i])){
      scoring.file$flip.strand[i] <- TRUE
    }
  }
}

cat(paste0("strand flips on chr ", chr.i))
print(table(scoring.file$flip.strand))
cat("\n\n")

# create column for a flipped strand with allele swap
scoring.file$flip.strand.swap.alleles <- FALSE
for(i in 1:nrow(scoring.file)){
  if(!scoring.file$alleles.match[i] & !scoring.file$alleles.flip[i] & !scoring.file$flip.strand[i]){
    if (scoring.file$alt[i] == strand_flip(scoring.file$other_allele[i]) & scoring.file$ref[i] == strand_flip(scoring.file$effect_allele[i])){
      scoring.file$flip.strand.swap.alleles[i] <- TRUE
    }
  }
}

cat(paste0("strand flips with allele swap on chr ", chr.i))
print(  table(scoring.file$flip.strand.swap.alleles))
cat("\n\n")

# flag to exclude a variant for mismatched alleles
scoring.file$drop <- FALSE
scoring.file$drop[!scoring.file$alleles.match & 
                    !scoring.file$alleles.flip & 
                    !scoring.file$flip.strand & 
                    !scoring.file$flip.strand.swap.alleles] <- TRUE

cat(paste0("variants with unreconcilable alleles to be dropped on chr ", chr.i))
print(table(scoring.file$drop))
cat("\n\n")


# create final scoring file for chr
final.scoring.file <- scoring.file %>% filter(!drop)

# subset dos to match the final scoring file
final.dos <- dos %>% dplyr::select(ID, as.character(final.scoring.file$var.id))

# check that they are lined up in the same variant order
cat(paste0("check variants align between scoring and dos files on chr ", chr.i))
print(table(names(final.dos)[2:ncol(final.dos)] == as.character(final.scoring.file$var.id)))
cat("\n\n") 

# Flip dosages for variants that have allele flips
cols.to.flip <- final.scoring.file %>% filter(alleles.flip | flip.strand.swap.alleles) %>% dplyr::pull(var.id) %>% as.character
final.dos[,cols.to.flip] <- (2 - final.dos[,cols.to.flip])

save(final.scoring.file, file=paste0("~/PRS/LDL/scoring-files/final.scoring.file.chr", chr.i, ".RData"))
save(final.dos, file=paste0("~/PRS/LDL/dosages/final.dos.chr", chr.i))


