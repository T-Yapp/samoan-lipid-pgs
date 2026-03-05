#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Chromosome number must be supplied", call.=FALSE)
} else {
  chr.i <- args[1]
}

library(tidyverse)

load(paste0("~/PRS/LDL_1990/scoring-files/final.scoring.file.chr", chr.i, ".RData"))
load(paste0("~/PRS/LDL_1990/dosages/final.dos.chr", chr.i))

samoa.scores <- tibble(ID = final.dos$ID,
                       score = 0)

# For each person, extract the vector of dosages and multiple by the vector of effect weights
weights <- final.scoring.file %>% dplyr::pull(effect_weight)

for(person.i in 1:nrow(samoa.scores)){
  samoa.scores$score[person.i] <- sum(weights * final.dos[person.i, -1 ]) 
}

cat(paste0("Summary of PRS on chr ", chr.i))
cat("\n")
print(summary(samoa.scores$score))
cat("\n\n") 

save(samoa.scores, file=paste0("~/PRS/LDL_1990/prs-samoan-scores/samoa.prs.chr", chr.i, ".RData"))

