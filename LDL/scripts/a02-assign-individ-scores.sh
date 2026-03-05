#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -l h_vmem=50G

# take chromosome number as argument and pass to Rscript 

echo Started scoring on chr $1: `date +%Y-%m-%dT%H:%M:%S` 

Rscript --vanilla /home/tjy19/PRS/LDL/scripts/a02-assign-individ-scores.R $1

echo Finished scoring on chr $1: `date +%Y-%m-%dT%H:%M:%S` 
