#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -l h_vmem=1G


# note: qsub having a problem with file names that start with numbers so I had to remove the prefixes
# reran because issue with creating variant include files the first time around

for i in {1..22} "X"
do 

  qsub /home/tjy19/PRS/LDL_2002/scripts/a02-assign-individ-scores.sh $i
done

