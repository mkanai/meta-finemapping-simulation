#!/bin/bash

set -eu

pop=$1
sim=$2
chr=$3

cd true_pheno/

for pheno in $(seq 100)
do

plink2 \
--pfile ../${pop}/${pop}_sim${sim}_chr${chr} \
--keep ../${pop}/${pop}_sim${sim}_chr${chr}.king.cutoff.in.id \
--score pheno${pheno}.beta.txt 1 3 4 header no-mean-imputation \
        cols=maybefid,nallele,dosagesum,scoresums \
--out ${pop}_sim${sim}_chr${chr}_pheno${pheno}

Rscript ../compute_true_pheno.R \
--pop $pop \
--sim $sim \
--chr $chr \
--pheno $pheno

done