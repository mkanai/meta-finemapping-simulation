#!/bin/bash

set -eu

pop=$1
sim=$2
chr=$3
pheno="${SGE_TASK_ID}"

mkdir -p /broad/hptmp/mkanai
basedir=/humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping
cd $basedir/new_simulations/true_pheno/

awk -v pop=$pop -v sim=$sim '
BEGIN {
    OFS = "\t"
}
NR == 1 {
    for (i = 1; i <= NF; i++) {
        col[$i] = i
    }
    print "variant", "ref", "alt", "beta"
}
NR > 1 {
    beta_col = "beta_"pop"_sim"sim
    print $col["variant"], $col["ref"], $col["alt"], $col[beta_col]
}
' pheno${pheno}.beta.txt > /broad/hptmp/mkanai/${pop}_sim${sim}_chr${chr}_pheno${pheno}.txt

plink2 \
--pfile $basedir/${pop}/${pop}_sim${sim}_chr${chr} \
--keep $basedir/${pop}/${pop}_sim${sim}_chr${chr}.king.cutoff.in.id \
--score /broad/hptmp/mkanai/${pop}_sim${sim}_chr${chr}_pheno${pheno}.txt 1 3 4 header no-mean-imputation \
        cols=maybefid,nallele,dosagesum,scoresums \
--out ${pop}_sim${sim}_chr${chr}_pheno${pheno}

Rscript $basedir/compute_true_pheno.R \
--pop $pop \
--sim $sim \
--chr $chr \
--pheno $pheno

bgzip -f ${pop}_sim${sim}_chr${chr}_pheno${pheno}.pheno