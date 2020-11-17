#!/bin/bash

set -eu

pop=$1
sim=$2
chr=$3
chip=$4
ref=$5

fname=${pop}_sim${sim}_chr${chr}
basedir=/humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping

for pheno in $(seq 100)
do

mkdir -p $basedir/assoc/pheno${pheno} && cd $_

plink2 \
--pfile $basedir/imputed/${fname}_${chip}_${ref} \
--maf 0.001 \
--pheno $basedir/true_pheno/${fname}_pheno${pheno}.pheno iid-only \
--pheno-name PHENO \
--covar $basedir/$pop/$fname.pca \
--covar-name PC1-PC10 \
--linear omit-ref hide-covar \
cols=a1count,nobs,orbeta,se,p,err \
--out ${fname}_pheno${pheno}_${chip}_${ref}

bgzip -f ${fname}_pheno${pheno}_${chip}_${ref}.PHENO.glm.linear

done
