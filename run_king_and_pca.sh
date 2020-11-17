#!/bin/bash

set -eu

pop=$1
sim=$2
chr=$3

cd $pop

fname=${pop}_sim${sim}_chr${chr}

# compute relatedness
plink2 \
--pfile $fname \
--make-king triangle bin \
--out $fname

# filter monozygotic twins, duplicated, or first-degree samples
plink2 \
--pfile $fname \
--king-cutoff $fname 0.177 \
--out $fname

# LD pruning for PCA
plink2 \
--pfile $fname \
--maf 0.01 \
--keep $fname.king.cutoff.in.id \
--indep-pairwise 1e3 1e3 0.5 \
--out $fname

# Run PCA on unrelated samples
plink2 \
--pfile $fname \
--extract $fname.prune.in \
--keep $fname.king.cutoff.in.id \
--freq counts \
--pca allele-wts 10 approx \
--out $fname

# Project PCA for all
plink2 \
--pfile $fname \
--read-freq $fname.acount \
--score $fname.eigenvec.allele 2 5 header-read no-mean-imputation \
        variance-standardize \
--score-col-nums 6-15 \
--out $fname

awk -v prefix=${pop}_sim${sim} '
NR == 1 {
  for (i = 1; i <= NF; i++) {
    col[$i] = i
  }
  header = "IID"
  for (i = 1; i <= 10; i++) {
    header = header"\tPC"i
  }
  print header
}
NR > 1 {
  line = prefix"_"(NR - 1)
  for (i = 1; i <= 10; i++) {
    line = line"\t"$col["PC"i"_AVG"]
  }
  print line
}
' $fname.sscore > $fname.pca
