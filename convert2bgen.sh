#!/bin/bash

set -eu

fname=$1

plink2 \
--vcf $fname.dose.vcf.gz dosage=DS \
--extract-if-info 'R2 > 0.6' \
--make-pgen \
--out $fname

cp $fname.psam $fname.psam.bak
awk -v fname=$fname '
BEGIN {
  OFS = "\t"
}
NR == 1 {
  print
}
NR > 1 {
  split(fname, a, "_")
  print a[1]"_"a[2]"_"(NR - 1), $2
}
' $fname.psam.bak > $fname.psam

plink2 \
--pfile $fname \
--export bgen-1.2 ref-first \
--out $fname

bgenix -g $fname.bgen -index
