#!/bin/bash

set -eu

pop=$1
sim=$2
chr=$3
chip=$4

basedir=/humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping/snplist

mkdir -p $pop/$chip && cd $_

plink2 \
--pfile ../${pop}_sim${sim}_chr${chr} \
--extract $basedir/$chip.snplist \
--mac 1 \
--export vcf-4.2 bgz id-paste=iid \
--out ${pop}_sim${sim}_chr${chr}_$chip
