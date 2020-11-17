#!/bin/bash

set -eu

pop=$1
sim=$2
chr=$3
n=${4:-'10000'}

if [[ $pop == 'EUR' ]]
then
  Ne=11418
elif [[ $pop == 'EAS' ]]
then
  Ne=14269
elif [[ $pop == 'AFR' ]]
then
  Ne=17469
elif [[ $pop == 'SAS' ]]
then
  Ne=14269
else
  echo "Unrecognized pop: $pop" > /dev/stderr
  exit 1
fi

sleep $sim

basedir=/humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping/1000GP_Phase3/chr$chr

mkdir -p $pop && cd $_

hapgen2 \
-h $basedir/1000GP_Phase3_chr${chr}.${pop}.haps \
-l $basedir/1000GP_Phase3_chr${chr}.${pop}.legend \
-m $basedir/../genetic_map_chr${chr}_combined_b37.txt \
-o ${pop}_sim${sim}_chr${chr}.gz \
-dl 60596 1 1 1 \
-n $n 0 \
-Ne $Ne

plink2 \
--gen ${pop}_sim${sim}_chr${chr}.controls.gen.gz ref-first \
--sample ${pop}_sim${sim}_chr${chr}.controls.sample \
--oxford-single-chr $chr \
--make-pgen \
--out ${pop}_sim${sim}_chr${chr}

# plink2 \
# --pfile ${pop}_sim${sim}_chr${chr} \
# --export vcf bgz id-paste=iid \
# --mac 1 \
# --out ${pop}_sim${sim}_chr${chr}
