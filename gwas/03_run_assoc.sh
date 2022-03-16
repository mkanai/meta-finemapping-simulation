#!/bin/bash

set -eu

pop=$1
sim=$2
chr=$3
chip=$4
ref=$5
downsample=${6:-""}
pheno="${SGE_TASK_ID}"

fname=${pop}_sim${sim}_chr${chr}
basedir=/humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping
liftover_variants=$basedir/imputed/TopMed_GRCh38.variants.liftover.txt

if [[ $downsample ]]
then
    suffix="_downsample${downsample}"
    keep="--keep $basedir/${pop}/downsample/${pop}_sim${sim}_chr${chr}.downsample${downsample}.keep"
else
    suffix=""
    keep=""
fi


mkdir -p $basedir/new_simulations/assoc${suffix}/pheno${pheno} && cd $_

if [[ ! -s ${fname}_pheno${pheno}_${chip}_${ref}.PHENO.glm.linear.gz ]]
then

    if [[ -s $basedir/imputed/${fname}_${chip}_${ref}.exclude ]]
    then
        exclude="--exclude $basedir/imputed/${fname}_${chip}_${ref}.exclude"
    else
        exclude=""
    fi

    plink2 \
    --pfile $basedir/imputed/${fname}_${chip}_${ref} \
    --maf 0.001 \
    $exclude \
    $keep \
    --pheno $basedir/new_simulations/true_pheno/${fname}_pheno${pheno}.pheno.gz iid-only \
    --pheno-name PHENO \
    --covar $basedir/$pop/$fname.pca \
    --covar-name PC1-PC10 \
    --linear omit-ref hide-covar \
    cols=chrom,pos,ref,alt,a1count,nobs,orbeta,se,p,err \
    --out ${fname}_pheno${pheno}_${chip}_${ref}

    awk '
    BEGIN {
        OFS = "\t"
    }
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            col[$i] = i
        }
        print "CHR", "BP", "SNP", "A1", "A2", "A1_CT", "OBS_CT", "BETA", "SE", "P"
    }
    NR > 1 {
        if ($col["ALT"] == $col["A1"]) {
            a2 = $col["REF"]
        } else {
            a2 = $col["ALT"]
        }
        print $col["#CHROM"], $col["POS"], $col["ID"], $col["A1"], a2, $col["A1_CT"], $col["OBS_CT"], $col["BETA"], $col["SE"], $col["P"]
    }
    ' ${fname}_pheno${pheno}_${chip}_${ref}.PHENO.glm.linear | bgzip -c > ${fname}_pheno${pheno}_${chip}_${ref}.PHENO.glm.linear.gz

    rm ${fname}_pheno${pheno}_${chip}_${ref}.PHENO.glm.linear
fi

if [[ $ref == "TopMed" && (! -s ${fname}_pheno${pheno}_${chip}_${ref}-liftover.PHENO.glm.linear.gz) ]]
then
    echo -n "Lifting over... "
    awk '
    BEGIN {
        OFS = "\t"
    }
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            col[$i] = i
        }
    }
    FNR == NR && FNR > 1 {
        variant[$col["original_variant"]] = $col["variant"]
        ref_alt_flip[$col["original_variant"]] = $col["ref_alt_flip"]
        next
    }
    FNR < NR && FNR == 1 {
        for (i = 1; i <= NF; i++) {
            col[$i] = i
        }
        print $0, "ref_alt_flip", "original_variant"
    }
    FNR < NR && FNR > 1 && $col["SNP"] in variant {
        original_variant = $col["SNP"]
        flip = ref_alt_flip[original_variant]
        new_variant = variant[original_variant]
        split(new_variant, a, ":")
        if ($col["A1"] != a[4]) {
            $col["BETA"] = -$col["BETA"]
            $col["A1_CT"] = 2 * $col["OBS_CT"] - $col["A1_CT"]
        }
        $col["SNP"] = new_variant        
        $col["CHR"] = a[1]
        $col["BP"] = a[2]
        $col["A2"] = a[3]
        $col["A1"] = a[4]
        print $0, flip, original_variant | "sort -k2,2 -n"
    }
    ' ${liftover_variants} <(zcat ${fname}_pheno${pheno}_${chip}_${ref}.PHENO.glm.linear.gz) \
    | bgzip -c > ${fname}_pheno${pheno}_${chip}_${ref}-liftover.PHENO.glm.linear.gz
    echo "Done."
fi