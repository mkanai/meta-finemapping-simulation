#!/bin/bash

set -eu

config=$1
downsample=${2:-""}
pheno="${SGE_TASK_ID}"

basedir=/humgen/atgu1/methods/mkanai/workspace/202011_meta_finemapping
config_file=${basedir}/config/meta_config10_grch38.tsv

if [[ $downsample ]]
then
    suffix="_downsample${downsample}"
else
    suffix=""
fi

mkdir -p $basedir/new_simulations/meta_analysis_mamba_grch38${suffix} && cd $_

if [[ ! -s meta.pheno${pheno}.config${config}.tsv.gz ]]
then

Rscript $basedir/meta_analysis.R \
--config $config \
--pheno $pheno \
--config-file ${config_file} \
--sumstats-format "$basedir/new_simulations/assoc${suffix}/pheno%d/%s.PHENO.glm.linear.gz" \
--out-format "$basedir/new_simulations/meta_analysis_mamba_grch38${suffix}/meta.pheno%d.config%d.tsv" \
--plink-path $basedir/plink

fi