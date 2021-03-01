#!/bin/bash

set -eu

# default number of bins is 10
nbins=${1:-10}

for config in $(seq 16)
do 

python eval_meta_finemapping.py --config $config --nbins $nbins

done
