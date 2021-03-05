#!/bin/bash

set -eu

for config in $(seq 16)
do 

python eval_meta_finemapping.py --config $config

done
