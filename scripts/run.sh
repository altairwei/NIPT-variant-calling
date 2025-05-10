#!/usr/bin/env bash

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

if [ "$CONDA_DEFAULT_ENV" != "nipt" ]; then
    echo "Activating conda environment: nipt" 1>&2
    source activate nipt
else
    echo "Environment nipt is already activated." 1>&2
fi

for (( i=$1; i<50; i++ ))
do
    snakemake \
        --rerun-triggers input mtime \
        --cores all \
        --batch Module_1_Alignment_MultiQC=$i/50 \
        --ri \
        report
done
