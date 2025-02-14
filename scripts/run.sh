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

snakemake \
    --rerun-triggers input mtime \
    --cores all \
    "$@"
