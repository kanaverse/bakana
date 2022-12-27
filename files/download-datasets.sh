#!/bin/bash

set -e
set -u

base=https://github.com/jkanche/random-test-files/releases/download
dir=datasets
mkdir -p ${dir}

download() {
    output=${dir}/$2
    if [ ! -e $output ]
    then
        curl -L ${base}/$1 > $output
    fi
}

# Download PBMC datasets.
download 10x-pbmc-v1.0.0/pbmc3k-matrix.mtx.gz pbmc3k-matrix.mtx.gz
download 10x-pbmc-v1.0.0/pbmc3k-features.tsv.gz pbmc3k-features.tsv.gz
download 10x-pbmc-v1.0.0/pbmc3k-barcodes.tsv.gz pbmc3k-barcodes.tsv.gz
download 10x-pbmc-v1.0.0/pbmc4k-tenx.h5 pbmc4k-tenx.h5

download 10x-pbmc-v1.0.0/combined-matrix.mtx.gz pbmc-combined-matrix.mtx.gz
download 10x-pbmc-v1.0.0/combined-features.tsv.gz pbmc-combined-features.tsv.gz
download 10x-pbmc-v1.0.0/combined-barcodes.tsv.gz pbmc-combined-barcodes.tsv.gz

# Download Zeisel datasets.
download zeisel-brain-v1.0.0/csc.h5ad zeisel-brain.h5ad
download zeisel-brain-v1.0.0/csc.rds zeisel-brain.rds
download zeisel-brain-v1.0.0/dense_plus_positions.rds zeisel-brain-dense.rds

# Download immune datasets.
download 10x-immune-v1.0.0/immune_3.0.0_sub-matrix.mtx.gz immune_3.0.0-matrix.mtx.gz
download 10x-immune-v1.0.0/immune_3.0.0_sub-features.tsv.gz immune_3.0.0-features.tsv.gz
download 10x-immune-v1.0.0/immune_3.0.0_sub-barcodes.tsv.gz immune_3.0.0-barcodes.tsv.gz

download 10x-immune-v1.0.0/immune_3.0.0_sub-tenx.h5 immune_3.0.0-tenx.h5
download 10x-immune-v1.0.0/immune_3.0.0_sub-tenx.rds immune_3.0.0-tenx.rds

# Download CRISPR dataset.
download 10x-crispr-v1.0.0/crispr_6.0.0-tenx.h5 crispr_6.0.0-tenx.h5

download 10x-crispr-v1.0.0/crispr_6.0.0-matrix.mtx.gz crispr_6.0.0-matrix.mtx.gz
download 10x-crispr-v1.0.0/crispr_6.0.0-features.tsv.gz crispr_6.0.0-features.tsv.gz
download 10x-crispr-v1.0.0/crispr_6.0.0-barcodes.tsv.gz crispr_6.0.0-barcodes.tsv.gz

# Download the paul dataset.
download paul-hsc-v1.0.0/se.rds paul-hsc.rds
