#!/bin/bash

set -e
set -u

base=https://github.com/kanaverse/random-test-files/releases/download
dir=datasets
mkdir -p ${dir}

download() {
    local output=${dir}/$2
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
download zeisel-brain-v1.0.0/with_results.rds zeisel-brain-with-results.rds
download zeisel-brain-v1.0.0/with_results.h5ad zeisel-brain-with-results.h5ad

# Download immune datasets.
download 10x-immune-v1.0.0/immune_3.0.0_sub-matrix.mtx.gz immune_3.0.0-matrix.mtx.gz
download 10x-immune-v1.0.0/immune_3.0.0_sub-features.tsv.gz immune_3.0.0-features.tsv.gz
download 10x-immune-v1.0.0/immune_3.0.0_sub-barcodes.tsv.gz immune_3.0.0-barcodes.tsv.gz

download 10x-immune-v1.0.0/immune_3.0.0_sub-tenx.h5 immune_3.0.0-tenx.h5
download 10x-immune-v1.0.0/immune_3.0.0_sub-tenx.rds immune_3.0.0-tenx.rds

# Download CRISPR dataset.
download 10x-crispr-v1.0.0/crispr_6.0.0_sub-tenx.h5 crispr_6.0.0-tenx.h5

download 10x-crispr-v1.0.0/crispr_6.0.0_sub-matrix.mtx.gz crispr_6.0.0-matrix.mtx.gz
download 10x-crispr-v1.0.0/crispr_6.0.0_sub-features.tsv.gz crispr_6.0.0-features.tsv.gz
download 10x-crispr-v1.0.0/crispr_6.0.0_sub-barcodes.tsv.gz crispr_6.0.0-barcodes.tsv.gz

# Download the paul dataset.
download paul-hsc-v1.0.0/se.rds paul-hsc.rds

download_and_untar() {
    local input=$1
    local fname=$(basename $input)
    local output=$2/$fname
    local trunc=$(echo $output | sed "s/.tar.gz$//")
    if [ ! -e ${dir}/${trunc} ]
    then
        download $input $output
        (cd ${dir}/$2 && tar -xvf $fname)
        rm ${dir}/$output
    fi
}

mkdir -p ${dir}/alabaster
download_and_untar alabaster-v1.0.0/zeisel-brain-dense.tar.gz alabaster
download_and_untar alabaster-v1.0.0/zeisel-brain-sparse.tar.gz alabaster
download_and_untar alabaster-v1.0.0/zeisel-brain-stripped.tar.gz alabaster
download_and_untar alabaster-v1.0.0/zeisel-brain-sparse-results-delayed-external.tar.gz alabaster
download_and_untar alabaster-v1.0.0/zeisel-brain-sparse-results-delayed.tar.gz alabaster
download_and_untar alabaster-v1.0.0/zeisel-brain-sparse-results.tar.gz alabaster
download_and_untar alabaster-v1.0.0/zeisel-brain-dense-multimodal-results-delayed-external.tar.gz alabaster
download_and_untar alabaster-v1.0.0/zeisel-brain-dense-multimodal-results-delayed.tar.gz alabaster
download_and_untar alabaster-v1.0.0/zeisel-brain-dense-multimodal-results.tar.gz alabaster

mkdir -p ${dir}/ArtifactDB
download_and_untar ArtifactDB-v1.0.0/zeisel-brain-dense.tar.gz ArtifactDB
download_and_untar ArtifactDB-v1.0.0/zeisel-brain-sparse.tar.gz ArtifactDB
download_and_untar ArtifactDB-v1.0.0/zeisel-brain-stripped.tar.gz ArtifactDB
download_and_untar ArtifactDB-v1.0.0/zeisel-brain-sparse-results.tar.gz ArtifactDB
download_and_untar ArtifactDB-v1.0.0/zeisel-brain-dense-multimodal-results.tar.gz ArtifactDB
