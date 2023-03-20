#!/bin/bash

set -e
set -u
dir=feature-sets
mkdir -p ${dir}

base=https://github.com/LTLA/gesel-feedstock/releases/download/indices-v0.2.1
download() {
    output=${dir}/$1
    if [ ! -e $output ]
    then
        curl -L ${base}/$1 > $output
    fi
}

download 10090_collections.tsv.gz
download 10090_gene2set.tsv.gz
download 10090_set2gene.tsv.gz
download 10090_sets.tsv.gz

download 9606_collections.tsv.gz
download 9606_gene2set.tsv.gz
download 9606_set2gene.tsv.gz
download 9606_sets.tsv.gz

gbase=https://github.com/LTLA/gesel-feedstock/releases/download/genes-v1.0.0
gdownload() {
    output=${dir}/$1
    if [ ! -e $output ]
    then
        curl -L ${gbase}/$1 > $output
    fi
}

gdownload 10090_ensembl.tsv.gz
gdownload 10090_symbol.tsv.gz
gdownload 10090_entrez.tsv.gz

gdownload 9606_ensembl.tsv.gz
gdownload 9606_symbol.tsv.gz
gdownload 9606_entrez.tsv.gz
