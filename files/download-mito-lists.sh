#!/bin/bash

set -e
set -u

base=https://github.com/kanaverse/kana-special-features/releases/download
dir=mito-lists
mkdir -p ${dir}

download() {
    output=${dir}/$2
    if [ ! -e $output ]
    then
        curl -L ${base}/$1 > $output
    fi
}

download v1.0.0/9606-mito-ensembl.txt.gz 9606-mito-ensembl.txt.gz 
download v1.0.0/10090-mito-ensembl.txt.gz 10090-mito-ensembl.txt.gz 
download v1.0.0/9606-mito-symbol.txt.gz 9606-mito-symbol.txt.gz 
download v1.0.0/10090-mito-symbol.txt.gz 10090-mito-symbol.txt.gz 
