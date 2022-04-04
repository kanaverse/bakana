#!/bin/bash

set -e
set -u

base=https://github.com/clusterfork/singlepp-references/releases/download
dir=references
mkdir -p ${dir}

download() {
    output=${dir}/$2
    if [ ! -e $output ]
    then
        curl -L ${base}/$1 > $output
    fi
}

download hs-latest/BlueprintEncode_genes.csv.gz BlueprintEncode_genes.csv.gz
download hs-latest/BlueprintEncode_matrix.csv.gz BlueprintEncode_matrix.csv.gz
download hs-latest/BlueprintEncode_label_names_fine.csv.gz BlueprintEncode_label_names_fine.csv.gz
download hs-latest/BlueprintEncode_labels_fine.csv.gz BlueprintEncode_labels_fine.csv.gz
download hs-latest/BlueprintEncode_markers_fine.gmt.gz BlueprintEncode_markers_fine.gmt.gz

download mm-latest/ImmGen_genes.csv.gz ImmGen_genes.csv.gz
download mm-latest/ImmGen_matrix.csv.gz ImmGen_matrix.csv.gz
download mm-latest/ImmGen_label_names_fine.csv.gz ImmGen_label_names_fine.csv.gz
download mm-latest/ImmGen_labels_fine.csv.gz ImmGen_labels_fine.csv.gz
download mm-latest/ImmGen_markers_fine.gmt.gz ImmGen_markers_fine.gmt.gz
