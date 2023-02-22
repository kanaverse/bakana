#!/bin/bash

set -e
set -u

base=https://github.com/LTLA/singlepp-references/releases/download/v2.0.0
dir=references
mkdir -p ${dir}

download() {
    output=${dir}/$2
    if [ ! -e $output ]
    then
        curl -L ${base}/$1 > $output
    fi
}

download BlueprintEncode_genes.csv.gz BlueprintEncode_genes.csv.gz
download BlueprintEncode_matrix.csv.gz BlueprintEncode_matrix.csv.gz
download BlueprintEncode_label_names_fine.csv.gz BlueprintEncode_label_names_fine.csv.gz
download BlueprintEncode_labels_fine.csv.gz BlueprintEncode_labels_fine.csv.gz
download BlueprintEncode_markers_fine.gmt.gz BlueprintEncode_markers_fine.gmt.gz

download ImmGen_genes.csv.gz ImmGen_genes.csv.gz
download ImmGen_matrix.csv.gz ImmGen_matrix.csv.gz
download ImmGen_label_names_fine.csv.gz ImmGen_label_names_fine.csv.gz
download ImmGen_labels_fine.csv.gz ImmGen_labels_fine.csv.gz
download ImmGen_markers_fine.gmt.gz ImmGen_markers_fine.gmt.gz

download MouseRNAseq_genes.csv.gz MouseRNAseq_genes.csv.gz
download MouseRNAseq_matrix.csv.gz MouseRNAseq_matrix.csv.gz
download MouseRNAseq_label_names_fine.csv.gz MouseRNAseq_label_names_fine.csv.gz
download MouseRNAseq_labels_fine.csv.gz MouseRNAseq_labels_fine.csv.gz
download MouseRNAseq_markers_fine.gmt.gz MouseRNAseq_markers_fine.gmt.gz
