#!/bin/bash

set -e
set -u

base=https://github.com/kanaverse/kana-feature-sets/releases/download
dir=feature-sets
mkdir -p ${dir}

download() {
    output=${dir}/$2
    if [ ! -e $output ]
    then
        curl -L ${base}/$1 > $output
    fi
}

download mouse-GO_3.16.0/mouse-GO_3.16.0_features.csv.gz mouse-GO_3.16.0_features.csv.gz 
download mouse-GO_3.16.0/mouse-GO_3.16.0_sets.txt.gz mouse-GO_3.16.0_sets.txt.gz 

download human-GO_3.16.0/human-GO_3.16.0_features.csv.gz human-GO_3.16.0_features.csv.gz 
download human-GO_3.16.0/human-GO_3.16.0_sets.txt.gz human-GO_3.16.0_sets.txt.gz 

