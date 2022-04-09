#!/bin/bash

set -e
set -u

base=https://github.com/jkanche/kana-formats/releases/download
dir=legacy
mkdir -p ${dir}

download() {
    output=${dir}/$2
    if [ ! -e $output ]
    then
        curl -L ${base}/$1 > $output
    fi
}

download zeisel-v0-latest/zeisel_tenx_20220307.kana zeisel_tenx_20220307.kana
download zeisel-v0-latest/zeisel_mtx_20220306.kana zeisel_mtx_20220306.kana
