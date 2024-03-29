#!/bin/bash

set -e
set -u

mode=$(echo $1 | sed "s/\/$//")
if [ $mode != "browser" ] && [ $mode != "main" ]
then
    echo "need to specify 'module' or 'main' as the first argument"
    exit 1
fi

if [ $mode == "main" ]
then
    toss=web
    keep=node
else
    toss=node
    keep=web
fi

rm -rf ${mode}
mkdir -p ${mode}
cp -r src/* ${mode}

version=$(npm pkg get version)
echo "export const bakana_version=${version};" > ${mode}/version.js

for abdirs in readers/abstract readers/utils/abstract steps/abstract steps/utils/abstract dump/abstract
do 
    rm ${mode}/${abdirs}/*_${toss}.js

    to_rename=$(ls ${mode}/${abdirs}/*_${keep}.js)
    for x in ${to_rename[@]}
    do
        newname=$(echo $x | sed "s/_${keep}\\.js$/.js/g")
        mv $x $newname
    done
done
