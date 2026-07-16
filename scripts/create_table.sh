#!/bin/bash

BUCKET=$1
COUNTS=$2
OUT=$3

FIX_BUCKET=$(echo $BUCKET | sed 's:/*$::')
rm -f meta_paths.txt
while read f; do tmp=$(basename $f .premunge.gz) && b=$(basename $tmp .gz) && echo -e $b"\t"$f  >> meta_paths.txt; done < <(gsutil ls $FIX_BUCKET/**/*gz )
join -t $'\t' <(sort -k1 meta_paths.txt) $COUNTS > $OUT
rm -f meta_paths.txt

