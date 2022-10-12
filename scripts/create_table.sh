#!/bin/bash

SS=$1
COUNT=$2
OUT=$3
i=0

rm -rf $OUT && touch $OUT
N=$(wc -l < $SS)

while read f
do
    i=$((i+1))
    printf "\r$i/$N"
    PHENO=$(basename $f .gz)
    echo -e "$PHENO\t$f" | join - -t $'\t' <(cut -f 1,2 $COUNT) >> $OUT
done < $SS 
