#!/bin/bash


#build table with N
cat ~/r13/pheno/R13_pheno_n.tsv | awk 'BEGIN {FS="\t"; OFS="\t"} NR > 1 {print $1, $3+$4}'  > n_all.txt
#PREMUNGE TABLE
join -t $'\t' <(while read f; do base=$(basename $f .gz) && echo -e $base"\t"$f ; done <  <(gsutil ls gs://r13-data/regenie/release/summary_stats/*gz  )) n_all.txt  > r13_premunge_table.txt


# get paths from json and do same as above
ID='d065a7c6-0880-4550-965c-0def5cb40dcd'
join -t $'\t' <( while read f; do base=$(basename $f .premunge.gz) && echo -e $base"\t"$f  ; done < <( jq -r '.outputs["munge_fg.paths"][]' ~/Dropbox/Projects/CromwellInteract/tmp/$ID.json ) | sort)  n_all.txt  > r13_meta.txt 
