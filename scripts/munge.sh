#!/bin/bash

sumstats=$1
path=$(dirname $sumstats)
pheno=$(basename $sumstats .gz)

echo ${sumstats} ${path} ${pheno}

gunzip -c ${sumstats} |  head |\
    awk 'BEGIN{FS=OFS="\t"} \
    NR==1{for(i=1;i<=NF;i++) a[$i]=i; print "SNP","A1","A2","BETA","P"} \
    NR>1 {if($a["rsids"]!="") print $a["rsids"],$a["alt"],$a["ref"],$a["beta"],$a["pval"]}' | \
    awk 'BEGIN{FS=OFS="\t"} {n=split($1,a,","); for(i=1;i<=n;i++) print a[i],$0}' | \
    cut -f1,3- | gzip > ${path}/${pheno}.premunge.gz
