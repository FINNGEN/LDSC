#!/bin/bash

pheno=GWAS_of_interest

gunzip ${pheno}.gz

head -n 1 ${pheno} > header
tail -n +2 ${pheno} > tmp && mv tmp ${pheno}
  
    # Sort and merge with finngen_rsids.tsv
join -1 2 -2 2 -t $'\t' <(awk 'BEGIN{FS=OFS="\t"} {print $1"_"$2,$0}' ${pheno} | nl -nln | sort -b -k2,2) <(sort -b -k2,2 /mnt/work_disk/FG_UKBIO_corr/finngen_rsids.tsv) | sort -k2,2g | cut -f3-13 > tmp && mv tmp ${pheno}


    # Add the RSID field and re-attach header
awk '{print $0, "\tRSID"}' header > tmp && mv tmp header
cat header ${pheno} > tmp && mv tmp ${pheno}
rm header

cut -f1,2,3,4,6,7,8,9,10,11 ${pheno} > tmp && mv tmp ${pheno}

# Munge summary stats
/mnt/work_disk/ldsc/munge_sumstats.py \
    --sumstats ${pheno} \
    --out /mnt/work_disk/FG_UKBIO_corr/r5_J01 \
    --snp RSID \
    --a2 alt \
    --a1 ref \
    --p pval \
    --N 188038 
    #--merge-alleles /mnt/work_disk/ldsc/w_hm3.snplist
