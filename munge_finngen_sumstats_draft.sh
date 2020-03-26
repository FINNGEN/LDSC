#!/bin/bash

pheno=GWAS_of_interest_from_FinnGen
sample_size=N
output=name_of_output

gunzip ${pheno}.gz

head -n 1 ${pheno} > header
tail -n +2 ${pheno} > tmp && mv tmp ${pheno}
  
    # Sort and merge with finngen_rsids_info.tsv
join -1 2 -2 1 -t $'\t' <(awk 'BEGIN{FS=OFS="\t"} {print $1"_"$2,$0}' ${pheno} | nl -nln | sort -b -k2,2) <(tail -n+2 /mnt/work_disk/FG_UKBIO_corr/finngen_rsids_info.tsv | sort -b -k1,1) | sort -k2,2g | cut -f3-14 > tmp && mv tmp ${pheno}

    # Add the INFO and RSID col names and re-attach header
awk '{print $0, "\tINFO", "\tRSID"}' header > tmp && mv tmp header
cat header ${pheno} > tmp && mv tmp ${pheno}
rm header

cut -f1,2,3,4,6,7,8,9,10,11,12 ${pheno} > tmp && mv tmp ${pheno}

    # Munge summary stats
/mnt/work_disk/ldsc/munge_sumstats.py \
    --sumstats ${pheno} \
    --out ${output} \
    --snp RSID \
    --a2 alt \
    --a1 ref \
    --p pval \
    --N ${sample_size} \
    --info INFO
    #--merge-alleles /mnt/work_disk/ldsc/w_hm3.snplist
