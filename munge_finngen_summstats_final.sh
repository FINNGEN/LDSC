#!/bin/bash
readarray -t trait < FG_PHENOLIST.txt
for pheno in "${trait[@]}"
do
        # download summstats
    gsutil -m cp gs://fg_drugs/r5_gwas_v0/${pheno}.gz .

    gunzip ${pheno}.gz 

    head -n 1 ${pheno} > header
    
        # get sample size column for each phenotype 
    grep -w ${pheno} /mnt/work_disk/FG_UKBIO_corr/R5_drug_gwass_sample_size.txt | cut -f6 > sample_size.tmp
    awk -v n_sample="$(cat sample_size.tmp)" '{print $0"\t"n_sample}' ${pheno} | tail -n+2 > tmp && mv tmp ${pheno}
    rm sample_size.tmp
      
        # Sort and merge with finngen_rsids_info.tsv
    join -1 2 -2 1 -t $'\t' <(awk 'BEGIN{FS=OFS="\t"} {print $1"_"$2,$0}' ${pheno} | nl -nln | sort -b -k2,2) <(tail -n+2 /mnt/work_disk/FG_UKBIO_corr/finngen_rsids_info.tsv | sort -b -k1,1) | sort -k2,2g | cut -f3-15 > tmp && mv tmp ${pheno}

        # Add the INFO and RSID col names and re-attach header
    awk '{print $0, "\tn_sample", "\tINFO", "\tRSID"}' header > tmp && mv tmp header
    cat header ${pheno} > tmp && mv tmp ${pheno}
    rm header

    cut -f1,2,3,4,6,7,8,9,10,11,12,13 ${pheno} > tmp && mv tmp ${pheno}

        # Munge summary stats
    /mnt/work_disk/ldsc/munge_sumstats.py \
        --sumstats ${pheno} \
        --out sumstats/${pheno} \
        --snp RSID \
        --a2 alt \
        --a1 ref \
        --p pval \
        --N-col n_sample \
        --info INFO

    rm ${pheno}

done
