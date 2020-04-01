#!/bin/bash

readarray -t trait < UKB_PHENOLIST.txt
for pheno in "${trait[@]}"
do

    gunzip -c ${pheno}.bgz > ${pheno}

    #if header has either expected_case_minor_AC or expected_min_category_minor_AC --> remove this optional column
    head -n 1 ${pheno} > header
    if egrep -w "expected_case_minor_AC|expected_min_category_minor_AC" header
        
        then 
            cut -f1-3,5-12 ${pheno} | tail -n +2 > tmp && mv tmp ${pheno}

        else
            tail -n +2 ${pheno} > tmp && mv tmp ${pheno}
    fi

        # Sort and merge with annotated variants file
    join -1 2 -2 1 -t $'\t' <(cut -f1,5,8,9,11 ${pheno} | nl -nln | sort -b -V -k2,2) <(cut -f1-6,10 /mnt/work_disk/FG_UKBIO_corr/UKBIO_GWAS/variants.tsv | sort -b -V -k1,1) | sort -k2,2g > tmp && mv tmp ${pheno}

        # Remove variant ID column
    cut -f3-12 ${pheno} > tmp && mv tmp ${pheno}

        # Add header
    echo -e "n_complete_samples\tbeta\tse\tpval\tchr\tpos\tref\talt\trsid\tinfo" > header
    cat header ${pheno} > tmp && mv tmp ${pheno}
    rm header

        # Munge summary stats
    /mnt/work_disk/ldsc/munge_sumstats.py \
            --sumstats ${pheno} \
            --out sumstats/${pheno} \
            --snp rsid \
            --a2 alt \
            --a1 ref \
            --p pval \
            --N-col n_complete_samples \
            --info info

done
