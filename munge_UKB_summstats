#!/bin/bash

pheno=GWAS_of_interest_from_UKB
output=name_of_output

gunzip -c ${pheno}.bgz > ${pheno}

tail -n +2 ${pheno} > tmp && mv tmp ${pheno}

    # Sort and merge with annotated variants file
join -1 2 -2 1 -t $'\t' <(cut -f1,6,9,10,12 ${pheno} | nl -nln | sort -b -V -k2,2) <(cut -f1-6,10 /mnt/work_disk/FG_UKBIO_corr/UKBIO_GWAS/variants.tsv | sort -b -V -k1,1) | sort -k2,2g > tmp && mv tmp ${pheno}

    # Remove variant ID column
cut -f3-12 ${pheno} > tmp && mv tmp ${pheno}

    # Add header
echo -e "n_complete_samples\tbeta\tse\tpval\tchr\tpos\tref\talt\trsid\tinfo" > header
cat header ${pheno} > tmp && mv tmp ${pheno}
rm header

    # Munge summary stats
/mnt/work_disk/ldsc/munge_sumstats.py \
        --sumstats ${pheno} \
        --out ${output} \
        --snp rsid \
        --a2 alt \
        --a1 ref \
        --p pval \
        --N-col n_complete_samples \
        --info info
        #--merge-alleles /mnt/work_disk/ldsc/w_hm3.snplist
