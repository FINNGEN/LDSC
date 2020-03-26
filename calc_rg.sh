#!/bin/bash

pheno1=GWAS1_sumstats.gz
pheno2=GWAS2_sumstats.gz
output=name_of_output


/mnt/work_disk/ldsc/ldsc.py \
        --rg ${pheno1},${pheno2} \
        --ref-ld-chr eur_w_ld_chr/ \
        --w-ld-chr eur_w_ld_chr/ \
        --out ${output}
