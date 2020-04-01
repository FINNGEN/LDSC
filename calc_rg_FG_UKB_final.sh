#!/bin/bash

awk '{ print $1 }' FG_UKB_PHENOLIST.txt  >  FG_UKB_PHENOLIST_column1.txt
awk '{ print $2 }' FG_UKB_PHENOLIST.txt  >  FG_UKB_PHENOLIST_column2.txt

readarray -t column1 < FG_UKB_PHENOLIST_column1.txt
readarray -t column2 < FG_UKB_PHENOLIST_column2.txt

i=0;
for trait in "${column1[@]}"
do

    #sumstats from LDSC munge pipeline (FG and UKB)
    pheno1=/mnt/work_disk/FG_UKBIO_corr/FG_R5_DRUG_GWAS/sumstats/${column1[$i]}.sumstats.gz
    pheno2=/mnt/work_disk/FG_UKBIO_corr/UKBIO_GWAS/sumstats/${column2[$i]}.sumstats.gz
    output=${column1[$i]}_${column2[$i]}

    gunzip -c ${pheno1} > pheno1_tmp
    gunzip -c ${pheno2} > pheno2_tmp

      # join summstats based on rsID
    join -1 5 -2 5 -t $'\t' <(sort -b -V -k5,5 pheno1_tmp) <(sort -b -V -k5,5 pheno2_tmp) > summstats_combined

      # combine allele columns from both summstats
    awk '{print $1"\t"$2":"$3"\t"$8":"$9}' summstats_combined > summstats_combined_tmp

      # compare allele columns and output SNPs that do not have matching alleles
    awk '$2!=$3{print $0}' summstats_combined_tmp > summstats_combined_tmp2

      # from the SNPs that don't match, flip alleles of second column
    awk -F ":" '{print $1,$2,$3}' summstats_combined_tmp2 | awk '{print $1"\t"$2":"$3"\t"$5":"$4}' > summstats_combined_tmp3

      # compare alleles again (to make sure non-matching alleles aren't due to reversal of alleles) and output list of rsIDs that have incompatible alleles
    awk '$2!=$3{print $1}' summstats_combined_tmp3 > snps_incompatible_alleles.txt

      # exclude these SNPs from munged sumstats
    awk 'FNR==NR { a[$1]; next } !($5 in a)' snps_incompatible_alleles.txt pheno1_tmp > pheno1_tmp2
    awk 'FNR==NR { a[$1]; next } !($5 in a)' snps_incompatible_alleles.txt pheno2_tmp > pheno2_tmp2

      # re-zip sumstats
    gzip pheno1_tmp2
    gzip pheno2_tmp2

    /mnt/work_disk/ldsc/ldsc.py \
        --rg pheno1_tmp2.gz,pheno2_tmp2.gz \
        --ref-ld-chr /mnt/work_disk/FG_UKBIO_corr/eur_w_ld_chr/ \
        --w-ld-chr /mnt/work_disk/FG_UKBIO_corr/eur_w_ld_chr/ \
        --out rg/${output}

            # remove intermediate files
    rm summstats_combined*
    rm snps_incompatible_alleles.txt 
    rm pheno1*
    rm pheno2*

    let "i=i+1"

done

rm FG_UKB_PHENOLIST_column*
