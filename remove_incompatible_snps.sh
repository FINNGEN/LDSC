#!/bin/bash

#sumstats from LDSC munge pipeline (FG and UKB)
pheno1=GWAS1_sumstats.gz
pheno2=GWAS2_.sumstats.gz
output_pheno1=name_of_output_GWAS1_sumstats
output_pheno2=name_of_output_GWAS2_sumstats

gunzip -c ${pheno1} > pheno1_tmp
gunzip -c ${pheno2} > pheno2_tmp

  # join summstats based on rsID
join -1 4 -2 5 -t $'\t' <(sort -b -V -k4,4 pheno1_tmp) <(sort -b -V -k5,5 pheno2_tmp) > summstats_combined

  # combine allele columns from both summstats
awk '{print $1"\t"$2":"$3"\t"$8":"$9}' summstats_combined > summstats_combined_tmp

  # compare allele columns and output SNPs that do not have matching alleles
awk '$2!=$3{print $0}' summstats_combined_tmp > summstats_combined_tmp2

  # from the SNPs that don't match, flip alleles of second column
awk -F ":" '{print $1,$2,$3}' summstats_combined_tmp2 | awk '{print $1"\t"$2":"$3"\t"$5":"$4}' > summstats_combined_tmp3

  # compare alleles again (to make sure non-matching alleles aren't due to reversal of alleles) and output list of rsIDs that have incompatible alleles
awk '$2!=$3{print $1}' summstats_combined_tmp3 > snps_incompatible_alleles.txt

  # exclude these SNPs from munged sumstats
awk 'FNR==NR { a[$1]; next } !($4 in a)' snps_incompatible_alleles.txt pheno1_tmp > ${output_pheno1}
awk 'FNR==NR { a[$1]; next } !($5 in a)' snps_incompatible_alleles.txt pheno2_tmp > ${output_pheno2}

  # re-zip sumstats
gzip ${output_pheno1}
gzip ${output_pheno2}
