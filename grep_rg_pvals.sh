#!/bin/bash

readarray -t output < RG_FG_UKB_OUTPUT.txt
for pheno in "${output[@]}"
do
                # for each log file 
                        grep '^pheno1_tmp2.gz' "$pheno".log | awk '{print $3, $6}' >> rg_pvals.txt

                done
