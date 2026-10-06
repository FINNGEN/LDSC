#!/bin/bash
# Build the PHENO<TAB>PATH<TAB>N input table for ldsc.wdl (meta_fg / meta_other).
# Only phenos with a sumstat in the bucket AND an entry in the pheno_n file are kept.
# N = num_cases + num_controls. No header (the wdl splits the table into chunks).
#
# Usage: build_pheno_table.sh <pheno_n.tsv> <gs://bucket/path/> <out_table>
# e.g.   build_pheno_table.sh R14_pheno_n.tsv \
#          gs://finngen-production-library-green/finngen_R14/finngen_R14_analysis_data/summary_stats/release/ \
#          r14_table.txt

set -euo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 <pheno_n.tsv> <gs://bucket/path/> <out_table>" >&2
    exit 1
fi

PHENO_N=$1
BUCKET=$(echo "$2" | sed 's:/*$::')
OUT=$3

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# PHENO\tN from pheno_n file, columns located by header name
awk 'BEGIN {FS=OFS="\t"}
     NR==1 {for (i=1;i<=NF;i++) col[$i]=i
            if (!("phenocode" in col) || !("num_cases" in col) || !("num_controls" in col)) {
                print "missing phenocode/num_cases/num_controls column" > "/dev/stderr"; exit 1}
            next}
     {print $col["phenocode"], $col["num_cases"]+$col["num_controls"]}' "$PHENO_N" | sort -t$'\t' -k1,1 > "$TMP/n.txt"

# PHENO\tPATH from bucket (sumstats only, skip .tbi)
gsutil ls "$BUCKET/*.gz" | grep "\.gz$" | while read -r f; do
    printf '%s\t%s\n' "$(basename "$f" .gz)" "$f"
done | sort -t$'\t' -k1,1 > "$TMP/paths.txt"

join -t$'\t' "$TMP/paths.txt" "$TMP/n.txt" > "$OUT"

# report mismatches
n_missing_n=$(join -t$'\t' -v1 "$TMP/paths.txt" "$TMP/n.txt" | wc -l)
n_missing_ss=$(join -t$'\t' -v2 "$TMP/paths.txt" "$TMP/n.txt" | wc -l)
echo "sumstats in bucket:        $(wc -l < "$TMP/paths.txt")" >&2
echo "phenos in pheno_n:         $(wc -l < "$TMP/n.txt")" >&2
echo "written to $OUT:  $(wc -l < "$OUT")" >&2
if [ "$n_missing_n" -gt 0 ]; then
    echo "sumstats without N (dropped): $n_missing_n" >&2
    join -t$'\t' -v1 "$TMP/paths.txt" "$TMP/n.txt" | cut -f1 | head -5 | sed 's/^/  /' >&2
fi
if [ "$n_missing_ss" -gt 0 ]; then
    echo "phenos without sumstats (dropped): $n_missing_ss" >&2
    join -t$'\t' -v2 "$TMP/paths.txt" "$TMP/n.txt" | cut -f1 | head -5 | sed 's/^/  /' >&2
fi
