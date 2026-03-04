#!/bin/bash
usage() {
    cat <<EOF
Usage:
  $0 --input <sumstats.gz> --outdir <output_dir> --beta-col <BETA_COL> --p-col <P_COL> --a1-col <A1_COL> --a2-col <A2_COL> [--chrom-col <CHR_COL> --pos-col <POS_COL> --rsid-map <rsid_map.tsv>] [--rsid <RSID_COL>] [--convert-script <path>]

  If --rsid is provided, RSID + alleles/beta/p columns are extracted and RSID mapping is skipped. 
  If --rsid is NOT provided, --chrom-col, --pos-col, --a1-col, --a2-col, --rsid-map are all required, and mapping is performed.

  --convert-script <path>   Optional: Override path to convert_rsids.py (default: ~/Dropbox/Projects/commons/rsid_map/scripts/convert_rsids.py)

  --chrom-col/--pos-col can always be provided, but are ignored if --rsid is present.
EOF
    exit 1
}

CONVERT_SCRIPT_DEFAULT=~/Dropbox/Projects/LDSC/rsid_map/scripts/convert_rsids.py
CONVERT_SCRIPT="$CONVERT_SCRIPT_DEFAULT"

INPUT="" ; OUTDIR="" ; BETA_COL="" ; P_COL="" ; A1_COL="" ; A2_COL="" ; RSID_COL="" ; CHROM_COL="" ; POS_COL="" ; RSID_MAP=""

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --input) INPUT="$2"; shift; shift ;;
        --outdir) OUTDIR="$2"; shift; shift ;;
        --beta-col) BETA_COL="$2"; shift; shift ;;
        --p-col) P_COL="$2"; shift; shift ;;
        --a1-col) A1_COL="$2"; shift; shift ;;
        --a2-col) A2_COL="$2"; shift; shift ;;
        --rsid) RSID_COL="$2"; shift; shift ;;
        --chrom-col) CHROM_COL="$2"; shift; shift ;;
        --pos-col) POS_COL="$2"; shift; shift ;;
        --rsid-map) RSID_MAP="$2"; shift; shift ;;
        --convert-script) CONVERT_SCRIPT="$2"; shift; shift ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT" || -z "$OUTDIR" || -z "$BETA_COL" || -z "$P_COL" || -z "$A1_COL" || -z "$A2_COL" ]]; then
    usage
fi

mkdir -p "$OUTDIR"
BASE=$(basename "$INPUT" .txt.gz)
OUTFILE="${OUTDIR}/${BASE}.premunge.gz"

if [[ -n "$RSID_COL" ]]; then
    python3 -c "import gzip,sys,csv;f,o,rs,a1,a2,b,p=sys.argv[1:8];r=csv.DictReader(gzip.open(f,'rt') if f.endswith('.gz') else open(f),delimiter='\t');w=csv.writer(open(o,'w'),delimiter='\t');w.writerow(['SNP','A1','A2','BETA','P']);[w.writerow([row[rs],row[a1],row[a2],row[b],row[p]]) for row in r]" \
        "$INPUT" "${OUTDIR}/${BASE}.premunge" "$RSID_COL" "$A1_COL" "$A2_COL" "$BETA_COL" "$P_COL"
    gzip -f "${OUTDIR}/${BASE}.premunge"
    echo "All done! Final file is: ${OUTFILE}"
    exit 0
fi

if [[ -z "$CHROM_COL" || -z "$POS_COL" || -z "$RSID_MAP" ]]; then
    echo "Error: --chrom-col, --pos-col, --rsid-map are required when --rsid is not specified."
    usage
fi

python3 -c "import gzip,sys,csv;f,o,c,p,a1,a2,b,pc=sys.argv[1:9];r=csv.DictReader(gzip.open(f,'rt') if f.endswith('.gz') else open(f),delimiter='\t');w=csv.writer(open(o,'w'),delimiter='\t');w.writerow(['CHR','POS','REF','ALT','BETA','P']);[w.writerow([row[c],row[p],row[a1],row[a2],row[b],row[pc]]) for row in r]" \
    "$INPUT" "${OUTDIR}/${BASE}.for_rsid.tsv" "$CHROM_COL" "$POS_COL" "$A1_COL" "$A2_COL" "$BETA_COL" "$P_COL"

awk -F"\t" 'BEGIN{OFS="\t"} NR==1{print "ID","A1","A2","BETA","P"} NR>1{gsub(/^chr/,"",$1); print $1 ":" $2 ":" $3 ":" $4, $3, $4, $5, $6}' \
  "${OUTDIR}/${BASE}.for_rsid.tsv" > "${OUTDIR}/${BASE}.pre_convert.tsv"

python3 "$CONVERT_SCRIPT" \
    --file "${OUTDIR}/${BASE}.pre_convert.tsv" \
    --map "$RSID_MAP" \
    --out "$OUTDIR" \
    --to-rsid \
    --metadata 0 1 2 \
    --columns 0 1 2 3 4

awk 'NR==1{print "SNP\tA1\tA2\tBETA\tP"} NR>1' "${OUTDIR}/${BASE}.pre_convert.rsid" > "${OUTDIR}/${BASE}.premunge"
gzip -f "${OUTDIR}/${BASE}.premunge"

echo "All done! Final file is: ${OUTFILE}"
