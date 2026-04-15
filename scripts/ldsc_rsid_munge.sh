#!/bin/bash

# Usage Info
usage() {
    cat <<EOF
Usage:
  $0 --input <sumstats.gz> --outdir <output_dir> --beta-col <BETA_COL> --p-col <P_COL> --a1-col <A1_COL> --a2-col <A2_COL>
       [--chrom-col <CHR_COL> --pos-col <POS_COL> --rsid-map <rsid_map.tsv>] [--rsid <RSID_COL>] [--convert-script <path>]

- If --rsid is provided, it uses the input RSID column directly, skipping mapping.
- If --rsid is NOT provided, it requires chrom/pos/alleles and maps to RSID using --rsid-map.

Other options:
  --convert-script <path>   (default: ~/Dropbox/Projects/commons/rsid_map/scripts/convert_rsids.py)
EOF
    exit 1
}

CONVERT_SCRIPT=${CONVERT_SCRIPT:-~/Dropbox/Projects/LDSC/rsid_map/scripts/convert_rsids.py}

# Parse args
while [[ $# -gt 0 ]]; do
    key="$1"; shift
    case $key in
        --input) INPUT="$1"; shift ;;
        --outdir) OUTDIR="$1"; shift ;;
        --beta-col) BETA_COL="$1"; shift ;;
        --p-col)   P_COL="$1"; shift ;;
        --a1-col)  A1_COL="$1"; shift ;;
        --a2-col)  A2_COL="$1"; shift ;;
        --rsid)    RSID_COL="$1"; shift ;;
        --chrom-col) CHROM_COL="$1"; shift ;;
        --pos-col)   POS_COL="$1"; shift ;;
        --rsid-map)  RSID_MAP="$1"; shift ;;
        --convert-script) CONVERT_SCRIPT="$1"; shift ;;
        *) usage ;;
    esac
done

# Required checks
if [[ -z "$INPUT" || -z "$OUTDIR" || -z "$BETA_COL" || -z "$P_COL" || -z "$A1_COL" || -z "$A2_COL" ]]; then usage; fi
mkdir -p "$OUTDIR"

# Basename stripping .gz and .txt
BASE=$(basename "$INPUT")
for ext in .txt.gz .tsv.gz .gz .txt .tsv; do
    BASE="${BASE%$ext}"
done
OUTFILE="${OUTDIR}/${BASE}.premunge.gz"
echo "BASE: $BASE"
echo "OUTFILE: $OUTFILE"

##################################
# If --rsid is given: Shortcut: #
##################################
if [[ -n "$RSID_COL" ]]; then
    # Directly write compressed premunge.gz with the correct header from the required columns
    python3 -c "
import gzip, sys, csv
f, o, rs, a1, a2, b, p = sys.argv[1:8]
with (gzip.open(f, 'rt') if f.endswith('.gz') else open(f)) as fin, gzip.open(o, 'wt') as fout:
    reader = csv.DictReader(fin, delimiter='\t')
    writer = csv.writer(fout, delimiter='\t', lineterminator='\n')
    writer.writerow(['SNP', 'A1', 'A2', 'BETA', 'P'])  # Output header
    for row in reader:
        writer.writerow([row[rs], row[a1], row[a2], row[b], row[p]])
" \
    "$INPUT" "$OUTFILE" "$RSID_COL" "$A1_COL" "$A2_COL" "$BETA_COL" "$P_COL"
    echo "All done! Final file is: ${OUTFILE}"
    exit 0
fi

#####################################################
# Otherwise: Build SNPIDs from fields and map to RSID #
#####################################################
if [[ -z "$CHROM_COL" || -z "$POS_COL" || -z "$RSID_MAP" ]]; then
    echo "ERROR: When --rsid is not used, you must supply --chrom-col, --pos-col, and --rsid-map."
    usage
fi

TMP_PRECONVERT="${OUTDIR}/${BASE}.pre_convert.tsv.gz"
TMP_RSID="${OUTDIR}/${BASE}.pre_convert.tsv.rsid.gz"

# 1. Build "SNPID" column as "CHR:POS:A1:A2"
python3 -c "
import gzip, sys, csv
infile, outfile, CHR, POS, A1, A2, BETA, P = sys.argv[1:9]
with (gzip.open(infile, 'rt') if infile.endswith('.gz') else open(infile)) as fin, gzip.open(outfile, 'wt') as fout:
    reader = csv.DictReader(fin, delimiter='\t')
    writer = csv.writer(fout, delimiter='\t', lineterminator='\n')
    writer.writerow(['SNP','A1','A2','BETA','P'])
    for row in reader:
        snpid = f'{row[CHR]}:{row[POS]}:{row[A1]}:{row[A2]}'
        writer.writerow([snpid, row[A1], row[A2], row[BETA], row[P]])
" \
    "$INPUT" "$TMP_PRECONVERT" "$CHROM_COL" "$POS_COL" "$A1_COL" "$A2_COL" "$BETA_COL" "$P_COL"

# 2. Map SNPID to RSIDs using convert_rsids.py (produces .rsid.gz)
python3 "$CONVERT_SCRIPT" \
    --file "$TMP_PRECONVERT" \
    --map "$RSID_MAP" \
    --out "$OUTDIR" \
    --to-rsid --gz \
    --metadata 0 \
    --columns 0 1 2 3 4

# 3. Finalize: Move result to .premunge.gz
mv "$TMP_RSID" "$OUTFILE"
echo "All done! Final file is: $OUTFILE"


#########################################################
# FILES CREATED (Example: PHENO.gz) and What They Are:  #
#########################################################
# --- If called with --input PHENO.gz --outdir OUT/ ...:
#  
# 1. OUT/PHENO.pre_convert.tsv.gz
#     - Contains: columns (SNPID, A1, A2, BETA, P)
#     - SNPID is CHR:POS:A1:A2 made from your requested columns.
#
# 2. OUT/PHENO.pre_convert.tsv.rsid.gz
#     - Like above, but with SNPID mapped to RSID wherever possible.
#
# 3. OUT/PHENO.premunge.gz
#     - Final output (identical to above, with correct naming)
#
#   (Intermediate files can be deleted if desired.)
