#!/bin/bash

DESTINATION_BUCKET="gs://finngen-production-library-green/ldsc/r13/premunge"
SSPATH="gs://r13-data/regenie/release/summary_stats_annotated/"
ID='b9606267-efd5-4f4b-bf91-c78d084d08de'
JSON_FILE="~/Dropbox/Projects/CromwellInteract/tmp/${ID}.json"
OUTPUT_FILE="r13_table.txt"
N_COUNT_FILE="n_all.txt"
PREMUNGE_TABLE_FILE="r13_premunge_table.txt"

# Default values for test mode
TEST_MODE="false"
MAX_PATHS_TO_PROCESS=0 # 0 means no limit by default

# Parse command-line arguments
while getopts "tn:" opt; do
    case ${opt} in
        t ) # Enable test mode
            TEST_MODE="true"
            ;;
        n ) # Set maximum paths to process
            MAX_PATHS_TO_PROCESS=$OPTARG
            ;;
        \? )
            echo "Usage: $0 [-t] [-n <number_of_paths>]"
            echo "  -t: Enable test mode (process only a limited number of paths)"
            echo "  -n <number>: Specify the maximum number of paths to process in test mode"
            exit 1
            ;;
    esac
done
shift $((OPTIND -1))

# Inform user about test mode status
if [ "$TEST_MODE" = "true" ] && [ "$MAX_PATHS_TO_PROCESS" -eq 0 ]; then
    MAX_PATHS_TO_PROCESS=10
    echo "Test mode enabled with default limit of 10 paths."
elif [ "$TEST_MODE" = "true" ]; then
    echo "Test mode enabled with limit of $MAX_PATHS_TO_PROCESS paths."
else
    echo "Test mode is disabled. All paths will be processed."
fi

# Generate n_all.txt if it does not exist
if [ ! -f "$N_COUNT_FILE" ]; then
    echo "Generating $N_COUNT_FILE..."
    cat ~/r13/pheno/R13_pheno_n.tsv | awk 'BEGIN {FS="\t"; OFS="\t"} NR > 1 {print $1, $3+$4}' > "$N_COUNT_FILE"
else
    echo "$N_COUNT_FILE already exists, skipping generation."
fi

# Generate r13_premunge_table.txt if it does not exist
if [ ! -f "$PREMUNGE_TABLE_FILE" ]; then
    echo "Generating $PREMUNGE_TABLE_FILE..."
    join -t $'\t' <(while read f; do base=$(basename "$f" .gz) && echo -e "$base\t$f" ; done < <(gsutil ls "$SSPATH"*gz )) "$N_COUNT_FILE" > "$PREMUNGE_TABLE_FILE"
else
    echo "$PREMUNGE_TABLE_FILE already exists, skipping generation."
fi

declare -A n_counts

while IFS=$'\t' read -r pheno_id n_val; do
    n_counts["$pheno_id"]="$n_val"
done < "$N_COUNT_FILE"

# Temporary file to store paths for batch copy and output generation
TEMP_PATHS_FILE=$(mktemp)
processed_count=0

# First pass: Collect paths and prepare data for batch copy and output
while read -r original_gcs_path; do
    if [ "$TEST_MODE" = "true" ] && [ "$processed_count" -ge "$MAX_PATHS_TO_PROCESS" ]; then
        echo "Test mode limit reached during path collection. Processing up to $MAX_PATHS_TO_PROCESS paths."
        break
    fi

    if [ -z "$original_gcs_path" ]; then
        continue
    fi

    full_filename=$(basename "$original_gcs_path")
    filename_base=$(echo "$full_filename" | sed 's/\.premunge\.gz$//')
    new_gcs_path="${DESTINATION_BUCKET}/${full_filename}"
    n_count_val=${n_counts["$filename_base"]}
    if [ -z "$n_count_val" ]; then
        n_count_val="NA"
    fi

    # Store relevant data for later processing: original_path, base_name, new_path, n_count
    echo -e "${original_gcs_path}\t${filename_base}\t${new_gcs_path}\t${n_count_val}" >> "$TEMP_PATHS_FILE"
    processed_count=$((processed_count + 1))
done < <(eval "jq -r '.outputs[\"munge_fg.paths\"][]' ${JSON_FILE}")

# Extract only original GCS paths from the temporary file for batch copy
echo "Starting batch copy to $DESTINATION_BUCKET..."
if [ -s "$TEMP_PATHS_FILE" ]; then # Check if temp file is not empty
    awk '{print $1}' "$TEMP_PATHS_FILE" | gsutil -m cp -I "$DESTINATION_BUCKET"
else
    echo "No paths to copy."
fi

# Generate the final OUTPUT_FILE from the collected data using awk for efficiency
echo "Generating $OUTPUT_FILE..."
awk -F'\t' '{print $2"\t"$3"\t"$4}' "$TEMP_PATHS_FILE" > "$OUTPUT_FILE"

# Clean up temporary file
rm "$TEMP_PATHS_FILE"
echo "Script finished. Output in $OUTPUT_FILE"
