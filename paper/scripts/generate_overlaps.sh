#!/usr/bin/env bash
set -euo pipefail

n_reads=5000
# List of input files (adjust paths as necessary)
input_files=($(fd -e fq.gz -a . "../data/reads/simplex/v4.3.0/dna_r10.4.1_e8.2_400bps_sup@v4.3.0"))

# Loop over each file in the input list
for input_file in "${input_files[@]}"; do
    # Extract base filename without extension
    base_name=$(basename "$input_file" .fq.gz)
    
    # Split the base name on double underscore and take the left part
    base_name_left=$(echo "$base_name" | awk -F'__' '{print $1}')
    
    if [ "$base_name_left" == "ATCC_14035" ]; then
        continue
    fi

    echo "Processing $base_name"

    fq="${base_name_left}.${n_reads}.fq"
    # Run rasusa to downsample reads
    rasusa reads -n 5000 -o "$fq" "$input_file"
    
    paf="${base_name_left}.5000.paf"
    # Run minimap2 to perform read overlap
    minimap2 -t 8 -x ava-ont "$fq" "$fq" > "$paf"
    
    # Calculate GC content and output to a TSV file
    seqkit fx2tab --name --only-id --gc "$fq" > "${base_name_left}.gc.tsv"
    
    # Extract genome size from reference file stats
    seqkit stats -T "../data/references/${base_name}.fa" | cut -f5 | tail -1 > "${base_name_left}.gsize"
done

echo "Processing complete for all files."
