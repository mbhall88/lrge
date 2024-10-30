#!/usr/bin/env bash
set -euo pipefail

exec 2> "${snakemake_log[0]}"

tmpdir=$(mktemp -d)
trap 'rm -rf $tmpdir' EXIT

run=${snakemake_wildcards[run]}

kingfisher get ${snakemake_params[opts]} --output-directory "$tmpdir" --ascp-ssh-key "${snakemake_params[ascp_ssh_key]}"

# Find all files matching the pattern <run>*.fastq.gz in the output directory
matches=$(find "$tmpdir" -type f -name "*.fastq.gz")

# Count how many files match the pattern
count=$(echo "$matches" | wc -l)

full_fastq="${tmpdir}/${run}.all.fq.gz"

if [ $count -eq 1 ]; then
    mv "$matches" "$full_fastq"
    echo "[SUCCESS]: Renamed $matches to $full_fastq" >&2
elif [ $count -gt 1 ]; then
    # If more than one match, concatenate them
    cat $matches > "$full_fastq"
    
    echo "[SUCCESS]: Concatenated the following files into $full_fastq:" >&2
    echo "$matches" >&2
else
    # If no matches found, print an error and exit
    echo "[ERROR]: No matching files found for $run" >&2
    ls -la "$tmpdir" >&2
    exit 1
fi

# downsample to the maximum number of bases
max_bases=${snakemake_params[max_bases]}
seed=${snakemake_params[seed]}
output="${snakemake_output[fastq]}"

echo "[INFO]: Downsampling to $max_bases bases" >&2
rasusa reads -s $seed -b $max_bases -o "$output" "$full_fastq"

# get stats for the fastq file
seqkit stats -aT "$output" > "${snakemake_output[stats]}"