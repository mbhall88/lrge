#!/usr/bin/env bash
set -euxo pipefail

exec 2> "${snakemake_log[0]}"

tmpdir=$(mktemp -d)
trap 'rm -rf $tmpdir' EXIT

run=${snakemake_wildcards[run]}

kingfisher get -m ena-ascp ena-ftp --check-md5sums -f fastq.gz --force -r $run --debug --output-directory "$tmpdir" --ascp-ssh-key "${snakemake_params[ascp_ssh_key]}"

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

# download the RefSeq assembly for the run
accession=${snakemake_params[asm_accession]}
zipfile="${tmpdir}/dataset.zip"
datasets download genome accession $accession --filename "$zipfile"
asm="${tmpdir}/${accession}.fna"
unzip -p "$zipfile" '*genomic.fna' > "$asm"

# map the reads to the assembly and keep only those that map
clean_fastq="${tmpdir}/${run}.clean.fq"
minimap2 -x ${snakemake_params[platform]} -t ${snakemake[threads]} "$asm" "$full_fastq" | 
    cut -f 1 | 
    sort -u | 
    seqkit grep -f - -o "$clean_fastq" "$full_fastq"

rm "$full_fastq"

# downsample to the maximum number of bases
max_bases=${snakemake_params[max_bases]}
min_bases=${snakemake_params[min_bases]}
seed=${snakemake_params[seed]}
output="${snakemake_output[fastq]}"

echo "[INFO]: Downsampling to $max_bases bases" >&2
rasusa reads -s $seed -b $max_bases -O u "$clean_fastq" | seqkit rename -o "$output" 

# get stats for the fastq file
tmpstats=$(mktemp)
tmpinfo=$(mktemp)
seqkit stats -aT "$output" > "$tmpstats"

# get the 5th column from the second row of the stats file, which is the number of bases
num_bases=$(awk 'NR==2 {print $5}' "$tmpstats")

# if the number of bases is less than the minimum, exit with an error
if [ $num_bases -lt $min_bases ]; then
    echo "[ERROR]: Cleaning and downsampling $run resulted in $num_bases bases, which is less than the minimum of $min_bases" >&2
    exit 1
fi

# get read length and quality for each read
seqkit fx2tab -nilqH "$output" > "$tmpinfo"

mv "$tmpstats" "${snakemake_output[stats]}"
mv "$tmpinfo" "${snakemake_output[info]}"
