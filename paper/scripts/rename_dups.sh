#!/usr/bin/env bash
set -euo pipefail

# this script take a list of accession IDs and a search directory
# it finds the fastq file for that accession in the search directory
# it then runs seqkit rename to rename reads with duplicated IDs and then overwrites the original fastq file

acc_list=$1
search_dir=$2

while read -r acc; do
    echo "[INFO]: Processing $acc" >&2
    # find the fastq file for the accession
    fastq=$(fd -e '.fastq.gz' $acc $search_dir)
    if [ -z "$fastq" ]; then
        echo "[ERROR]: No fastq file found for $acc" >&2
        exit 1
    fi
    # make sure there is only one fastq file
    if [ $(echo "$fastq" | wc -l) -ne 1 ]; then
        echo "[ERROR]: More than one fastq file found for $acc" >&2
        exit 1
    fi

    echo "[INFO]: Found $fastq" >&2

    # rename reads with duplicated IDs
    tmp_fastq=$(mktemp --suffix=.fastq.gz)

    echo "[INFO]: Renaming reads in $fastq" >&2

    seqkit rename -o $tmp_fastq $fastq
    mv $tmp_fastq $fastq
    echo "[SUCCESS]: Renamed reads in $fastq" >&2
done <$acc_list
