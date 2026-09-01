#!/usr/bin/env bash
# Issue #31: does the corrected internal-match filter recover the other two X. oryzae runs?
# Prep mirrors paper/workflow/scripts/download.sh; estimates use a binary built from a tree
# containing ONLY the #31 change.
set -euxo pipefail

run=$1
asm_accession=$2
threads=${3:-8}

WORKDIR=/scratch/user/uqmhal11/lrge-issue29
LRGE=/scratch/temp/27749624/claude-1100000562/-scratch-user-uqmhal11-lrge/01a9cc44-050e-460c-9be2-8cb707278dae/scratchpad/lrge-31/target/release/lrge
ASCP_KEY="$HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh"
DL_SEED=324; MAX_BASES=1G; MIN_BASES=27000000
LRGE_SEED=4556; LRGE_TARGET=10000; LRGE_QUERY=5000

outdir="$WORKDIR/$run"; mkdir -p "$outdir"; cd "$outdir"

source /home/uqmhal11/sw/miniforge3/etc/profile.d/conda.sh
conda activate kingfisher
export PATH="$PATH:/home/uqmhal11/sw/miniforge3/envs/binf/bin:/home/uqmhal11/sw/miniforge3/envs/Cstriatum/bin"

tmpdir=$(mktemp -d -p "$outdir"); trap 'rm -rf "$tmpdir"' EXIT

full_fastq="${outdir}/${run}.all.fq.gz"
if [ ! -s "$full_fastq" ]; then
    dldir="${tmpdir}/dl"; mkdir -p "$dldir"
    kingfisher get -m ena-ascp ena-ftp --check-md5sums -f fastq.gz --force \
        -r "$run" --output-directory "$dldir" --ascp-ssh-key "$ASCP_KEY"
    matches=$(find "$dldir" -type f -name "*.fastq.gz")
    count=$(echo "$matches" | wc -l)
    if [ "$count" -eq 1 ]; then mv "$matches" "$full_fastq"
    elif [ "$count" -gt 1 ]; then cat $matches > "$full_fastq"
    else echo "[ERROR]: no files for $run" >&2; exit 1; fi
else
    echo "[INFO]: reusing cached download" >&2
fi

if [ ! -s "${run}.fastq.gz" ]; then
    zipfile="${tmpdir}/dataset.zip"
    datasets download genome accession "$asm_accession" --filename "$zipfile"
    asm="${tmpdir}/${asm_accession}.fna"
    unzip -p "$zipfile" '*genomic.fna' > "$asm"
    clean_fastq="${tmpdir}/${run}.clean.fq"
    minimap2 -x map-pb -t "$threads" "$asm" "$full_fastq" \
        | cut -f 1 | sort -u | seqkit grep -f - -o "$clean_fastq" "$full_fastq"
    # rasusa 4.x dropped the 2.x `-O u` value; stdout is uncompressed by default
    rasusa reads -s "$DL_SEED" -b "$MAX_BASES" "$clean_fastq" | seqkit rename -o "${run}.fastq.gz"
fi

seqkit stats -aT "${run}.fastq.gz" > "${run}.stats.tsv"
num_bases=$(awk 'NR==2 {print $5}' "${run}.stats.tsv")
[ "$num_bases" -ge "$MIN_BASES" ] || { echo "[ERROR]: $run has $num_bases bases" >&2; exit 1; }

variant () {
  local tag=$1; shift
  "$LRGE" -t "$threads" -vv -P pb -s "$LRGE_SEED" -T "$LRGE_TARGET" -Q "$LRGE_QUERY" "$@" \
      -o "${run}.${tag}.size" "${run}.fastq.gz" 2> "${run}.${tag}.log"
  grep -F 'Estimate for' "${run}.${tag}.log" | awk '{print $NF}' > "${run}.${tag}.perread.txt"
  echo "[RESULT] $run $tag estimate=$(cat "${run}.${tag}.size") n_inf=$(grep -c inf "${run}.${tag}.perread.txt")"
}

variant clean-default
variant clean-F -F
variant clean-F-strict -F --max-overhang-ratio 0.05
echo "[DONE] $run"
