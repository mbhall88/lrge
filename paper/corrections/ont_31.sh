#!/usr/bin/env bash
# The 13 remaining ONT sub-0.5x failures: how much of each is mechanism 2?
# Also measures the internal-match fraction of each run's overlap set, which is the
# first half of the "depth skew generates internal matches" hypothesis (NOTES.md 4.4).
set -euxo pipefail

run=$1; asm_accession=$2; threads=${3:-8}

WORKDIR=/scratch/user/uqmhal11/lrge-issue29
LRGE=/scratch/user/uqmhal11/lrge/target/release/lrge
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
fi

if [ ! -s "${run}.fastq.gz" ]; then
    zipfile="${tmpdir}/dataset.zip"
    datasets download genome accession "$asm_accession" --filename "$zipfile"
    asm="${tmpdir}/${asm_accession}.fna"
    unzip -p "$zipfile" '*genomic.fna' > "$asm"
    clean_fastq="${tmpdir}/${run}.clean.fq"
    minimap2 -x map-ont -t "$threads" "$asm" "$full_fastq" \
        | cut -f 1 | sort -u | seqkit grep -f - -o "$clean_fastq" "$full_fastq"
    rasusa reads -s "$DL_SEED" -b "$MAX_BASES" "$clean_fastq" | seqkit rename -o "${run}.fastq.gz"
fi

seqkit stats -aT "${run}.fastq.gz" > "${run}.stats.tsv"
num_bases=$(awk 'NR==2 {print $5}' "${run}.stats.tsv")
echo "[INPUT] $run num_seqs=$(awk 'NR==2{print $4}' "${run}.stats.tsv") num_bases=$num_bases"
[ "$num_bases" -ge "$MIN_BASES" ] || { echo "[ERROR]: $run has $num_bases bases" >&2; exit 1; }

variant () {
  local tag=$1; shift
  "$LRGE" -t "$threads" -vv -P ont -s "$LRGE_SEED" -T "$LRGE_TARGET" -Q "$LRGE_QUERY" "$@" \
      -o "${run}.${tag}.size" "${run}.fastq.gz" 2> "${run}.${tag}.log"
  grep -F 'Estimate for' "${run}.${tag}.log" | awk '{print $NF}' > "${run}.${tag}.perread.txt"
  echo "[RESULT] $run $tag estimate=$(cat "${run}.${tag}.size") n_inf=$(grep -c inf "${run}.${tag}.perread.txt")"
}

# default run keeps its overlaps so the internal-match fraction can be measured
paftmp="${tmpdir}/paf"; mkdir -p "$paftmp"
"$LRGE" -t "$threads" -vv -P ont -s "$LRGE_SEED" -T "$LRGE_TARGET" -Q "$LRGE_QUERY" \
    --keep-temp --temp "$paftmp" -o "${run}.fixed-default.size" "${run}.fastq.gz" \
    2> "${run}.fixed-default.log"
grep -F 'Estimate for' "${run}.fixed-default.log" | awk '{print $NF}' > "${run}.fixed-default.perread.txt"
echo "[RESULT] $run fixed-default estimate=$(cat "${run}.fixed-default.size") n_inf=$(grep -c inf "${run}.fixed-default.perread.txt")"

paf=$(find "$paftmp" -name 'overlaps.paf' | head -1)
if [ -n "$paf" ]; then
  awk -F'\t' -v OFS='\t' '{
    qs=$3; qe=$4; ql=$2; ts=$8; te=$9; tl=$7;
    if ($5=="+") { a=(qs<ts?qs:ts); b=((ql-qe)<(tl-te)?(ql-qe):(tl-te)) }
    else         { a=(qs<(tl-te)?qs:(tl-te)); b=((ql-qe)<ts?(ql-qe):ts) }
    oh=a+b; ml=((qe-qs)>(te-ts)?(qe-qs):(te-ts));
    if (ml>0) { n++; if (oh > 0.2*ml) i++; if (oh > 0.05*ml) i5++ }
  } END{ printf "[PAF] '"$run"' n_overlaps=%d internal_at_0.2=%d (%.2f%%) internal_at_0.05=%d (%.2f%%)\n", n, i, 100*i/n, i5, 100*i5/n }' "$paf"
fi
rm -rf "$paftmp"

variant fixed-F -F
variant fixed-F-strict -F --max-overhang-ratio 0.05
echo "[DONE] $run"
