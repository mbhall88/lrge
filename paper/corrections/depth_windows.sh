#!/usr/bin/env bash
# Windowed depth to distinguish intra-contig depth skew from flat coverage.
# Usage: depth_windows.sh <RUN> <ASSEMBLY> <THREADS> <ont|pb> [WINDOW]
set -euxo pipefail
run=$1; asm_accession=$2; threads=${3:-8}; platform=${4:-ont}; win=${5:-1000}
case "$platform" in ont) preset=map-ont ;; pb) preset=map-pb ;; *) exit 1 ;; esac

outdir=/scratch/user/uqmhal11/lrge-issue29/$run
cd "$outdir"
export PATH="$PATH:/home/uqmhal11/sw/miniforge3/envs/LRE/bin:/home/uqmhal11/sw/miniforge3/envs/Cstriatum/bin"

tmpdir=$(mktemp -d -p "$outdir"); trap 'rm -rf "$tmpdir"' EXIT
datasets download genome accession "$asm_accession" --filename "$tmpdir/ds.zip"
unzip -p "$tmpdir/ds.zip" '*genomic.fna' > "$tmpdir/asm.fna"

minimap2 -ax "$preset" -t "$threads" "$tmpdir/asm.fna" "${run}.fastq.gz" \
  | samtools sort -@ 4 -o "${run}.bam" -
samtools index "${run}.bam"

# mean depth per window, keeping the BAM this time
samtools depth -a "${run}.bam" \
  | awk -v w="$win" '{b=int($2/w); sum[$1"\t"b]+=$3; n[$1"\t"b]++}
      END{for(k in sum) printf "%s\t%d\n", k, sum[k]/n[k]}' \
  | sort -k1,1 -k2,2n > "${run}.windepth.tsv"

python3 - "$run" "$win" <<'PY'
import sys
run,win=sys.argv[1],int(sys.argv[2])
d=[int(l.split('\t')[2]) for l in open(f'{run}.windepth.tsv')]
d.sort()
n=len(d)
q=lambda p: d[min(int(p*(n-1)),n-1)]
med=q(0.5)
print(f'[WINDEPTH] {run} windows={n} size={win}')
for p in (0.01,0.1,0.5,0.9,0.99,0.999,1.0):
    print(f'   p{p:<6} {q(p):>8}  ({q(p)/max(med,1):.1f}x median)')
top=sum(x for x in d if x > 10*med)
print(f'   bases-weighted share in windows >10x median: {top/max(sum(d),1):.1%}')
PY
echo "[DONE-WIN] $run"
