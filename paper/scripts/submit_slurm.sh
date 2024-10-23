#!/usr/bin/env bash
set -eu

JOB_NAME="snakemake_master_process"
LOG_DIR="logs"

if [[ ! -d "$LOG_DIR" ]]; then
    echo "Error: Log directory $LOG_DIR does not exist"
    exit 1
fi

MEMORY="32G"
TIME="${TIME:-1d}"
THREADS=2
BINDS="/data/scratch/projects/punim2009/"
SINGULARITY_ARGS="-B $BINDS"
DEFAULT_TMP="slurm_account=punim2009"
CMD="snakemake --sdm conda apptainer --executor slurm --jobs 2000 --default-resources $DEFAULT_TMP --slurm-init-seconds-before-status-checks=20 --rerun-incomplete --local-cores $THREADS $* --singularity-args '$SINGULARITY_ARGS'"

ssubmit -t "$TIME" -m "$MEMORY" -o "$LOG_DIR"/"$JOB_NAME".o \
    -e "$LOG_DIR"/"$JOB_NAME".e "$JOB_NAME" "$CMD" -- -c "$THREADS"