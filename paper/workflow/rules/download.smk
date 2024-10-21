rule download:
    output:
        fastq=temp(RESULTS / "downloads/{dir1}/{dir2}/{dir3}/{run}/{run}.fastq.gz"),
    log:
        LOGS / "download/{dir1}/{dir2}/{dir3}/{run}.log"
    container:
        "docker://quay.io/biocontainers/fastq-dl:2.0.4--pyhdfd78af_0"
    resources:
        mem="1GB",
        runtime="12h"
    shadow:
        "shallow"
    params:
        outdir=lambda wildcards, output: Path(output.fastq).parent,
        opts="--provider sra"
    shell:
        """
        exec 2> {log}
        tmpdir=$(mktemp -d)
        trap 'rm -rf $tmpdir' EXIT

        fastq-dl {params.opts} -a {run} -o "$tmpdir"

        # Find all files matching the pattern <run>*.fastq.gz in the output directory
        matches=$(find "$tmpdir" -type f -name "{run}*.fastq.gz")

        # Count how many files match the pattern
        count=$(echo "$matches" | wc -l)

        if [ $count -eq 1 ]; then
            mv "$matches" {output.fastq}
            echo "[SUCCESS]: Renamed $matches to {output.fastq}" >&2
        elif [ $count -gt 1 ]; then
            # If more than one match, concatenate them
            cat $matches > {output.fastq}
            
            echo "[SUCCESS]: Concatenated the following files into {output.fastq}:" >&2
            echo "$matches" >&2
        else
            # If no matches found, print an error and exit
            echo "[ERROR]: No matching files found for {run}" >&2
            ls -la "$tmpdir" >&2
            exit 1
        fi
        """
