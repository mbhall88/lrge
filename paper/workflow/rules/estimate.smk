def lrge_opts(wildcards):
    opts = ""
    if wildcards.strategy == "ava":
        opts = "-n 25000 -s rand -vv"
    elif wildcards.strategy == "2set":
        opts = "-L 5000 -O 10000 -r -v"
    else:
        raise UnimplementedError(f"Unknown strategy: {wildcards.strategy}")

    platform = samplesheet.loc[wildcards.run, "Instrument Platform"]
    if "OXFORD" in platform:
        opts += " -P ont"
    elif "PACBIO" in platform:
        opts += " -P pb"
    else:
        raise UnimplementedError(f"Unknown platform: {platform}")

    return opts


rule estimate_lrge:
    input:
        fastq=rules.download.output.fastq,
        script=SCRIPTS / "lrge-{strategy}.py",
    output:
        size=RESULTS / "estimates/lrge-{strategy}/{dir1}/{dir2}/{dir3}/{run}/{run}.size",
    log:
        LOGS / "estimate_lrge/{strategy}/{dir1}/{dir2}/{dir3}/{run}.log",
    group:
        "estimate"
    benchmark:
        BENCH / "estimate/lrge-{strategy}/{dir1}/{dir2}/{dir3}/{run}.bench.tsv"
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: (4**attempt) * 1_000,
        runtime=lambda wildcards, attempt: f"{attempt}h",
    conda:
        ENVS / "lrge.yaml"
    shadow:
        "shallow"
    params:
        opts=lrge_opts,
    shell:
        "python {input.script} -t {threads} {params.opts} {input.fastq} > {output.size} 2> {log}"


rule estimate_mash:
    input:
        fastq=rules.download.output.fastq,
    output:
        size=RESULTS / "estimates/mash/{dir1}/{dir2}/{dir3}/{run}/{run}.size",
    log:
        LOGS / "estimate_mash/{dir1}/{dir2}/{dir3}/{run}.log",
    group:
        "estimate"
    benchmark:
        BENCH / "estimate/mash/{dir1}/{dir2}/{dir3}/{run}.bench.tsv"
    resources:
        mem_mb=lambda wildcards, attempt: (4**attempt) * 1_000,
        runtime=lambda wildcards, attempt: f"{attempt}h",
    conda:
        ENVS / "mash.yaml"
    shadow:
        "shallow"
    params:
        opts="-r",
        min_copies=config["mash"]["min_copies"],
        sketch_size=config["mash"]["sketch_size"],
    shell:
        r"""
        > {log}
        mash sketch {params.opts} -m {params.min_copies} -s {params.sketch_size} -o $(mktemp) {input.fastq} 2>&1 | \
            tee -a {log} | \
            rg -o -r '$1' 'genome size: (\d+.*)$' | \
            python -c "import sys;print(float(sys.stdin.read().strip()))" > {output.size} 2> {log}
        """


rule estimate_genomescope:
    input:
        fastq=rules.download.output.fastq,
        script=SCRIPTS / "genomescope.py",
    output:
        size=RESULTS / "estimates/genomescope/{dir1}/{dir2}/{dir3}/{run}/{run}.size",
    log:
        LOGS / "estimate_genomescope/{dir1}/{dir2}/{dir3}/{run}.log",
    group:
        "estimate"
    benchmark:
        BENCH / "estimate/genomescope/{dir1}/{dir2}/{dir3}/{run}.bench.tsv"
    resources:
        mem_mb=lambda wildcards, attempt: (4**attempt) * 1_000,
        runtime=lambda wildcards, attempt: f"{attempt}h",
    conda:
        ENVS / "genomescope.yaml"
    shadow:
        "shallow"
    params:
        ploidy=config["genomescope"]["ploidy"],
        min_copies=config["genomescope"]["min_copies"],
        kmer_size=config["genomescope"]["kmer_size"],
        max_copies=config["genomescope"]["max_copies"],
    shell:
        r"""
        # we do this as we append to the log file in the pipeline to avoid overwrites
        > {log}

        (python {input.script} -p {params.ploidy} -m {params.min_copies} \
            -k {params.kmer_size} -M {params.max_copies} --tmp $(mktemp -d) {input.fastq} 2>&1 | \
            tee -a {log} | \
            rg -o -r '$1' 'len:(\d+)') > {output.size} 2>> {log}
        """


rule combine_estimates:
    input:
        estimates=combine_estimate_paths,
        benchmarks=combine_benchmark_paths,
        stats=combine_stats_paths,
    output:
        RESULTS / "estimates/estimates.tsv",
    log:
        LOGS / "combine_estimates.log",
    resources:
        mem="1GB",
        runtime="30m",
    conda:
        ENVS / "combine.yaml"
    params:
        samplesheet=samplesheet,
    script:
        SCRIPTS / "combine.py"
