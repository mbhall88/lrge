def lrge_opts(wildcards):
    opts = ""
    if wildcards.strategy == "ava":
        opts = "-n 25000 -s rand"
    elif wildcards.strategy == "2set":
        opts = "-L 5000 -O 10000 -r"
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
    benchmark:
        BENCH / "estimate/lrge-{strategy}/{dir1}/{dir2}/{dir3}/{run}.bench.tsv"
    threads: 4
    resources:
        mem="4GB",
        runtime="15m",
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
    benchmark:
        BENCH / "estimate/mash/{dir1}/{dir2}/{dir3}/{run}.bench.tsv"
    resources:
        mem="4GB",
        runtime="15m",
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
    benchmark:
        BENCH / "estimate/genomescope/{dir1}/{dir2}/{dir3}/{run}.bench.tsv"
    resources:
        mem=lambda wildcards, attempt: f"{8* attempt}GB",
        runtime="15m",
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


expand_methods = [s for s in methods for _ in range(len(runs))]


rule combine_estimates:
    input:
        estimates=expand(
            RESULTS / "estimates/{strategy}/{dir1}/{dir2}/{dir3}/{run}/{run}.size",
            zip,
            strategy=expand_methods,
            dir1=dir1s * len(expand_methods),
            dir2=dir2s * len(expand_methods),
            dir3=dir3s * len(expand_methods),
            run=runs * len(expand_methods),
        ),
        benchmarks=expand(
            BENCH / "estimate/{strategy}/{dir1}/{dir2}/{dir3}/{run}.bench.tsv",
            zip,
            strategy=expand_methods,
            dir1=dir1s * len(expand_methods),
            dir2=dir2s * len(expand_methods),
            dir3=dir3s * len(expand_methods),
            run=runs * len(expand_methods),
        ),
    output:
        RESULTS / "estimates/estimates.tsv",
    log:
        LOGS / "combine_estimates.log",
    resources:
        mem="1GB",
        runtime="10m",
    conda:
        ENVS / "combine.yaml"
    params:
        samplesheet=samplesheet,
    script:
        SCRIPTS / "combine.py"
