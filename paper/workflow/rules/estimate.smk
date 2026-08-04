def lrge_opts(wildcards):
    opts = ["-vv"]

    seed = config["lrge"].get("seed")
    if seed is not None:
        opts.extend(["-s", f"{seed}"])

    if wildcards.strategy == "ava":
        num = config["lrge"]["ava"]
        opts.extend(["-n", f"{num}"])
    elif wildcards.strategy == "2set":
        target = config["lrge"]["twoset"]["target"]
        query = config["lrge"]["twoset"]["query"]
        opts.extend(["-T", f"{target}", "-Q", f"{query}"])
    else:
        raise UnimplementedError(f"Unknown strategy: {wildcards.strategy}")

    platform = samplesheet.loc[wildcards.run, "Instrument Platform"]
    if "OXFORD" in platform:
        platform = "ont"
    elif "PACBIO" in platform:
        platform = "pb"
    else:
        raise UnimplementedError(f"Unknown platform: {platform}")

    opts.extend(["-P", platform])

    return " ".join(opts)


rule estimate_lrge:
    input:
        fastq=rules.download.output.fastq,
        bin=WORKFLOW / config["lrge"]["bin"],
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
        runtime=lambda wildcards, attempt: f"{attempt*2}h",
    shadow:
        "shallow"
    params:
        opts=lrge_opts,
    shell:
        "{input.bin} -t {threads} {params.opts} -o {output.size} {input.fastq} 2> {log}"


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

rule estimate_raven:
    input:
        fastq=rules.download.output.fastq,
    output:
        size=RESULTS / "estimates/raven/{dir1}/{dir2}/{dir3}/{run}/{run}.size",
    log:
        LOGS / "estimate_raven/{dir1}/{dir2}/{dir3}/{run}.log",
    group:
        "estimate"
    benchmark:
        BENCH / "estimate/raven/{dir1}/{dir2}/{dir3}/{run}.bench.tsv"
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: (30*attempt) * 1_000,
        runtime=lambda wildcards, attempt: f"{6*attempt}h",
    container:
        "docker://quay.io/biocontainers/raven-assembler:1.8.3--h43eeafb_1"
    shadow:
        "shallow"
    shell:
        r"""
        (raven -p 0 -t {threads} {input.fastq} | \
            grep -v '^>' | \
            tr -d '\n' | \
            wc -c > {output.size}) 2> {log}
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
