rule download:
    output:
        fastq=temp(RESULTS / "downloads/{dir1}/{dir2}/{dir3}/{run}/{run}.fastq.gz"),
        stats=RESULTS / "stats/{dir1}/{dir2}/{dir3}/{run}.stats.tsv",
        info=RESULTS / "stats/{dir1}/{dir2}/{dir3}/{run}.info.tsv",
    log:
        LOGS / "download/{dir1}/{dir2}/{dir3}/{run}.log",
    # group:
    #     "estimate"
    threads: 4
    conda:
        ENVS / "download.yaml"
    resources:
        mem_mb=8_000,
        runtime="12h",
    shadow:
        "shallow"
    params:
        ascp_ssh_key=config["ascp_ssh_key"],
        max_bases=config["download"]["max_bases"],
        min_bases=config["download"]["min_bases"],
        seed=config["download"]["seed"],
        asm_accession=lambda wildcards: samplesheet.loc[wildcards.run, "Assembly Accession"],
        platform=lambda wildcards: "map-ont" if samplesheet.loc[wildcards.run, "Instrument Platform"] == "OXFORD_NANOPORE" else "map-pb",
    script:
        SCRIPTS / "download.sh"
