rule download:
    output:
        fastq=temp(RESULTS / "downloads/{dir1}/{dir2}/{dir3}/{run}/{run}.fastq.gz"),
        stats=RESULTS / "stats/{dir1}/{dir2}/{dir3}/{run}.stats.tsv",
    log:
        LOGS / "download/{dir1}/{dir2}/{dir3}/{run}.log",
    # group:
    #     "estimate"
    conda:
        ENVS / "download.yaml"
    resources:
        mem_mb=6_000,
        runtime="12h",
    shadow:
        "shallow"
    params:
        ascp_ssh_key=config["ascp_ssh_key"],
        opts="-m ena-ascp ena-ftp --check-md5sums -f fastq.gz --force -r {run} --debug",
        max_bases=config["download"]["max_bases"],
        seed=config["download"]["seed"],
    script:
        SCRIPTS / "download.sh"
