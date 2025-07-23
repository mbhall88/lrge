import sys

sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path
import pandas as pd

samplesheet = snakemake.params.samplesheet.copy(deep=True)

drop_cols = [
    "Assembly BioSample Project name  ",
    "Assembly BioSample Strain ",
    "Assembly Name",
]
samplesheet = samplesheet.drop(columns=drop_cols)
rename_cols = {
    "Assembly BioSample Accession": "biosample",
    "Assembly Stats Total Sequence Length": "true_size",
    "Assembly Accession": "asm_accession",
    "Organism Name": "organism",
    "Organism Taxonomic ID": "taxid",
    "Assembly BioProject Accession": "bioproject",
    "Organism Infraspecific Names Strain": "strain",
    "Assembly Sequencing Tech": "asm_seq_tech",
    "Assembly Stats Total Number of Chromosomes": "n_chromosomes",
    "Assembly Stats Genome Coverage": "asm_coverage",
    "Instrument Platform": "platform",
    "CheckM completeness": "checkm_completeness",
    "CheckM contamination": "checkm_contamination",
    "CheckM completeness percentile": "checkm_completeness_percentile",
    "Assembly Release Date": "release_date",
    "Library Selection": "library_selection",
    "Library Source": "library_source",
    "Library Strategy": "library_strategy",
}
samplesheet = samplesheet.rename(columns=rename_cols)
samplesheet.index.names = ["run"]

estimates = dict()
for p in map(Path, snakemake.input.estimates):
    run = p.name.split(".")[0]
    method = p.parts[-6]
    estimates[(run, method)] = p

stats = dict()
keep_stats_cols = ["num_seqs", "sum_len", "avg_len", "Q2", "AvgQual"]
for p in map(Path, snakemake.input.stats):
    run = p.name.split(".")[0]
    df = pd.read_csv(p, sep="\t")
    # remove columns not in keep_stats_cols
    df = df.loc[:, keep_stats_cols]
    # the stats df only has one row, extract it as a list
    stats[run] = df.iloc[0].to_list()

rows = []
bad_runs = []
failed = False

for p in map(Path, snakemake.input.benchmarks):
    run = p.name.split(".")[0]
    method = p.parts[-5]
    estimate_file = estimates[(run, method)]
    est = estimate_file.read_text().strip()
    try:
        est = float(est)
    except ValueError:
        print(f"estimate is empty for {run}", file=sys.stderr)
        failed = True
        bad_runs.append(run)
        continue

    true_size = samplesheet.loc[run, "true_size"]

    # samplesheet.loc[run, "estimate"] = est
    # samplesheet.loc[run, "method"] = method
    # get the row from the samplesheet
    row = samplesheet.loc[run]
    row = row.to_list()
    row.extend([run, est, method])

    relative_size = est / true_size
    relative_error = (est - true_size) / true_size * 100

    # samplesheet.loc[run, "relative_size"] = relative_size
    # samplesheet.loc[run, "relative_error"] = relative_error
    row.extend([relative_size, relative_error])

    bench = pd.read_csv(p, sep="\t")
    cpu_time = bench["cpu_time"].mean()
    memory = bench["max_rss"].mean()
    # samplesheet.loc[run, "cpu_time"] = cpu_time
    # samplesheet.loc[run, "memory_mb"] = memory
    row.extend([cpu_time, memory])

    row.extend(stats[run])

    rows.append(row)

if failed:
    print("Failed to read estimates for the following runs:", file=sys.stderr)
    for r in bad_runs:
        print(r, file=sys.stderr)
    sys.exit(1)

outsheet = pd.DataFrame(
    rows,
    columns=list(samplesheet.columns)
    + [
        "run",
        "estimate",
        "method",
        "relative_size",
        "relative_error",
        "cpu_time",
        "memory_mb",
        "stats_num_seqs",
        "stats_sum_len",
        "stats_avg_len",
        "stats_median_len",
        "avg_qual",
    ],
)

outsheet.to_csv(snakemake.output[0], sep="\t", index=False)
