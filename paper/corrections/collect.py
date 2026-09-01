#!/usr/bin/env python3
"""Collect one row per accession/mode from the PR #45 benchmark run directories."""
import re, sys, os, csv

BENCH = "/scratch/user/uqmhal11/lrge-issue34-benchmark"
RUNS = os.path.join(BENCH, "runs")
TRUTH_TSV = "/scratch/user/uqmhal11/lrge/paper/corrections/rerun_estimates.tsv"

OUTLIERS = """DRR213976 SRR10259778 SRR10353548 SRR10388020 SRR10861747 SRR10861751
SRR24489322 SRR26465521 SRR26465523 SRR26465524 SRR26465526 SRR26465560 SRR26465563
SRR26715165 SRR26715166""".split()
CONTROLS = {"SRR12247681": "OXFORD_NANOPORE", "SRR8618952": "PACBIO_SMRT"}

truth, platform, rerun_default = {}, {}, {}
with open(TRUTH_TSV) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        truth[row["run"]] = row["truth_size"]
        platform[row["run"]] = row["platform"]
        rerun_default[row["run"]] = row["rerun_default"]

UNITS = {"bp": 1, "kbp": 1e3, "Mbp": 1e6, "Gbp": 1e9}

SKEW = re.compile(
    r"Depth skew detected \(99\.9th percentile minimizer count is ([\d.]+)x the median; "
    r"sampled reads: (\d+)\); depth normalization retained (\d+) of (\d+) reads")
FORCED = re.compile(r"Depth normalization forced; retained (\d+) of (\d+) reads")
NOOVL = re.compile(r"(\d+) \([\d.]+%\) query read\(s\) did not overlap")
EST = re.compile(r"Estimated genome size: ([\d.]+) (\w+)(?: \(IQR: ([\d.]+) (\w+) - ([\d.]+) (\w+)\))?")


def to_bp(value, unit):
    return str(int(round(float(value) * UNITS[unit])))


def read(path, default="NA"):
    try:
        with open(path) as fh:
            return fh.read().strip()
    except OSError:
        return default


def parse_time(path):
    elapsed, rss = "NA", "NA"
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith("Elapsed (wall clock) time"):
                    elapsed = line.split(": ", 1)[1]
                elif line.startswith("Maximum resident set size"):
                    rss = line.rsplit(" ", 1)[1]
    except OSError:
        pass
    return elapsed, rss


FIELDS = ["accession", "class", "platform", "truth_bp", "mode", "commit", "estimate_bp",
          "estimate_over_truth", "interval_low_bp", "interval_high_bp",
          "normalization_engaged", "skew_score", "detection_sample_reads", "retained_reads",
          "total_reads", "queries_without_overlaps", "elapsed", "max_rss_kb", "exit_status",
          "stdout_path", "stderr_path"]

rows = []
for acc in sorted(os.listdir(RUNS)):
    for mode in sorted(os.listdir(os.path.join(RUNS, acc))):
        d = os.path.join(RUNS, acc, mode)
        if not os.path.isdir(d):
            continue
        err = read(os.path.join(d, "lrge.err"), "")
        est = read(os.path.join(d, "estimate.out"))
        row = dict.fromkeys(FIELDS, "NA")
        row.update(
            accession=acc,
            cls="outlier" if acc in OUTLIERS else "control",
            platform=platform.get(acc, "NA"),
            truth_bp=truth.get(acc, "NA"),
            mode=mode,
            commit=read(os.path.join(d, "commit.txt")),
            estimate_bp=est if est else "NA",
            exit_status=read(os.path.join(d, "exit_status.txt")),
            stdout_path=os.path.join(d, "estimate.out"),
            stderr_path=os.path.join(d, "lrge.err"),
        )
        row["class"] = row.pop("cls")
        if est.isdigit() and row["truth_bp"] not in ("NA", ""):
            row["estimate_over_truth"] = f"{int(est) / int(row['truth_bp']):.3f}"
        m = SKEW.search(err)
        if m:
            row.update(normalization_engaged="yes", skew_score=m.group(1),
                       detection_sample_reads=m.group(2), retained_reads=m.group(3),
                       total_reads=m.group(4))
        elif FORCED.search(err):
            m = FORCED.search(err)
            row.update(normalization_engaged="forced", retained_reads=m.group(1),
                       total_reads=m.group(2))
        else:
            row["normalization_engaged"] = "no"
        m = NOOVL.search(err)
        row["queries_without_overlaps"] = m.group(1) if m else "0"
        m = EST.search(err)
        if m and m.group(3):
            row["interval_low_bp"] = to_bp(m.group(3), m.group(4))
            row["interval_high_bp"] = to_bp(m.group(5), m.group(6))
        row["elapsed"], row["max_rss_kb"] = parse_time(os.path.join(d, "time.txt"))
        rows.append(row)

order = {"never": 0, "auto": 1, "auto_F": 2, "always": 3}
rows.sort(key=lambda r: (r["class"], r["accession"], order.get(r["mode"], 9)))
w = csv.DictWriter(sys.stdout, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
w.writeheader()
w.writerows(rows)
