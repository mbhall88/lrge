"""Assemble the v1.0.0 estimate for every accession in the paper's benchmark.

No new sequencing runs were needed. Nanopore accessions come from the #36 harness `norm` arm, which
is this code with `--normalize auto` and no filter, and nanopore cannot move under the `-P` fix
because `Platform::Nanopore` selects the preset those runs already got. PacBio accessions come from
the #38 harness `pb_preset` arm, which is the same code with `-P` reaching the aligner. The `v030`
arm of the #36 harness gives the 0.3.0 estimate for the same reads, so the two versions can be
compared without the read set being a variable.
"""

import csv
from pathlib import Path

FILTER = Path("/scratch/user/uqmhal11/lrge-filter/results/filter_benchmark.tsv")
PLATFORM = Path("/scratch/user/uqmhal11/lrge-platform/results/platform_benchmark.tsv")
OUT = Path(__file__).resolve().parents[1] / "v1_estimates.tsv"


def rows(path, arm):
    with open(path) as handle:
        return {
            row["run"]: row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["arm"] == arm
        }


def main():
    legacy = rows(FILTER, "v030")
    nanopore = rows(FILTER, "norm")
    pacbio = rows(PLATFORM, "pb_preset")

    written = {"OXFORD_NANOPORE": 0, "PACBIO_SMRT": 0}
    with open(OUT, "w") as out:
        print(
            "run\tplatform\ttrue_size\tv030_estimate\tv1_estimate"
            "\tv030_relative_error\tv1_relative_error\tsource",
            file=out,
        )
        for run, old in sorted(legacy.items()):
            platform = old["platform"]
            if platform == "PACBIO_SMRT":
                new, source = pacbio.get(run), "pb_preset"
            else:
                new, source = nanopore.get(run), "norm"
            if new is None:
                raise SystemExit(f"no v1 estimate for {run}")
            truth = float(old["truth_bp"])
            old_bp, new_bp = float(old["estimate_bp"]), float(new["estimate_bp"])
            print(
                f"{run}\t{platform}\t{truth:.0f}\t{old_bp:.0f}\t{new_bp:.0f}"
                f"\t{(old_bp - truth) / truth * 100:.4f}\t{(new_bp - truth) / truth * 100:.4f}"
                f"\t{source}",
                file=out,
            )
            written[platform] += 1

    print({k: v for k, v in written.items()}, "->", OUT)


if __name__ == "__main__":
    main()
