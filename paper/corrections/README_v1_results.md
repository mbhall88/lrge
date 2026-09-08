# Redrawing the paper's accuracy figures for v1.0.0

The [paper][doi]'s benchmark is 3,370 bacterial read sets with a matched RefSeq assembly, and its
accuracy figures compare LRGE against GenomeScope2, Mash and Raven on all of them. This is those
figures redrawn with the v1.0.0 estimate, and what it took to get there without running the
benchmark again.

## Where the estimates come from

No new sequencing runs were needed. Both halves were already on disk:

- **Nanopore, 2,468 accessions.** The `norm` arm of the #36 harness, which is this code with
  `--normalize auto` and no filter. Nanopore cannot move under the `-P` fix, because
  `Platform::Nanopore` selects the preset those runs already got, and all 21 nanopore accessions
  carried through the #38 rerun reproduced their #36 estimate to the base.
- **PacBio, 902 accessions.** The `pb_preset` arm of the #38 harness, which is the same code with
  `-P` reaching the aligner. See [`README_issue38.md`](README_issue38.md).
- **The 0.3.0 baseline, 3,370 accessions.** The `v030` arm of the #36 harness, so the two versions
  are compared on identical reads.

The other three methods are the paper's own numbers, unchanged, from
[`results/estimates/estimates.tsv`](../results/estimates/estimates.tsv).

`scripts/build_v1_estimates.py` assembles [`v1_estimates.tsv`](v1_estimates.tsv) and
`scripts/plot_v1_results.py` draws the figures and writes the summary tables. The style, the
one-way ANOVA and the Tukey HSD annotations follow
[`../notebooks/results.ipynb`](../notebooks/results.ipynb) so the figures can be read against the
published ones.

## The one caveat, and why it does not bite

LRGE's read sets were rebuilt from ENA for the #36 and #38 harnesses, while the other three methods'
numbers were measured on the paper's own sets. The two are equivalent: 0.3.0 on the rebuilt reads
reproduces the paper's own two-set figures to within a rounding error.

| | paper's `lrge-2set` | 0.3.0 on rebuilt reads |
|---|---|---|
| nanopore, median absolute relative error | 4.84% | 4.84% |
| nanopore, within 10% of the truth | 75.8% | 75.1% |
| PacBio, median absolute relative error | 27.1% | 27.2% |
| PacBio, within 10% of the truth | 13.2% | 13.3% |

The `lrge 0.3.0` column of `version_absolute_relative_error` is that same rebuilt baseline, so the
agreement is visible in the figure rather than only asserted here.

## What the figures say

| | median absolute relative error | within 10% of the truth |
|---|---|---|
| 0.3.0, nanopore | 4.84% | 75.1% |
| **1.0.0, nanopore** | **4.77%** | **76.4%** |
| 0.3.0, PacBio | 27.2% | 13.3% |
| **1.0.0, PacBio** | **16.1%** | **26.2%** |

Over both platforms the median goes from 7.1% to 6.4% and the runs within 10% from 58.6% to 62.9%,
which moves LRGE ahead of GenomeScope2 and Mash on both measures. Raven is still far more accurate
than any of them, and still assembles the reads to get there.

The split by cause is clean. The whole PacBio gain is the `-P` fix: with the preset still wrong,
normalizing leaves the PacBio median at 27.21%, exactly where `v030` and `never` leave it. The whole
nanopore gain is depth normalization, which engages only on the inputs it finds skew in.

## What is not redrawn

The all-vs-all strategy, because it was not re-run, so the figures carry one LRGE column where the
paper's carry two.

The CPU and memory figure. Time and memory are properties of the machine as much as of the tool, and
these runs were made on different hardware from the paper's, so drawing them beside the paper's
measurements of GenomeScope2, Mash and Raven would compare the two clusters. What v1.0.0 costs
against v0.3.0 on one machine is in the README's
[Uneven read depth](../../README.md#uneven-read-depth) section and in
[`issue36_mode_cost.tsv`](issue36_mode_cost.tsv).

## Files

- `v1_estimates.tsv` - one row per accession: the 0.3.0 and 1.0.0 estimates and their relative errors.
- `v1_method_accuracy.tsv` - the method comparison summarised by method and platform.
- `v1_version_accuracy.tsv` - the same for the two versions.
- `figures/method_absolute_relative_error.{png,pdf}` - LRGE 1.0.0 against the three other methods.
- `figures/platform_relative_error.{png,pdf}` - the same, signed, so over and under estimation separate.
- `figures/version_absolute_relative_error.{png,pdf}` - 0.3.0 against 1.0.0 on identical reads.

[doi]: https://doi.org/10.1093/bioinformatics/btaf593
