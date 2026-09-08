# Re-deriving the reported interval quantiles

Issue [#38]. LRGE prints an interval either side of its estimate, taken as two percentiles of the
per-read estimates. The paper fitted one pair for both platforms, the 15th and 65th percentiles, and
reported that it contained the true genome size on about 92% of its benchmark. Depth-aware read
selection reshapes the per-read distribution, so the ticket was to re-derive the pair on the
post-fix benchmark and measure what it covers now.

The short answer is that one pair cannot serve both platforms. The shipped pair covers 87.2% of the
benchmark, and that average is made of 97.0% on nanopore and 60.2% on PacBio. There are now two
pairs, each the narrowest that covers 95% of its own platform's runs. On nanopore the interval comes
out 18% narrower than the one shipped today.

Two bugs turned up on the way. `-P/--platform` never reached the aligner from the command line, so
every release up to 1.0.0 overlapped PacBio reads with the nanopore preset, and the PacBio half of
the benchmark had to be re-run before anything could be fitted on it. The per-read estimates `-vv`
prints were the ones with internal matches removed, whatever the run was asked to do with them. Both
are fixed here.

## What the paper did, and what it picked

The procedure is in `paper/notebooks/per_read_estimates.ipynb`. Per-read estimates were pulled out
of `-vv` logs, divided by the true genome size, and reduced to percentiles per accession. The
interval's width was fixed at 50 percentile points and the offset swept: for each `n`, how often does
`[n, n+50]` contain the truth. The notebook's stored output puts the two-set peak at an offset of
13, covering 92.6% at a mean width of 0.531 times the genome, and the narrowest interval at an
offset of 16, covering 92.3% at 0.524. The shipped offset of 15 sits between them, at 92.5% and
0.525, which is the trade the paper took.

## The runs

The [#36] harness ran every arm at `-vv` and archived the stderr, so the whole per-read distribution
for all 3,370 accessions is on disk and the nanopore fit needed no new sequencing runs. Its `norm`
arm is this code as it ships, `--normalize auto` with no filter.

Reading the distribution back out of the logs was checked against the tool rather than trusted. For
the `v030`, `never` and `norm` arms the percentiles recovered from the logs reproduce the median
LRGE wrote to its output file and both interval bounds it logged, on all 3,370 accessions, with no
mismatches. Coverage computed the other way round, from the Mbp bounds in the tool's own log line,
comes to 87.18% against the 87.15% from the recovered percentiles, the difference being the two
decimal places the log rounds to.

The two `-F` arms are left out. Their trace lines carry the estimate with internal matches removed
whether or not the run removed them, which is the second bug above, so their logs do not describe
the answers those runs gave.

## What the shipped pair covers now

| arm | what it is | all 3,370 | 2,468 nanopore | 902 PacBio |
|---|---|---|---|---|
| `v030` | the last release | 86.8% | 96.5% | 60.3% |
| `never` | this code, both new mechanisms off | 86.9% | 96.6% | 60.2% |
| `norm` | this code as it ships | 87.2% | 97.0% | 60.2% |

Depth-aware selection is not what moved the number. It is worth 0.4 points, and the benchmark was
already not giving 92% before any of this work started. The PacBio column here is the buggy preset,
which is what those releases actually did.

The gap against the paper's 92% is unexplained. It is not the estimator: the paper's own
per-accession two-set estimates and the rebuilt `v030` run agree to a median of 0.75% per accession,
with median relative sizes of 1.0518 and 1.0514. The per-read table the notebook computed its figure
from is not on disk, so it could not be traced further.

## `-P` never reached the aligner

`lrge/src/main.rs` parsed `--platform`, printed it in the debug line, and never passed it to a
builder. `Builder::platform` defaults to nanopore, so every command-line run since the workspace
split overlapped with `ava-ont`. Library callers were unaffected. `git log -S".platform(" --
lrge/src/main.rs` returns nothing: the call was never written.

That makes the 902 PacBio runs in the #36 benchmark unusable for fitting a PacBio interval, because
they are PacBio reads mapped as if they were nanopore. All 902 were re-run on the fixed binary,
under both presets, from one copy of the reads, along with 21 nanopore accessions carried as a check
on the claim that nanopore cannot move. The `-P ont` arm is the control: it has to reproduce the
estimate already on disk, or the read set was rebuilt differently and the two benchmarks do not
compare.

921 of the 923 accessions rebuilt to the same read count and base count as in #36, and 920
reproduced its estimate to the base, the 21 nanopore accessions among them. Nanopore cannot move
under this fix and does not. The one accession that matched on reads and bases but not on
the estimate, `SRR10672468`, differs by 426 bp in 4.87 Mbp, which is minimap2 giving slightly
different mappings across threads. The two that differ, `SRR5066061` and `SRR6515889`, came back
from ENA with a different number of reads for near enough the same number of bases, so the archive
served a different file composition and their read sets are not the #36 ones.

Giving PacBio reads the PacBio preset is worth a lot on its own:

| | median estimate over truth | within 10% | over 2x | mean abs log2 error |
|---|---|---|---|---|
| `-P ont`, which is what shipped | 1.266 | 13.4% | 10.4% | 0.520 |
| `-P pb`, with the flag wired up | 1.155 | 26.2% | 4.5% | 0.351 |

over the 902 PacBio accessions. It does not close the gap to nanopore, which sizes 76.3% of its runs
within 10%, but the flag was doing nothing at all.

## The fit

Coverage and width trade against each other smoothly, so the fit needs a target rather than an
optimum. The paper's fixed width of 50 points is one way to pick one, and `issue38_width50_scan.tsv`
keeps that scan so the shipped pair can be compared on the rule that chose it.
`issue38_coverage_grid.tsv` drops the fixed width and reports every pair from 1-51 to 49-99.

The pairs below are the narrowest reaching 95% coverage on their own platform. Every third accession
was held out as a check that the fit is not reading noise.

| | pair | coverage | mean width | median width | fitted on two thirds | held-out third |
|---|---|---|---|---|---|---|
| nanopore, shipped | 15-65 | 97.0% | 0.471 | 0.431 | 96.5% | 97.9% |
| nanopore, fitted | 20.5-63.5 | 95.2% | 0.385 | 0.354 | 94.6% | 96.5% |
| PacBio, shipped | 15-65 | 80.2% | 0.616 | 0.481 | | |
| PacBio, fitted | 1-61.5 | 95.0% | 0.988 | 0.851 | 94.8% | 95.4% |

Widths are multiples of the true genome size. Both PacBio rows are on the fixed preset; under the
one that shipped, the same 15-65 pair covers 60.2% at a mean width of 0.832.

Neither fit sits on a knife edge or depends on genome size. Bootstrapping over accessions puts
nanopore at [94.4%, 96.0%] and PacBio at [93.5%, 96.5%]. Across genome-size quartiles nanopore runs
94.5%, 95.3%, 95.3% and 95.8%, and PacBio 93.4%, 94.7%, 96.4% and 95.6%. Moving either end of either
pair by a percentile point moves coverage by less than a point.

What the fit is really measuring is where the truth falls in a run's own estimates:
`issue38_truth_percentile.tsv` has it per accession. The median nanopore run puts the true size at
the 46th percentile of its per-read estimates, and the median PacBio run at the 29th. Under the
preset that shipped, the median PacBio run put it at the 20th.

## What this changes

- `liblrge::estimate` gains `NANOPORE_LOWER_QUANTILE`, `NANOPORE_UPPER_QUANTILE`,
  `PACBIO_LOWER_QUANTILE` and `PACBIO_UPPER_QUANTILE`, and `Platform::interval_quantiles` hands back
  the pair for a platform. `LOWER_QUANTILE` and `UPPER_QUANTILE` keep their values and are
  deprecated.
- `--q1` and `--q3` take their defaults from `-P`.
- The interval is labelled `95% CI` when it is the fitted pair for the platform, and named as
  `q<lower>-q<upper>` when the caller gives a pair of their own, because an arbitrary pair has no
  measured coverage.

The PacBio lower bound is the first percentile, which is about forty estimates in at the default
read counts but falls between the two smallest estimates a run makes at `-Q 100`. A small PacBio
run's lower bound is one read's estimate and should be read as such.

## Tables

- `issue38_width50_scan.tsv` - the paper's scan, width fixed at 50 percentile points.
- `issue38_coverage_grid.tsv` - every pair from 1-51 to 49-99, with coverage and width.
- `issue38_target_fit.tsv` - the narrowest pair reaching 90, 92, 95 and 97%, with held-out coverage.
- `issue38_truth_percentile.tsv` - where the true size falls in each run's per-read estimates.
- `issue38_preset_change.tsv` - every accession in the re-run, under both presets.

The first four carry a `pacbio_nanopore_preset` group alongside `pacbio`, which is the same 902
accessions under the preset that shipped, so the two can be compared directly.

[#36]: https://github.com/mbhall88/lrge/issues/36
[#38]: https://github.com/mbhall88/lrge/issues/38
