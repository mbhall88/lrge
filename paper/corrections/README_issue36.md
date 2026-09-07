# Fitting the depth-normalization constants on the benchmark

Issue [#36]. Two constants in `liblrge/src/depth_skew.rs` were picked by argument rather than
measurement: `SKEW_THRESHOLD`, how far the 99.9th percentile minimizer count has to run ahead of
the median before an input is called skewed, and `RETENTION_TARGET_MULTIPLIER`, the multiple of the
input's median depth that reads are kept down to. This is how they were fitted, and what the
benchmark said.

The short answer: neither moves. Normalizing at all is a large, reproducible win, and where the
constants sit inside a wide plateau is not.

What came out of the same runs is that the failures normalization leaves behind have a different
cause, that the cause is measurable during the overlap pass, and that `-F` can therefore be switched
on per run rather than per user. The second half of this document is that mechanism, and the
benchmark's verdict on it: it belongs in the tool, but not on by default, and the threshold it needs
is not the one a sample of hard cases suggested. From
[The next mechanism, and a limit on this fit](#the-next-mechanism-and-a-limit-on-this-fit).

## The runs

The paper's benchmark is 3,370 bacterial long-read runs with a matched RefSeq assembly. Its
workflow marks the FASTQ files `temp()`, so none survived and all 3,370 were rebuilt from ENA:
fetched with kingfisher, cleaned by keeping only reads that map to the run's own assembly, and
downsampled to 1 Gbp with the paper's seed. `workflow/envs/benchmark.yaml` in the harness pins the
same versions as `paper/workflow/envs/download.yaml`, which matters more than it sounds: `rasusa`'s
subsampling RNG changed between 2.x and 4.x, and under 4.1 the accessions kept on disk from the #29
work do not reproduce their published read counts.

Every accession completed. Twelve failed once with a truncated download and succeeded on retry.

## Why seven runs per accession, and not sixty-six

`SKEW_THRESHOLD` enters the code in exactly one place, `skewed: score >= SKEW_THRESHOLD`, and the
score itself does not depend on it or on the retention multiplier. A run scoring below a candidate
threshold behaves exactly as `--normalize never` does. So

```
estimate(T, C) = normalized[C]  if score >= T
                 legacy         otherwise
```

reconstructs the whole grid from one unnormalized run plus one normalized run per candidate
multiplier. That is seven runs per accession rather than one per cell. The normalized runs use
binaries built with `SKEW_THRESHOLD = 0` so that every run reports its score and normalizes.

This was checked before it was relied on, and again afterwards. Four accessions straddling the
shipped threshold reproduced the stock binary exactly, including `SRR13183064`, which scores 4.0 and
so must fall back to the unnormalized estimate. Then 23 further accessions were run through the
stock binary in `auto` mode and compared against their reconstructed values: 23 exact, 0 differing.

## The grid

Thresholds against multipliers (1, 2, 3, 4, 6, 8), scored on all 3,370 runs. The grid sweeps
seventeen thresholds from 0 to 100 plus a no-normalization baseline; the selection table narrows to
the eleven from 4 to 48 that are plausible choices:

- `issue36_benchmark_sweep.tsv` — one row per accession per variant, 23,590 rows. The raw material.
- `issue36_constant_grid.tsv` — the grid over all 3,370 runs.
- `issue36_constant_grid_holdout.tsv` — the same grid, holding out every third accession.
- `issue36_selection.tsv` — rescues against regressions, with two bootstraps over accessions:
  `boot_net_*` on the rescue-minus-regression count, and `mean_err_vs_shipped` with `boot_err_*`
  pairing each cell's mean |log2| error against the shipped pair.
- `issue36_low_depth_multiplier.tsv` — #56's low-depth ladder, rebuilt across the six multipliers.
- `issue36_timing_measured.tsv` — the stock binary timed in both modes, three replicates.
- `issue36_failure_outcomes.tsv` — every run under half its true size, before or after.
- `issue36_movers.tsv` — the 186 runs the detector fires on, classed by what happened to them.
- `issue36_filter_contained.tsv` — 27 accessions with and without `-F`, crossed with normalization.
- `issue36_internal_match_probe.tsv` — what `-F` would remove, measured without removing it.

## What it says

**Normalizing works.** Against no normalization, the shipped configuration takes the mean |log2|
error over the benchmark from 0.2372 to 0.2188, the runs estimating under half their true size from
25 to 13, and the runs landing within 10% of the truth from 1,978 to 2,006.

**The constants sit on a plateau.** The best cell by mean error, a threshold of 24 with a multiplier
of 1, beats the shipped pair by 0.00041, on a paired 95% bootstrap interval of [-0.00157, +0.00068]
(`issue36_selection.tsv`, the `mean_err_vs_shipped` and `boot_err_*` columns). Which
cell wins at all depends on the accuracy band you score on: the tightest bands favour a threshold of
4 with a multiplier of 8, the widest a threshold of 6 with a multiplier of 4, and none of them beats
the shipped pair by more than about one percent of the runs.

**The benchmark alone would raise the multiplier, and it should not.** Going from 2 to 4 halves the
already-correct runs that normalization pushes out of the band, 13 down to 6, with no loss of
rescues. But all thirteen are estimates drifting up between 2 and 18 percent, none of them landing
more than 19 percent above the truth, and a larger multiplier buys that by normalizing less. The
benchmark holds three skewed runs with a profile median depth under four, so it cannot see what
normalizing less costs where depth is thin. #56's ladder can, and rebuilt across the same six
multipliers it is monotone. Over the 60 ladder runs sitting at a profile median depth of 3 or less,
the median estimate falls from 0.938 of the truth at a multiplier of 1, to 0.901 at 2, to 0.853 at
8. Across all 96 engaged ladder runs the same slide is 0.927 to 0.904 to 0.865. The shallow side is
both better sampled and worse harmed, so the multiplier stays.

**The threshold and the multiplier are independent.** Every run in the low-depth ladder scores at
least 58, so no candidate threshold would stop reaching them. The threshold's real effect is how
much of the benchmark is disturbed to fix the few runs that need it: 186 runs at 16, against 73 at
24 for the same rescues. It was left at 16 because raising it also gives up the borderline runs
normalization lifts into the band, which is the half of the trade a count of regressions does not
show.

## Against the published run

The estimates here are not the paper's, because the code has moved since: #31 alone changed which
overlaps count. Holding normalization off isolates that drift, and it is close to neutral. Of the
3,370 runs, 9 reproduce the published estimate exactly and 273 land within 0.1% of it, but in
aggregate the published run puts 1,989 within 10% of the truth against 1,978 at HEAD, with 25 runs
under half their true size either way. The two sets of 25 overlap in 24: `SRR26465563` was published
at 0.4917 and sits at 0.5367 now.

Turning normalization on is what moves the aggregate, to 2,006 within 10% and 13 under half.

## The runs that were failing

`issue36_failure_outcomes.tsv` has one row per run that estimates under half its true size in any of
the three columns, 27 in all. Of the 25 failing at HEAD without normalization:

- 11 come back to within 10% of the truth, the largest being `SRR26465526` at 0.0122 to 0.9928, and
  `SRR10353548` at 0.0167 to 1.0274.
- 2 more clear 0.5x without reaching the band.
- 7 improve but stay under 0.5x, `SRR30162149` furthest at 0.115 to 0.343.
- 5 are unchanged, three of them because they score under the threshold and are never normalized:
  `SRR13183064` and `SRR13183067` both score 4.0, `SRR30357565` scores 15.33.
- 1 is made worse. `SRR13009132` goes from 0.6607 to 0.4813; normalization keeps 20.8% of its reads.

Thirteen still fail, and depth normalization is not the mechanism for them. Nine keep more than 90%
of their reads through normalization, so there is almost nothing to remove, and raven sizes eight of
those nine within 5% from the same reads. Two, `SRR13183064` and `SRR13183067`, are thin enough at a
median depth of 2 that genomescope and raven miss them as badly.

## The next mechanism, and a limit on this fit

Nine of the thirteen keep more than 90% of their reads, so depth is not their problem. The
estimator divides by the number of overlaps a query read finds, so an overlap that is not really an
overlap drives the estimate down, and an internal match, two reads sharing a repeat rather than a
locus, is exactly that. `-F` / `--filter-contained` drops them. That was tested on the 27
accessions still on disk from the #29 work, in four arms: normalization on or off, crossed with
`-F` on or off, stock binary throughout. The runs are in `issue36_filter_contained.tsv`.

`-F` is the mechanism the residual failures need. Of the 12 residuals in this set, `-F` takes 11 out
of the sub-0.5x band and 7 of them to within 10% of the truth: `SRR16631313` 0.121 to 0.901,
`SRR30357568` 0.220 to 0.998, `SRR30162149` 0.347 to 0.872. Only `SRR13183064` stays under half, and
it is one of the two that genomescope and raven also miss.

It cannot simply be turned on, because it raises every estimate. The per-run effect ranges from
1.08x to 7.46x and is never below 1, which is what removing overlaps from a denominator does. On the
13 accessions here that normalization alone already put within 10% of the truth, `-F` pushes all 13
above 1.1x, by 8 to 24 percent. So the filter is not miscalibrated in some runs and right in others:
it shifts the whole scale, and the estimator's constants were derived without it. Making it the
default is a recalibration, which is #38's subject, not a flag flip.

### Whether a run needs the filter can be read off one unfiltered pass

`-F` cannot be the default because it moves every estimate, but it does not have to be all or
nothing. Whether a run needs it is measurable before deciding, and cheaply, because `is_internal`
is arithmetic on fields the mapping already carries. A throwaway build counted what `-F` would have
discarded on runs that did not discard it (`filter/probe.patch` in the harness; estimates verified
unchanged against the stock binary). The statistic that matters is the fraction of a run's unique
overlaps that internal matches account for, `overlap_drop_frac` in
`issue36_internal_match_probe.tsv`.

It separates the two populations. Every run that normalization alone already put inside the band
sits at 0.61 or below, and six of the seven the filter rescues sit at 0.78 or above, reaching 0.945
on `SRR30162149`. The seventh, `SRR24489322`, sits at 0.424, in among the quiet runs, and no
threshold can take it without taking correct runs too. Across the whole of the gap between the two
populations the outcome is the same, so the constant first went in the middle of it, at 0.7. The
benchmark later moved it, for reasons this set could not have shown: see
[What the whole benchmark said about that rule](#what-the-whole-benchmark-said-about-that-rule).

## Switching the filter on from that measurement

The rule was then implemented rather than left as a prototype. `-F` takes a mode, as `--normalize`
does: `never` (the default, and what LRGE did before), `always` (what bare `-F` has always meant),
and `auto`. Under `auto` the mapping pass carries two overlap counts per query read, one with the
internal matches and one without, and picks between them once it has seen how much of the evidence
they account for. Both counts come out of the one pass, so nothing is mapped twice.

On the 27, that reproduced the prototype exactly: `auto` engaged on the eight runs the probe said it
would, took 19 of the 27 within 10% of the truth against 13 for never filtering and 9 for always,
and disturbed none of the thirteen already inside the band. `never` and `always` reproduced the
stock binary's estimate to the base pair on all 27, which is the check that mattered, because the
counting loops had been rewritten to carry both counts. Those runs are in
`issue36_dynamic_filter.tsv`.

## What the whole benchmark said about that rule

The 27 were the wrong 27, and the benchmark says so. All 3,370 accessions were rebuilt and run in
five arms: the last release, this branch with both mechanisms off, `--normalize auto`, and that plus
`-F auto` and `-F always`. The per-accession rows are in
`issue36_filter_benchmark_summary.tsv`.

`auto` does not fire too often. At the prototype's threshold it fired on 57 of 3,370, and at the
threshold it now carries, on 30. The problem is which 30.

| where the run started | fires at 0.7 | fires at 0.8 |
|---|---|---|
| under 0.5x | 9 | 7 |
| 0.5 to 0.9x | 3 | 1 |
| within 0.9 to 1.1x | 4 | 0 |
| 1.1 to 1.5x | 7 | 2 |
| over 1.5x | 34 | 20 |

Most of what it fires on was already reading high, and filtering can only push an estimate up, so
most of what it fires on gets worse. At 0.7 the median run it touches goes from 1.755x of the truth
to 3.413x, and 44 of the 57 end above 1.5x.

The reason is that the statistic does not mark underestimation. It marks repeats, and over a
representative sample the two point opposite ways:

| internal-match share | median estimate / truth | within 10% |
|---|---|---|
| 0.03 to 0.13 | 0.995 | 319 / 337 |
| 0.17 to 0.20 | 1.032 | 214 / 337 |
| 0.27 to 0.30 | 1.060 | 199 / 337 |
| 0.40 to 0.48 | 1.134 | 133 / 337 |
| 0.48 to 0.98 | 1.245 | 82 / 337 |

On the 27 the association ran the other way, because those accessions were selected for being the
failures depth normalization leaves behind. Repeat-rich genomes that read high were not in the set
to be seen. That is the whole of the error, and it is worth naming plainly: the prototype's
threshold was fitted on a sample chosen by the outcome it was meant to predict.

## Where the threshold went, and why the mode is not the default

Scoring every threshold against the same reads:

| threshold | fires | within 10% | under 0.5x | mean \|log2\| | taken out of the band | brought in |
|---|---|---|---|---|---|---|
| no filter | 0 | 2006 | 13 | 0.2195 | - | - |
| 0.70 | 57 | 2010 | 4 | 0.2357 | 4 | 8 |
| 0.76 | 41 | 2012 | 5 | 0.2324 | 1 | 7 |
| **0.80** | **30** | **2012** | **6** | **0.2286** | **0** | **6** |
| 0.85 | 17 | 2008 | 9 | 0.2264 | 0 | 2 |
| 0.90 | 12 | 2008 | 9 | 0.2237 | 0 | 2 |

No threshold beats not filtering on mean error. That is the finding, and it is why `never` stays the
default: a mechanism that raises the average error of a benchmark cannot be switched on for
everybody, however many catastrophes it fixes.

The threshold moved to 0.8 because that is the first value that disturbs nothing already correct.
The highest share among runs the estimator already puts within 10% of the truth is 0.796, on
`SRR13170267`, which filtering would take from 1.095x to 4.169x. Above 0.796 the count landing
within 10% also peaks. Below it the rule starts trading correct runs for rescues, which is a trade
an opt-in mode should not be making on the user's behalf.

What `auto` is for, then, is the failure it was found in: a genome you have reason to think is
repeat-rich, estimating far too low. It takes the runs under half their true size from 13 to 6 and
costs no run that was right. On the 27 outlier accessions the higher threshold gives 18 within 10%
and 5 under half, against 19 and 4 at 0.7; the two runs that trade places are the price of not
disturbing four correct benchmark runs.

## What the modes cost

Five arms on the same reads, per-run ratios rather than ratios of medians, in
`issue36_mode_cost.tsv`.

| | median | p90 | p99 | max | faster than 0.3.0 |
|---|---|---|---|---|---|
| this branch, both modes off | 1.00x | 1.03x | 1.36x | 3.62x | 1778 / 3370 |
| `--normalize auto` | 1.02x | 1.10x | 1.60x | 3.01x | 911 / 3370 |
| `--normalize auto -F auto` | 1.02x | 1.10x | 1.61x | 2.71x | 949 / 3370 |
| `--normalize auto -F always` | 1.01x | 1.10x | 1.61x | 2.65x | 996 / 3370 |

Three things come out of that. Everything added since the release costs nothing when the modes are
off, at a median of 0.998x. `-F auto` costs nothing on top of `--normalize auto`, at a median of
1.000x wall clock and 1.000x peak memory, because the second overlap count rides along with the
first. And normalizing shortens the tail rather than lengthening it: the slowest benchmark run takes
943 seconds unnormalized and 452 normalized, and the p99 falls from 190 to 163 seconds.

Cut by the condition each run was in, against the same binary with the mode off:

| condition | n | wall | peak memory | faster |
|---|---|---|---|---|
| detector does not fire | 3184 | 1.02x | 1.02x | 671 |
| detector fires | 186 | 1.09x | 1.06x | 74 |
| filter fires | 30 | 1.00x | 1.00x | 20 |

This corrects something the README claimed from a seven-accession sample: that a run the detector
fires on is usually quicker. Over 186 such runs the median is 1.09x and only 74 come out faster. The
spread is wide in both directions, 0.06x to 3.01x, because whether normalization saves time depends
on how many reads it removes.

### A limit this puts on the fit above

The threshold and the multiplier were swept on the default path, without `-F`, because that is the
invocation the published benchmark used. The plateau is a plateau for that path. If `-F auto` ever
becomes the default, both constants have to be swept again.

## What separates a rescue from a nuisance

The issue asked whether anything distinguishes the runs normalization rescues from the
already-correct
runs it merely moves, and named retained fraction, skew score, the shift in selected read lengths
and
the size of the enriched element as candidates. Taking the 13 rescued against the 13 pushed out of
the
band, the answer is that none of them separates cleanly.

| | rescued (13) | regressed (13) |
|---|---|---|
| median skew score | 180 | 22 |
| skew score range | 28 to 426 | 16 to 117 |
| median retained fraction | 0.335 | 0.447 |
| retained fraction range | 0.140 to 0.699 | 0.159 to 0.696 |
| median target read length shift | 1.199 | 0.997 |
| median estimate before | 0.027 | 1.044 |

Retained fraction is useless: the two ranges are almost the same interval. The skew score is the
best
of them and still only works at the ends. Every run scoring under 28 is a regression and every run
over 117 is a rescue, but between those the two interleave, four rescues against five regressions.

`SRR12247681`, the control the issue points at, is why. It scores 117, its target reads come out 25%
longer, and normalization drops 55% of them, which is the profile of a rescue on every axis the
issue
proposed. It was already right at 0.9398, and normalization takes it to 1.1101. The only thing that
tells it apart from a genuine rescue is where its estimate already was, which is the one thing the
estimator cannot know. So a second gate on the skew score would not buy anything here.

## Cost, on the sample this was first measured on

Superseded by [What the modes cost](#what-the-modes-cost), which measures the same thing on all
3,370 accessions. Kept because it is where the README's original figures came from, and because one
of them turned out to be wrong. Twenty-three accessions, three replicates each, stock binary: where
the detector does not fire the median cost is 1.01x wall clock and 1.02x peak memory, close to the
1.02x and 1.02x the benchmark gives. Where it does fire, the seven engaged accessions here ranged
from 0.32x to 2.74x and mostly came out faster, which is what the README was written from. Over 186
engaged runs the median is 1.09x and only 74 are faster, so "usually quicker" was an artefact of
seven runs.

## Reproducing it

The harness is `/scratch/user/uqmhal11/lrge-issue36`: a Snakemake workflow whose `download` rule is
the paper's, an `estimate` rule making the seven runs per accession, and `analyse.py`, `select.py`
and `detail.py` over the result. Reads are `temp()`, so the benchmark never sits on disk at once.

The filter work runs off the 27 accessions kept on disk from #29 instead, so it needs no downloads:
`filter/sweep.sh` for the four stock arms, `filter/probe.sh` and `filter/probe.patch` for the
instrumented pass, and `dynamic/sweep.sh` with `dynamic/summarise.py` for the three arms of the
shipped binary.

The five-arm benchmark is a second harness, `/scratch/user/uqmhal11/lrge-filter`, with the same
download rule and `workflow/analyse.py` over the result. Its `v030` arm is built from the
`lrge-0.3.0` tag. Only four of its five arms are measured: `auto` is derived from the unfiltered and
filtered arms and the share, which is exact because a run under `auto` takes one estimate or the
other, and was checked against a measured `-F auto` arm on all 3,370 accessions. That is what lets
the threshold move without the benchmark being rebuilt.

[#36]: https://github.com/mbhall88/lrge/issues/36
