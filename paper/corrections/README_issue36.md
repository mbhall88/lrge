# Fitting the depth-normalization constants on the benchmark

Issue [#36]. Two constants in `liblrge/src/depth_skew.rs` were picked by argument rather than
measurement: `SKEW_THRESHOLD`, how far the 99.9th percentile minimizer count has to run ahead of
the median before an input is called skewed, and `RETENTION_TARGET_MULTIPLIER`, the multiple of the
input's median depth that reads are kept down to. This is how they were fitted, and what the
benchmark said.

The short answer: neither moves. Normalizing at all is a large, reproducible win, and where the
constants sit inside a wide plateau is not.

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

It separates the two populations. The eight runs where `-F` gains most all sit above 0.77, reaching
0.945 on `SRR30162149`; every run that normalization alone already put inside the band sits at 0.61
or below. Switching `-F` on above a threshold in that gap beats both fixed choices on these 27
accessions:

| rule | within 0.9-1.1x | under 0.5x |
|---|---|---|
| `auto`, never `-F` (ships today) | 13 | 12 |
| `auto`, always `-F` | 9 | 1 |
| `auto`, `-F` above a 0.62 to 0.75 drop fraction | 19 | 4 |

No run that was already inside the band gets filtered at any threshold from 0.62 up, which is what
keeps the overshoot away, and the result is flat across that range rather than balanced on a point.
The four runs still failing under the rule are the ones that need `-F` but do not look like it:
`SRR10259778` gains 3.79x from filtering on a drop fraction of only 0.26.

This is a prototype on 27 outlier-enriched accessions, not a fitted constant. It says the signal
exists and is cheap, not where the threshold belongs; that needs the full benchmark, and the
benchmark's reads would have to be rebuilt again to get it.

This is also the limit on the fit above. The threshold and the multiplier were swept on the default
path, without `-F`, because that is the invocation the published benchmark used. The plateau is a
plateau for that path. If `-F` ever becomes the default, both constants have to be swept again.

Two caveats on this table specifically. The 27 accessions are the outlier-enriched set collected for
#29, not a sample of the benchmark, so the counts above describe them and not LRGE in general. And 8
of the 27 hit the 1 Gbp cap, where the read set on disk is not quite the one the benchmark rebuilt,
so their unnormalized values differ slightly from `issue36_benchmark_sweep.tsv`.

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

## Cost

Measured on 23 accessions, three replicates each, stock binary. Where the detector does not fire the
median cost is 1.01x wall clock and 1.02x peak memory, so detection is close to free, which is the
work of #46, #47, #52 and #53. Where it does fire the run is usually *faster*, because
normalization leaves fewer reads to overlap: the seven engaged accessions range from 0.32x to
2.74x.

## Reproducing it

The harness is `/scratch/user/uqmhal11/lrge-issue36`: a Snakemake workflow whose `download` rule is
the paper's, an `estimate` rule making the seven runs per accession, and `analyse.py`, `select.py`
and `detail.py` over the result. Reads are `temp()`, so the benchmark never sits on disk at once.

[#36]: https://github.com/mbhall88/lrge/issues/36
