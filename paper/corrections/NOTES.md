# Post-publication corrections: working notes

**Status: working notes, not the post.** This file accumulates everything needed to write a
technical update for lrge users once both mechanisms behind [#29][i29] are addressed. Mechanism 2 is
fixed and measured. Mechanism 1 is implemented in [PR #45][pr45] and measured on 17 runs (§4.6), but
that PR is unmerged and unreleased, and nothing yet bounds how many of the 3,370 benchmark runs its
default would change. That gap is the main thing still blocking the post.

Paper: Hall MB, Zhou C, Coin LJM. *Genome size estimation from long read overlaps.* Bioinformatics
41(11), November 2025, btaf593. [doi:10.1093/bioinformatics/btaf593][doi]

No republication is planned. The intended output is a standalone technical post explaining what
changed, why, and which published results it affects.

---

## 1. The headline, in one paragraph

Two independent defects made lrge under-estimate genome size on a small minority of runs. They have
different causes, different fixes, and — crucially — **different consequences for the published
results**. Mechanism 2 (an inverted overlap filter) never touched the default code path, so it
changes no number in the paper; it made an opt-in flag do the opposite of what it documented, and it
explains 10 of the benchmark's 25 worst failures. Mechanism 1 (depth skew) *does* affect default
output and *will* change published numbers once fixed.

## 2. How large is the problem, in benchmark terms

From `paper/results/estimates/estimates.tsv`, two-set strategy, relative to truth:

| Band | Runs | Share |
|---|---|---|
| < 0.5x | 25 | 0.74% |
| < 0.8x | 46 | 1.36% |
| 0.8–1.2x | 2,534 | 75.19% |
| > 1.2x | 790 | 23.44% |
| **Total** | **3,370** | |

**All 25 have now been re-run with the mechanism-2 fix.** They do not split by platform — they split
by whether `-F` does anything at all:

| Class | Runs | `-F` response | What they are |
|---|---|---|---|
| **Responsive** | 18 | 1.4x–11.3x improvement | all 10 PacBio, plus 8 ONT |
| **Inert** | 7 | no change (<=1%) | *N. gonorrhoeae*, all collapsed to 0.012x–0.027x |

The inert class looked at first like a third failure mode. It is not: all seven are heavily
depth-skewed by high-copy plasmids, and `-F` is blind to them only because their reads are shorter
than the plasmid, which makes the bad overlaps look like dovetails rather than internal matches
(§4.5). So the tail is fully attributed — mechanism 1 and mechanism 2, with `-F` incidentally
rescuing part of the mechanism-1 population.

The mechanism-1 fix has since been benchmarked on all 15 ONT sub-0.5x runs (§4.6). Eleven land in
0.8x–1.2x where none did before, including all seven inert runs. Of the four that do not, three are
internal-match failures that normalization is not aimed at, and one was missed by the skew detector.

---

## 3. Mechanism 2 — inverted internal-match filtering (fixed)

Issue [#31][i31]. PR [#40][pr40], merged 2026-08-31 as `fc09157` (`fix!: correct inverted
internal-match filtering`) + `7af167a` (`refactor: single definition of the default overhang ratio`).

### 3.1 What was wrong

Three defects, all in the handling of `-F`/`--filter-contained`:

1. **The two two-set code paths filtered in opposite directions.** `align_reads` (the default path)
   and `ava` called `PafRecord::is_internal()`, which ended `overhang_ratio < max_overhang_ratio`,
   and the callers skipped the mapping when it returned true — discarding the *small*-overhang
   alignments. `align_reads_inverse` (reached via `--use-min-ref`) inlined the identical formula
   with the opposite comparison, discarding the *large*-overhang alignments. Same inputs, same
   formula, inverted comparison: `-F` meant opposite things depending on which path a run took.
2. **The default path was the wrong one.** Under the miniasm overhang formula, a dovetail overlap
   and a containment both leave near-zero overhang; an internal match — two reads sharing a repeat
   rather than the same locus — leaves large unaligned tails. So internal matches are the
   large-overhang case. The default path discarded the genuine overlaps and kept the internal
   matches, contradicting the flag's own help text.
3. **`--max-overhang-ratio` was silently a no-op without `-F`.** Both builders assigned the ratio
   only inside `if filter_contained`. Passing it alone produced byte-identical output with no
   warning.

Worth stating plainly in the post: the flag is named `--filter-contained` but documented as
filtering internal matches, and **the predicate cannot distinguish the two** — containment and
dovetail both give near-zero overhang. The fix implements the documented behaviour (internal-match
removal). Contained-read removal would need a different test.

### 3.2 What changed, and what did not

- `is_internal` now returns true on *large* overhang, and there is one implementation shared by all
  three call sites.
- `--max-overhang-ratio` now requires `-F` and errors rather than being ignored.
- **`-F` is off by default, so default output is unchanged.** Verified empirically: five accessions
  (SRR16767125, SRR26465560, SRR8618952, SRR26465526, SRR12247681) reproduce their pre-fix estimates
  exactly. **This is why no figure or table in the paper changes because of mechanism 2.**

**Release status: unreleased.** The fix is on `main`, 7 commits ahead of tag `lrge-0.3.0`. Because
it is a `fix!`, the next release is a breaking one. The post cannot tell users to upgrade until that
release exists — check before publishing.

### 3.3 Measured impact

All 10 PacBio sub-0.5x runs, re-run from the fixed binary. Prep mirrors
`paper/workflow/scripts/download.sh`; settings `-P pb -s 4556 -T 10000 -Q 5000`. Full data in
[`rerun_estimates.tsv`](./rerun_estimates.tsv).

| Run | Organism | Truth | Paper 2-set | Default (re-run) | `-F` | `-F` @0.05 | Target depth |
|---|---|---|---|---|---|---|---|
| SRR30357565 | *C. belfantii* | 2,791,411 | 0.249x | 0.238x | **0.992x** | 1.011x | 78.3x |
| SRR30357566 | *C. belfantii* | 2,834,287 | 0.404x | 0.391x | **1.003x** | 1.019x | 73.6x |
| SRR30357567 | *C. belfantii* | 2,897,507 | 0.239x | 0.253x | **1.005x** | 1.021x | 73.7x |
| SRR30357568 | *C. belfantii* | 2,759,957 | 0.207x | 0.273x | **1.001x** | 1.016x | 77.0x |
| SRR30162149 | *Bradyrhizobium* sp. 2S1 | 11,093,179 | 0.122x | 0.105x | **0.880x** | 0.928x | 8.6x |
| SRR17188519 | *X. oryzae* | 4,929,712 | 0.118x | 0.107x | 0.900x | 0.943x | 6.2x |
| SRR16767125 | *X. oryzae* | 4,963,479 | 0.051x | 0.061x | 0.695x | 0.730x | 0.70x |
| SRR16631313 | *X. oryzae* | 4,966,657 | 0.055x | 0.052x | 0.528x | 0.570x | 0.48x |
| SRR13183064 | *B. dolosa* | 6,409,097 | 0.323x | 0.302x | 0.434x | 1.288x | 2.7x |
| SRR13183067 | *C. violaceum* | 4,808,040 | 0.400x | 0.380x | 0.629x | 1.584x | 2.4x |

Re-run default figures differ from the paper's by up to ~30% (worst: SRR30357568, 0.207x vs 0.273x)
because inputs were re-downloaded under rasusa 4.x, whose sampling differs from the 2.x run used for
the paper. These are collapsed estimates and are unstable to resampling; the point of the column is
that both are far below truth, not that they agree. Do **not** cite the paper-vs-re-run difference as
evidence about the fix — the rigorous test of "default is unchanged" is §3.2, which compares two
binaries on *identical* input and gets exact equality.

#### Direct before/after on identical input

The table above compares the fixed binary against the *paper's* numbers, which also differ by
resampling. The cleanest demonstration is the same input run through both binaries — `SRR16767125`,
where the pre-fix runs from the original #31 investigation were kept:

| Variant | Pre-fix | Infinite est. | Post-fix | Infinite est. |
|---|---|---|---|---|
| default | 304,556 | 428 | 304,556 | 428 |
| `-F` | 197,672 (0.040x) | 1,883 | 3,447,360 (0.695x) | 997 |
| `-F --max-overhang-ratio 0.05` | 198,720 | 1,864 | 3,625,604 | 1,098 |
| `--max-overhang-ratio 0.05` alone | 304,556 | 428 | *(now rejected by the CLI)* | — |

Three things are visible at once: the default column is unchanged; `-F` went from making the estimate
**5x worse** to making it **11x better**; and the last row is defect 3 — passing the ratio without
`-F` produced byte-identical output to the default run.

#### Threshold sweep

Post-fix, on `SRR16767125` (0.70x target depth, truth 4,963,479):

| `--max-overhang-ratio` | Estimate | rel |
|---|---|---|
| 1.0 | 3,030,010 | 0.610x |
| 0.5 | 3,260,360 | 0.657x |
| 0.2 (default) | 3,447,360 | 0.695x |
| 0.05 | 3,625,604 | 0.730x |

Monotone: stricter filtering raises the estimate on a thin-target run. Read alongside §3.6 — the same
sweep on a well-sampled run moves the estimate the *other* side of truth.

**The strongest single result for the post** is the *C. belfantii* group: four near-replicates (same
species, same BioProject, same platform, assemblies within 5% of each other) whose published
estimates spread 0.207x–0.404x — a factor of two apart on effectively the same input — collapse to
0.992x–1.005x. That spread was not biology and not sampling noise.

### 3.4 Recovery tracks target-set *depth*

The quantity that governs how completely `-F` recovers an estimate is target-set depth:
`total target bases / genome size`, logged at `-vv` as `Total target bases`.

| Run | Targets | Target depth | `-F` | Infinite estimates (of 5,000) |
|---|---|---|---|---|
| SRR30357565 | 10,000 | 78.3x | 0.992x | 0 |
| SRR30357568 | 10,000 | 77.0x | 1.001x | 0 |
| SRR30357567 | 10,000 | 73.7x | 1.005x | 0 |
| SRR30357566 | 10,000 | 73.6x | 1.003x | 0 |
| SRR30162149 | 10,000 | 8.6x | 0.880x | 3 |
| SRR17188519 | 2,992 | 6.2x | 0.900x | 0 |
| SRR13183064 | 10,000 | 2.7x | 0.434x | 1,745 |
| SRR13183067 | 10,000 | 2.4x | 0.629x | 2,488 |
| SRR16767125 | 338 | 0.70x | 0.695x | 997 |
| SRR16631313 | 238 | 0.48x | 0.528x | 1,674 |

Two regimes:

- **Well-sampled target set (>= ~70x): `-F` is exact.** Four runs within 0.8% of truth, zero
  infinite estimates.
- **Thin target set: `-F` removes most of the deficit but under-estimates in proportion.**

Note the residual bias is *not* only lost query reads. SRR17188519 has **zero** infinite estimates
at 6.2x depth and still lands at 0.900x; SRR30162149 has 3 at 8.6x and lands at 0.880x. Losing
queries to infinity makes it worse, but thin targets bias downward even when every query gets an
estimate.

The three *X. oryzae* runs looked like poor recoveries until this was worked out: they are small read
sets (5,238–7,992 reads), so after the 5,000-read query split only 238–2,992 targets remained. An
earlier note on [#29][i29] attributed recovery to target *count*; count was a proxy for depth.

### 3.5 `-F` should not become the default

On `SRR8618952` — a healthy *X. oryzae* run, 14.9x target depth — `-F` moves the estimate from
1.200x to 1.467x. Turning it on unconditionally would trade a rare under-estimate for a common
over-estimate.

**The ONT batch makes this sharper: `-F` overshoots on runs it is supposed to rescue.** Three of the
13 land above truth (`SRR26465563` 1.380x, `SRR26715165` 1.294x, `DRR213976` 1.132x), and at
`--max-overhang-ratio 0.05` five do, two of them past 2x (`SRR26465563` 2.220x, `SRR26715165`
2.040x). So `-F` is not a safe default even restricted to the runs that respond to it: on this
evidence it converts roughly a third of its successes into over-estimates.

**The read-count control.** Before the mechanism was identified, small read count was the leading
hypothesis, since every *X. oryzae* failure had one. `SRR8618952` was subsampled to test it directly
(truth 4,907,789, default flags):

| Reads | Estimate | rel |
|---|---|---|
| 5,338 | 3,178,520 | 0.648x |
| 8,000 | 5,534,898 | 1.128x |
| 15,000 | 5,818,243 | 1.186x |
| 23,000 | 5,843,991 | 1.191x |
| 40,000 | 5,768,675 | 1.175x |
| 136,128 (all) | 5,890,324 | 1.200x |

At 5,338 reads — matching `SRR16767125` exactly — a healthy run still returns 0.648x, not 0.061x. So
a small read set degrades an estimate but does not collapse it, which is what ruled the hypothesis
out and sent the investigation to the filter. It also shows the degradation is real, which is the
§3.4 target-depth effect seen from the other direction.

### 3.6 A fixed `--max-overhang-ratio` is wrong in both directions

The stricter 0.05 threshold is harmless at high depth (1.011x–1.021x) and overshoots badly at low
depth (**1.288x**, **1.584x**). Same flag, same value, error changing sign across the depth range.
This is the empirical case for a threshold that adapts to achieved overlap density — the same signal
the mechanism-1 design already computes. Folding it into that work is the current recommendation;
it is not yet a ticket.

### 3.7 Weak evidence to flag honestly

`SRR13183064` and `SRR13183067` are FDAARGOS raw subread sets: 58 bp minimum read length, ~1.1–1.7 kb
mean, published median lengths 696 bp and 373 bp. A third to a half of their query reads produce no
estimate at all, before and after the fix. Overlap-based estimation on those inputs is marginal
regardless of filtering. Treat their partial recovery as "the mechanism is present", not as
calibration.

---

## 4. Mechanism 1 — depth skew (implemented in [PR #45][pr45], unmerged)

Issue [#29][i29] proper. Ticket chain: [#33][i33] (detect skew, observation only) → [#34][i34] (the
fix) → [#35][i35] (low-memory fallback) → [#36][i36] (fit constants) → [#37][i37] (release) →
[#38][i38] (re-derive quantiles). [#32][i32] and [#39][i39] are closed.

### 4.1 What is wrong

The median of per-read estimates collapses when a subpopulation of query reads comes from
high-depth sequence. `SRR26465560` (*E. faecium*, the run #29 was reported against) returns 315 kb
against a truth of 3.20 Mb. Six plasmids totalling 81 kb sit at up to 167x chromosomal depth and
supply the majority of the reads, so the estimator returns something closer to the size of the
over-represented element than the genome.

**Caveat added 2026-09-01:** this run is not a clean example of the mechanism. Most of its error
turns out to be mechanism 2 — see §4.3. The depth-skew diagnosis stands as a description of what
happens under skew; it is no longer safe to use this particular run as the demonstration of it.

### 4.2 The planned fix

Depth-aware read selection. Each read's depth is estimated as the **median** of its minimizer
counts — a median suppresses reads lying wholly inside a high-depth element while sparing reads that
merely span one, and spanning reads carry the most genome-size information. Reads are retained with
probability `min(1, C/depth)`. Target and query reads are drawn from a single normalized pool.
Normalization and subsampling happen together via weighted reservoir sampling, so the pass count is
unchanged. `--normalize auto|always|never`, defaulting to `auto`; an unskewed run takes the existing
path and produces byte-identical output.

### 4.3 The two mechanisms overlap on the reported run

`-F` was run post-fix on the two most extreme ONT failures (truth-relative, default flags in the
first column):

| Run | Organism | Truth | Target depth | Default | `-F` | `-F` @0.05 |
|---|---|---|---|---|---|---|
| SRR26465560 | *E. faecium* | 3,198,172 | 17.6x | 0.098x | **0.643x** | **0.855x** |
| SRR26465526 | *N. gonorrhoeae* | 2,230,369 | 25.1x | 0.012x | 0.013x | 0.023x |

**`SRR26465560` is the run #29 was reported against and the run [#34][i34] is specified around, and
most of its error is mechanism 2, not depth skew.** It recovers from 0.098x to 0.643x with no
depth-aware selection at all, and to 0.855x at the stricter threshold. Its residual gap is consistent
with genuine depth skew on top — at 17.6x target depth the §3.4 relationship predicts roughly 0.9x
from mechanism 2 alone, and it reaches 0.643x — but the plasmid story is no longer the whole story.

`SRR26465526` is the opposite: 25x target depth, and `-F` moves it from 0.012x to 0.013x. Whatever
collapses that run, mechanism 2 is not it.

Consequences to carry into the post and the tickets:

- The #34 acceptance criterion ("`SRR26465560` estimates substantially closer to 3.20 Mb than its
  current 0.098x") is still literally correct, because #31 did not change the default path. But it no
  longer isolates what #34 is supposed to fix: measure #34 against the post-#31 `-F` figure, or pick
  a run like `SRR26465526` where mechanism 2 demonstrably does nothing.
- Done 2026-09-01: `-F` was run on the remaining 13 ONT sub-0.5x runs, which is what made the
  attribution in §4.4 and §4.5 possible.

### 4.4 Hypothesis: mechanism 1 *generates* mechanism 2

Raised by MBH, 2026-09-01. The two mechanisms may not be independent. If a region is enriched in
depth, a large share of the read set comes from that region; those reads overlap each other heavily,
and if the enriched region is a repeat or a multi-copy element, many of those overlaps will be
internal matches rather than genuine dovetails. On that account mechanism 2 is partly a *symptom* of
mechanism 1, not a separate defect that happens to co-occur.

This fits the only skewed run measured so far. `SRR26465560` carries six plasmids at up to 167x
chromosomal depth, and `-F` alone recovers it from 0.098x to 0.643x (§4.3) — which is what you would
expect if the plasmid-derived reads were contributing a large population of internal matches that the
filter now removes.

**Predictions, in rough order of how cheaply they can be tested:**

1. Skewed runs should show a higher internal-match fraction in their overlap set than unskewed runs
   of comparable depth. **Partially tested, and the direction holds.** Measured on the 13 ONT runs
   (`[PAF]` lines, `internal_frac_0.2` in the TSV), internal-match fraction predicts whether `-F`
   helps: every run at >=25% responds (gain 2.1x–4.8x), and 6 of the 7 below 21% are completely inert
   (§4.5). The one exception is SRR10388020 at 9.9%, which still gains 1.5x.

   | Run | Internal frac. | `-F` gain |
   |---|---|---|
   | DRR213976 | 79.9% | 2.7x |
   | SRR26715165 | 47.7% | 3.2x |
   | SRR24489322 | 47.1% | 2.1x |
   | SRR26465563 | 44.1% | 2.6x |
   | SRR10259778 | 34.7% | 4.2x |
   | SRR26715166 | 25.2% | 4.8x |
   | SRR26465523 | 20.6% | 1.0x |
   | SRR26465524 | 19.5% | 1.0x |
   | SRR26465521 | 18.8% | 1.0x |
   | SRR10861751 | 12.3% | 1.0x |
   | SRR10861747 | 12.3% | 1.0x |
   | SRR10388020 | 9.9% | 1.5x |
   | SRR10353548 | 9.2% | 1.0x |

   **Now tested against depth (§4.5), and the result refines the hypothesis rather than confirming
   it.** All ten *N. gonorrhoeae* runs are heavily skewed by high-copy plasmids, yet their
   internal-match fractions range 9%–48%. So skew magnitude does **not** predict internal-match
   fraction — if anything the most extremely skewed runs have the *fewest* internal matches.

   What decides it is read length against the size of the enriched element. The coupling is real —
   depth skew is what generates the bad overlaps — but it does not have a fixed signature: the same
   skew presents as dovetails or as internal matches depending on that ratio, and only the latter is
   visible to `-F`. Internal fraction was not measured on the PacBio runs; adding it means re-running
   them with `--keep-temp`.
2. `-F` should help skewed runs more than unskewed ones. Partially supported already: it recovers
   `SRR26465560` substantially while *degrading* the healthy `SRR8618952` (§3.5).
3. After #34 lands, the residual benefit of `-F` on a skewed run should **shrink**, because
   normalization removes the reads that were generating the internal matches. If `-F` still helps
   just as much post-normalization, the mechanisms are independent after all. **Tested (§4.6):
   it shrinks on skew-driven failures and holds on internal-match-driven ones.** `SRR26465560`
   6.53x -> 1.28x, `SRR26715165` 3.16x -> 1.20x, `SRR26465563` 2.57x -> 1.26x; but `SRR24489322`
   2.06x -> 2.04x and `DRR213976` 2.74x -> 2.21x.
4. The two fixes should be **sub-additive** on skewed runs. Worth stating as a risk rather than a
   nicety: `-F` already overshoots on healthy input (1.200x -> 1.467x, §3.5), so normalization plus
   filtering could over-correct a skewed run past truth. **Confirmed (§4.6): `auto -F` overshoots on
   14 of 15 outliers.** `-F` must not be switched on automatically under detected skew.

Prediction 3 decides the framing of the post, and §4.6 now answers it: **neither framing alone is
right.** For the eleven runs normalization fixes, the two mechanisms are coupled and `-F`'s residual
benefit collapses — one problem, two symptoms. For the three runs it does not fix, `-F`'s benefit is
untouched and the failure is internal matches with no meaningful skew contribution — a genuinely
separate bug. The post should say the mechanisms overlap on most of the tail without claiming they
are the same defect.

§4.5 tilts this toward the first reading for the ONT tail — every one of those failures is depth
skew, and `-F` rescues the subset where the geometry happens to expose it. It says nothing yet about
the PacBio tail, where no depth profiling has been done.

Note this also complicates §3.4: the target-depth relationship there was fitted on PacBio runs with
no known skew. Whether recovery-vs-depth looks the same on skewed input is unknown.

### 4.5 The inert runs are mechanism 1 after all — read length decides whether `-F` sees it

**This section supersedes an earlier draft that called these a third failure mode. They are not.**

Depth profiles were run on all ten *N. gonorrhoeae* sub-0.5x runs — the seven inert ones and, as a
same-species control, the three that respond to `-F`. Every one of the ten is heavily depth-skewed,
and in exactly the same way: small, very high-copy plasmids.

| Run | Class | Chromosome | Small plasmids | Large plasmid |
|---|---|---|---|---|
| SRR10353548 | inert | 2,219 kb @ 1x | 6 kb @ 257x, 5 kb @ 226x | 42 kb @ 2x |
| SRR26465526 | inert | 2,179 kb @ 1x | 6 kb @ 606x, 5 kb @ 333x | 43 kb @ 8x |
| SRR26465524 | inert | 2,236 kb @ 1x | 6 kb @ 338x, 5 kb @ 177x | 43 kb @ 17x |
| SRR10388020 | responsive | 2,167 kb @ 1x | 6 kb @ 246x, 5 kb @ 229x | 43 kb @ 5x |
| SRR26715165 | responsive | 2,228 kb @ 1x | 5 kb @ 63x | 40 kb @ 14x |
| SRR26715166 | responsive | 2,229 kb @ 1x | 5 kb @ 126x | 40 kb @ 20x |

Roughly half of all mapped depth sits in ~11 kb of plasmid. So the earlier reading in these notes —
"nothing looks wrong with the inputs" — was wrong: the inputs are as skewed as anything in the
benchmark. What was actually observed is that *the overlaps they produce do not look like internal
matches*, which is a different claim.

**Depth skew does not predict whether `-F` helps. Read length does.**

| Run | Class | Median read | Top plasmid | read/plasmid | Internal frac. | `-F` gain |
|---|---|---|---|---|---|---|
| SRR26715165 | responsive | 11,824 | 5 kb | **2.36** | 47.7% | 3.2x |
| SRR26715166 | responsive | 9,406 | 5 kb | **1.88** | 25.2% | 4.8x |
| SRR26465524 | inert | 6,867 | 6 kb | 1.14 | 19.5% | 1.0x |
| SRR26465521 | inert | 6,110 | 6 kb | 1.02 | 18.8% | 1.0x |
| SRR10861751 | inert | 4,986 | 5 kb | 1.00 | 12.3% | 1.0x |
| SRR10861747 | inert | 4,867 | 5 kb | 0.97 | 12.3% | 1.0x |
| SRR26465523 | inert | 5,563 | 6 kb | 0.93 | 20.6% | 1.0x |
| SRR26465526 | inert | 5,607 | 6 kb | 0.93 | — | 1.0x |
| SRR10353548 | inert | 3,928 | 6 kb | 0.65 | 9.2% | 1.0x |
| SRR10388020 | responsive | 3,687 | 6 kb | 0.61 | 9.9% | 1.5x |

Every run whose reads are shorter than ~1.2 plasmid-lengths is inert. Both runs whose reads are
~2 plasmid-lengths respond strongly, and they carry the highest internal-match fractions.
`SRR10388020` is the one that does not fit — short reads, low internal fraction, yet a 1.5x gain, the
weakest of the responders.

**The geometric account.** Reads piled on a small circular plasmid overlap each other. If a read is
*shorter* than the plasmid, two such reads align end to end with almost nothing hanging off — the
geometry of a genuine dovetail, near-zero overhang, invisible to any overhang threshold. If a read is
*longer* than the plasmid it wraps past the origin, so two such reads share only part of their length
and leave large unaligned flanks — an internal match, which `-F` removes. Same biological cause,
opposite overlap geometry, decided by read length against element size.

**Consequences:**

- The seven inert runs are **mechanism 1**, so #34 should fix them, and the benchmark tail does not
  contain an unexplained third mechanism. All 25 sub-0.5x runs are now attributed.
- `-F` rescuing a skewed run is **incidental**, not a fix for skew. It happens to catch the subset
  where read length exceeds the enriched element. That is worth saying plainly in the post, because
  the *X. oryzae* and *C. belfantii* recoveries otherwise invite the reading that `-F` addresses
  depth skew generally. It does not.
- These runs are the sharpest acceptance test for #34 available: `-F` provably does nothing to them,
  so any improvement is attributable to normalization alone. Better than `SRR26465560` (§4.3), whose
  error is now known to be substantially mechanism 2. **Run against [PR #45][pr45]: all seven recover,
  0.012x–0.027x to 0.92x–1.00x (§4.6).**

### 4.6 Impact — measured on 17 runs

[PR #45][pr45] at commit `254e146` was benchmarked on all 15 ONT sub-0.5x outliers and two controls,
at paper settings (seed 4556, 10,000 targets, 5,000 queries, 8 threads, correct platform preset).
One release binary was built and every job invoked it. 37 cells, no failures. Per-row detail is in
[`depth_normalization_estimates.tsv`](./depth_normalization_estimates.tsv); the `old` column below is
`rerun_default` from [`rerun_estimates.tsv`](./rerun_estimates.tsv), which is the same input through
the post-#31 binary at defaults.

`--normalize never` reproduced `rerun_default` exactly on all five runs it was checked against
(`SRR10861751`, `SRR12247681`, `SRR26465560`, `SRR26715165`, `SRR8618952`), so `old` and `auto` are
comparable without qualification.

| Run | Skew score | Kept | old | `auto` | `auto -F` |
|---|---|---|---|---|---|
| SRR10353548 | 280x | 49.7% | 0.017x | **1.003x** | 1.124x |
| SRR10388020 | 205x | 52.9% | 0.351x | **1.025x** | 1.137x |
| SRR26465523 | 142x | 47.0% | 0.027x | **0.980x** | 1.237x |
| SRR26465563 | 82x | 64.5% | 0.537x | **0.980x** | 1.236x |
| SRR10861747 | 340x | 36.1% | 0.019x | **0.970x** | 1.056x |
| SRR26465526 | 425x | 20.3% | 0.012x | **0.961x** | 1.200x |
| SRR26465521 | 252x | 28.9% | 0.019x | **0.937x** | 1.151x |
| SRR10861751 | 203x | 33.1% | 0.018x | **0.936x** | 1.034x |
| SRR26465560 | 165x | 55.2% | 0.098x | **0.932x** | 1.197x |
| SRR26715165 | 38x | 61.3% | 0.409x | **0.927x** | 1.111x |
| SRR26465524 | 180x | 29.5% | 0.023x | **0.924x** | 1.166x |
| DRR213976 | 32x | 96.0% | 0.414x | 0.498x | 1.102x |
| SRR24489322 | 21x | 91.5% | 0.205x | 0.269x | 0.550x |
| SRR26715166 | 12x | — | 0.215x | 0.215x | 1.036x |
| SRR10259778 | 39x | 95.1% | 0.105x | 0.162x | 1.688x |
| *SRR12247681* (control) | 104x | 54.7% | 0.937x | *1.113x* | — |
| *SRR8618952* (control) | 3x | — | 1.200x | *1.200x* | — |

Outliers landing in 0.8x–1.2x: 0 of 15 before, 11 of 15 after. Mean absolute relative error
across the 15 falls from 0.835 to 0.222.

**The seven inert *N. gonorrhoeae* runs of §4.5 all recover**, from 0.012x–0.027x to 0.92x–1.00x.
That was the acceptance test those notes nominated, on the grounds that `-F` provably does nothing to
them, so the improvement is attributable to normalization alone. It passes.

#### The four that do not recover

Three of them — `DRR213976`, `SRR24489322`, `SRR10259778` — carry the highest internal-match
fractions in the ONT set (79.9%, 47.1%, 34.7%). Their skew scores are also the lowest of the
engaging runs (32x, 21x, 39x), and normalization drops only 4%–9% of their reads. These are
mechanism-2 failures, and normalization is not aimed at them: `-F` alone took `DRR213976` to 1.132x
before any of this. Every ONT run below 21% internal fraction lands in 0.92x–1.03x after
normalization. The converse does not hold — `SRR26715165` (47.7%) and `SRR26465563` (44.1%) recover
fully — so internal fraction bounds what normalization can do rather than predicting it.

`SRR26715166` is the one implementation finding here. **The detector missed it.** It scores 12.00x
against the 16x threshold, computed from 69 sampled reads, because the 1% detection sample was
applied to a 7,473-read input. Forcing the work it declined recovers the
run: `--normalize always` retains 4,140 of 7,473 reads and returns 1,900,525 bp, or 0.837x, against
0.215x under `auto`. So normalization would have fixed it and detection stopped it. The run is also
depth-profiled in §4.5 (5 kb plasmid at 126x median, 39% of mapped bases above 10x median), which is
independent confirmation that it is skewed. This argues for a floor on the detection sample rather
than a change to the 16x threshold, and belongs to [#36][i36].

#### The controls

`SRR8618952` is untouched: skew score 3.00x, detector abstains, and the selected `target.fa` and
`query.fa` are **byte-identical** to `--normalize never` (md5 `7cb387cb…` and `6ee97d83…`). Its depth
profile is flat: highest 1 kb window 1.9x the median, nothing above 10x. That is the criterion
"an evenly-covered input produces output identical to the pre-change behaviour", met.

`SRR12247681` regressed. The detector engages at 104x, normalization retains 172,655 of 315,509
reads, and the estimate moves from 0.937x to 1.113x, so absolute error grows from 6.3% to 11.3%.
**This is not a false positive.** Its depth profile is as skewed as anything in the failing tail:
a 3 kb plasmid at 697x chromosomal depth and a 6 kb plasmid at 197x, against a 5,060 kb chromosome
at 125x, with 39.3% of mapped bases in windows above 10x median. The run was estimating well
*despite* severe skew.

The gap generalizes past this one run. The detector answers "is this input depth-skewed", and the
set of skewed runs is larger than the set of runs whose estimates were wrong. Every run in
the second set is in the first, but not the reverse, and normalization moves runs in both. The 3,345
benchmark runs that are already correct are exposed to that difference. **Seventeen runs cannot say
how many of the 3,370 would move.** This sample was picked for being broken, so it yields no
population rate. Getting one means running the detector across the benchmark, which is its own
experiment.

#### Cost

On `SRR8618952`, where the detector abstains, three alternating repeats of each mode give
6.99s / 7.05s / 7.21s for `never` against 71.18s / 72.23s / 72.25s for `auto` — means of 7.08s
and 71.89s, a 10.2x wall-clock cost, at 1.18 GiB versus 1.19 GiB mean peak RSS. All six returned the
same estimate, 5,890,324 bp. The modes were interleaved so page-cache warmth cannot favour either.
The cost is the minimizer depth profile, which `auto` builds over every read during the counting
pass before it knows whether it will need it; `never` skips profiling entirely. Where the overlap
stage is heavier the ratio falls. `SRR26465560` runs 62s against 90s and `SRR26715165` 50s against
89s, but those pairs came from different Slurm jobs and are indicative only. Peak RSS never rose
materially in any pair; on several skewed runs `auto` used *less* memory, because normalization
shrinks the pool. Across all 37 cells peak RSS ranged 0.71–4.63 GiB.

#### What this settles from §4.4

**Prediction 3 holds where it was meant to.** On the runs where `-F` used to deliver a large gain and
normalization now engages hard, the residual benefit of `-F` collapses: `SRR26465560` 6.53x -> 1.28x,
`SRR26715165` 3.16x -> 1.20x, `SRR26465563` 2.57x -> 1.26x. Normalization is removing the reads that
were generating those internal matches, which is what the coupling hypothesis predicted. On the three
mechanism-2-dominated runs it does not: `SRR24489322` 2.06x -> 2.04x, `DRR213976` 2.74x -> 2.21x,
`SRR10259778` 4.18x -> 10.42x. So the two mechanisms are coupled on skew-driven failures and
independent on internal-match-driven ones, which is the first reading of §4.4 for most of the tail
but not a clean one-problem story.

**Prediction 4 holds.** `auto -F` overshoots truth on 14 of 15 outliers, 11 of
them past 1.1x and four at or beyond 1.2x, topping out at 1.688x (`SRR10259778`). The single
exception is `SRR24489322` at 0.550x. Mean absolute error, 0.195, is barely better than `auto` alone
at 0.222, and the errors now point the same way. `-F` must not be switched on automatically under
detected skew.

#### Follow-up measurements, 2026-09-02

Three questions MBH raised against the §4.6 results, each answered by a further experiment.

**The detector's sample is too small to give a stable verdict.** `SRR26715166` was re-run under
twelve seeds. Re-seeding redraws the 1% detection sample, so the spread across seeds is the sampling
noise the 16x threshold has to survive.

| | |
|---|---|
| Sampled reads | 64–87 (mean 77) |
| Skew score | 12–22, mean 18.0, sd 2.9 |
| Below the 16x threshold | 2 of 12 seeds (17%) |
| Estimate when detected | 0.710x–0.795x |
| Estimate when missed | 0.215x |

The threshold sits 0.7 sd below the mean score, so the verdict is decided by which ~70 reads are
drawn rather than by the input. The paper seed, 4556, happens to draw one of the unlucky samples.
A floor on **sampled reads** is the fix. Note the detector already has a floor,
`MIN_DISTINCT_MINIMIZERS = 128`, but it guards distinct minimizers, and 69 reads of ~9.4 kb clear it
easily — the noise is in read count, not minimizer count, so this is a new guard rather than a
retuned one. Since the sd of a sampling statistic falls as 1/sqrt(n), moving from ~77 to ~1,000 reads
puts the threshold roughly 2.5 sd out. Treat 1,000 as a starting point for [#36][i36] to fit: the sd
comes from twelve seeds on one run, and a 99.9th percentile only roughly obeys the sqrt(n) scaling.
The floor binds only on small inputs — `SRR8618952` already samples 1,368 reads and `SRR12247681`
3,172 — so it costs almost nothing. It also will not fully fix this run: even when detected it
reaches 0.71x–0.80x rather than the 0.837x of `--normalize always`, because 7,473 reads cannot fill
the requested 15,000.

**Nothing separates the `SRR12247681` regression from the runs that needed fixing.** Its score of
104x sits mid-range among the engaging runs, which span 21x to 425x, so raising the threshold above
it would lose five genuine recoveries (`SRR24489322`, `DRR213976`, `SRR26715165`, `SRR10259778`,
`SRR26465563`) to save one run five points of error. The regression is also not the visible edge of a
trend: the eleven recovered runs average 0.962x, so normalization slightly *under*-corrects on
average, and 1.113x is a lone outlier the other way. One loose thread if more regressions appear —
normalization enriched this run's target set for longer reads, 40.0 MB of FASTA against 30.9 MB at
the same read count. What would settle the question is the engagement rate across the 3,370
benchmark runs, not a mechanism inferred from one control.

**The 10.2x cost is an implementation problem, not an inherent one.** `perf stat` on `SRR8618952`,
1.00 Gbp:

| | `never` | `auto` |
|---|---|---|
| task-clock | 16.9s | 65.2s |
| instructions | 121G | 345G |
| cache-misses | 197M (11.9%) | 1,042M (26.6%) |
| IPC | 2.09 | 1.44 |
| CPU utilisation | 242% | 113% |

`perf record` on a build with `CARGO_PROFILE_RELEASE_STRIP=false` puts **67.8% of the whole run in
`DepthSkewDetector::observe`**, plus 5.6% in the `values_u64()` fold inside it. Minimizer hashing is
2.3% and the entire minimap2 overlap computation — the work lrge exists to do — is about 13%. So the
detector costs roughly five times the overlap stage, and the cost is the CountMinSketch bookkeeping,
not the minimizers: `SKETCH_ROWS * SKETCH_WIDTH` is 4 x 2^20 u32, a 16 MB table, and `increment`
scatters four writes into it per minimizer. At roughly 200M minimizers that predicts ~800M extra
misses against the 844M measured.

Two separable defects. `count_records` is a sequential `while r.next()` loop, so the pass is
single-threaded (113% CPU) while the overlap stage uses eight cores. And `depth_sketch.increment`
runs for every minimizer of every read *before* the detector has decided anything, so on an unskewed
run every bit of that work is discarded. Building only the 1% detection sketch in the first pass and
the full depth profile in a second pass, taken only when skew fires, would cut the unskewed case by
about a hundredfold; the roughly 1% of runs that are skewed would pay one extra read of the input.
Parallelising the profiling pass is an independent second win.

Ticketed 2026-09-02 as [#46][i46] (defer the profile until the detector fires, then parallelise it),
[#47][i47] (cut the per-minimizer sketch cost) and [#48][i48] (the detection floor, at 500 reads).
All three are parented on [#34][i34] and block [#36][i36]: the full-benchmark re-run must not start
while `auto` costs 10.2x on unskewed input, and the floor changes which runs engage, so it has to
land before engagement is measured. [PR #45][pr45] itself is not held up by any of them.

Raw output under `/scratch/user/uqmhal11/lrge-issue34-benchmark/{seedvar,perf}/`; harness in
`seedvar.sh` and `perf.sh`.

#### Still open

- [ ] Run the detector across the 3,370 two-set runs to find how many engage. Nothing here bounds it.
- [ ] Decide which published figures move. Needs the item above.
- [ ] Re-derive the reported interval ([#38][i38]): every `auto` row here reports a wide IQR, and the
      quantiles were fitted on the old selection.
- [ ] Depth-profile the 10 PacBio sub-0.5x runs (carried from §6).
---

## 5. What the post needs to say about the paper

Draft position, to revisit once mechanism 1 lands:

- **Mechanism 2 changes nothing in the paper.** `-F` is off by default and the benchmark was run at
  defaults. What it changes is the *explanation* of 10 reported failures, and the usability of a
  documented flag. Say this explicitly — readers will assume a bug fix invalidates results.
- **Mechanism 1 will change reported numbers** on skewed runs. Scope TBD.
- The headline accuracy claim (75% of runs within 0.8–1.2x) is not at risk from either mechanism;
  both concern the ~1% tail.
- Users who reported this independently: [#8][i8], [#28][i28], [#29][i29]. Check their threads before
  publishing so the post answers what they actually asked.

---

## 6. Open questions

- [x] **Are the ONT sub-0.5x failures immune to mechanism 2?** No, and not uniformly. All 25
      sub-0.5x runs have now been tested: 18 respond to `-F`, 7 are completely inert (§2, §4.5).
- [x] **Are the 7 inert *N. gonorrhoeae* runs depth-skewed?** Yes, severely — 5–6 kb plasmids at
      226x–606x chromosomal depth. They are mechanism 1, and #34 should fix them (§4.5).
- [x] **Does internal-match fraction track depth skew?** No. It tracks read length against the size
      of the enriched element (§4.5).
- [ ] Are the PacBio sub-0.5x runs depth-skewed too? No profiling has been done on them, so the
      attribution in §3 rests on `-F` response alone. Same experiment, ten more runs.
- [ ] Does the read-length/element-size ratio predict `-F` response outside *N. gonorrhoeae*?
      `SRR26465560` (§4.3) is the obvious first check.
- [x] **Are the two mechanisms independent, or is mechanism 2 downstream of mechanism 1?** Both,
      depending on the run. Coupled on the eleven that normalization fixes, independent on the three
      it does not (§4.6).
- [ ] **How many of the 3,370 two-set runs does the detector engage on?** The 17 runs in §4.6 cannot
      bound this — the sample was chosen for being broken. Needs the detector run across the
      benchmark, and it gates any statement about which published figures move. Owned by
      [#36][i36], which also now asks whether anything separates the runs normalization improves
      from the already-correct ones it moves.
- [x] **Should the detection sample have a floor?** Yes. Twelve seeds on `SRR26715166` score
      12x–22x (sd 2.9) from 64–87 reads, and 2 of 12 fall under the 16x threshold. Floor set at 500
      reads, which puts the threshold ~1.75 sd out on that run (§4.6 follow-ups). Ticketed as
      [#48][i48]; constant to be refitted by [#36][i36].
- [x] **Fix the 10.2x cost on unskewed input.** Ticketed as [#46][i46] (defer the profile until the
      detector fires, then parallelise it) and [#47][i47] (cut the per-minimizer sketch cost). The
      detection floor is [#48][i48]. All three block [#36][i36], so the full-benchmark re-run cannot
      start until the cost is fixed.
- [ ] Should `--max-overhang-ratio` adapt to overlap density rather than being a constant? (§3.6)
- [ ] Should `-F` ever default on, for some detectable class of input? Current answer: no. (§3.5)
- [ ] Does the target-depth relationship in §3.4 hold on ONT data, and on skewed input? (§4.4)
- [ ] What is the right guidance for users whose read set is too small to fill `-T`? The 338- and
      238-target runs are pathological, and lrge warns but proceeds.

---

## 7. Reproduction

Working directories `/scratch/user/uqmhal11/lrge-issue29` (mechanism 2, §3) and
`/scratch/user/uqmhal11/lrge-issue34-benchmark` (§4.6; run logs, `/usr/bin/time -v` output and kept
temp dirs under `runs/<accession>/<mode>/`). **Both are scratch and not backed up.** The durable
artefacts have been copied here:

- [`rerun_estimates.tsv`](./rerun_estimates.tsv) — one row per run: read stats, target bases and
  depth, published lrge estimates, re-run estimates for all three variants, infinite-estimate counts.
- [`perread/`](./perread) — per-read estimates for every variant, gzipped. These are what the
  infinite-estimate counts are computed from.
- `nottested_31.sh` / `xoryzae_31.sh` — download, prep and run three variants for one accession.
  Prep mirrors `paper/workflow/scripts/download.sh`.
- `verify_31_clean.sh` — the default-behaviour-unchanged check (§3.2).
- `ont_F.sh` — `-F` on ONT accessions already prepped locally.
- [`depth/`](./depth) — 1 kb windowed depth for the ten *N. gonorrhoeae* runs of §4.5 and the two
  §4.6 controls (`SRR12247681`, `SRR8618952`), gzipped:
  `contig`, `window index`, `mean depth`. Produced by `depth_windows.sh`, which also prints the
  quantile summary the section quotes.
- `depth_windows.sh` — maps a prepped read set to its reference and emits the windowed profile.
- [`depth_normalization_estimates.tsv`](./depth_normalization_estimates.tsv) — the §4.6 benchmark of
  [PR #45][pr45] at commit `254e146`: one row per accession and mode, with estimate, truth ratio,
  reported interval, detector verdict and score, retained and total reads, wall time, peak RSS and
  exit status. `skew_score_source` records whether the score came from the run's own log or from a
  `-v` rerun, since the not-detected score is only emitted at debug level. Committed identically to
  this branch and to the PR #45 branch, so it lands once whichever merges first.
- `run_one.sh`, `batch.sh`, `build.sh` — the §4.6 harness: build the pinned binary once, then run one
  accession and mode per cell under `/usr/bin/time -v`.
- `identity_check.sh` / `identity_md5.sh` — the byte-identity check of §4.6, comparing the selected
  `target.fa` and `query.fa` between `auto` and `never`. lrge nests these in a randomly named
  subdirectory, so the comparison matches on basename.
- `diagnose.sh` — reruns an accession at `-v` to surface the debug-level skew score on runs where
  `auto` abstains, and under `--normalize always` to separate a detector miss from a normalization
  failure.
- `timing.sh` — alternates `never` and `auto` across three repeats on one input so page-cache warmth
  cannot favour either mode.
- `seedvar.sh` — re-runs one accession under twelve seeds to measure the sampling noise in the skew
  score, since re-seeding redraws the 1% detection sample.
- `perf.sh` / `perf_symbols.sh` — `perf stat` for `never` against `auto`, and a symbolised
  `perf record`. The release profile sets `strip = true`, so the symbolised build needs
  `CARGO_PROFILE_RELEASE_STRIP=false` as well as `CARGO_PROFILE_RELEASE_DEBUG=true`.
- `collect.py` / `summarise.py` — build `depth_normalization_estimates.tsv` from the run directories,
  and print the §4.6 comparisons (band crossings, mean error, `-F` benefit before and after
  normalization, control preservation, cost).
- `ont_31.sh` — the ONT equivalent of `nottested_31.sh` (`-x map-ont`, `-P ont`). Additionally runs
  the default variant under `--keep-temp` and reports the internal-match fraction of its overlap set
  as a `[PAF]` line, for prediction 1 of §4.4. The PAF is deleted afterwards.

Not copied: raw FASTQs (~17 GB) and lrge logs (~48 MB). The logs hold `Total target bases`, which is
the only field in the tables above not recoverable from the TSV.

Operational note: seven simultaneous ENA transfers tripped ascp auth throttling and two downloads
failed; re-running them serially succeeded against unchanged paths. Serialise if repeating.

---

## 8. Running log

- **2026-09-02** — the §4.6 follow-ups ticketed: [#46][i46] (defer depth profiling until the detector
  fires, then parallelise the pass), [#47][i47] (cut the per-minimizer sketch cost) and [#48][i48]
  (floor the detection sample at 500 reads). Parented on [#34][i34], all three blocking [#36][i36].
  [#36][i36] also gained a criterion asking whether anything separates the runs normalization
  improves from the already-correct ones it moves, which is the question one control could not answer.
- **2026-09-02** — three follow-ups to the §4.6 benchmark, prompted by MBH. The detection sample is
  too small to give a stable verdict: twelve seeds on `SRR26715166` score 12x–22x (sd 2.9) from
  64–87 reads, and 2 of 12 fall under the 16x threshold, so a read-count floor around 1,000 is
  warranted. The `SRR12247681` regression cannot be separated by threshold and is not a systematic
  bias, so no code change is proposed. The 10.2x cost is an implementation problem: `perf` puts 67.8%
  of the run in `DepthSkewDetector::observe` against ~13% for the whole overlap stage, and the pass
  is single-threaded. Recorded under §4.6.
- **2026-09-02** — [PR #45][pr45] (mechanism 1) benchmarked at commit `254e146` on all 15 ONT
  sub-0.5x runs plus two controls, 37 cells. Outliers in 0.8x–1.2x go from 0 of 15 to 11 of 15; mean
  absolute relative error 0.835 -> 0.222. All seven inert *N. gonorrhoeae* runs recover, which is the
  acceptance test §4.5 nominated. `--normalize never` reproduces `rerun_default` exactly on five runs,
  and the PacBio control's selected reads are byte-identical under `auto`. Three findings for the PR:
  the detector missed `SRR26715166` on a 69-read sample, where `--normalize always` recovers it;
  the healthy-looking ONT control `SRR12247681` is genuinely skewed (3 kb plasmid at 697x chromosomal
  depth), so the detector fires correctly but the estimate moves 0.937x -> 1.113x; and `auto` costs
  10.2x wall time on unskewed input. §4.6 written from this, §4.4 predictions 3 and 4 resolved.
- **2026-08-31** — #31 implemented, reviewed, merged (PR #40). Default behaviour verified unchanged
  on five accessions.
- **2026-08-31** — three *X. oryzae* runs re-run; mechanism 2 confirmed. Posted to #29.
- **2026-08-31** — remaining seven PacBio sub-0.5x runs re-run; all improve, four land within 0.8%
  of truth. Posted to [#29][c2].
- **2026-09-01** — MBH raised the coupling hypothesis now recorded as §4.4: depth enrichment should
  itself produce internal matches, making mechanism 2 partly a symptom of mechanism 1. Reframes what
  the post is about; cannot be settled until #34 lands.
- **2026-09-01** — depth profiles on all 10 *N. gonorrhoeae* sub-0.5x runs (7 inert + 3 responsive
  controls). All ten are plasmid-skewed; the inert ones are mechanism 1, not a third mechanism, and
  §4.5 was rewritten. `-F` response tracks read length against plasmid size, not skew magnitude.
- **2026-09-01** — remaining 13 ONT sub-0.5x runs completed, with internal-match fractions measured
  from the overlap set. All 25 sub-0.5x runs are now tested. Result: 18 respond to `-F`, 7 (all
  *N. gonorrhoeae*, all ~0.02x) are inert — provisionally read as a third failure mode, corrected the
  same day by the depth profiles above (§4.5). Internal-match
  fraction predicts response (§4.4, prediction 1). `-F` overshoots truth on 3 of 13 ONT runs, which
  strengthens §3.5.
- **2026-09-01** — `-F` tested post-fix on the two most extreme ONT failures. `SRR26465560`
  0.098x -> 0.643x (0.855x at ratio 0.05); `SRR26465526` 0.012x -> 0.013x. This refuted the
  platform split drafted in §2 and forced the rewrite of §2, §4.3 and §6. The run the mechanism-1
  spec is built around is substantially a mechanism-2 failure.

[doi]: https://doi.org/10.1093/bioinformatics/btaf593
[i8]: https://github.com/mbhall88/lrge/issues/8
[i28]: https://github.com/mbhall88/lrge/issues/28
[i29]: https://github.com/mbhall88/lrge/issues/29
[i31]: https://github.com/mbhall88/lrge/issues/31
[i32]: https://github.com/mbhall88/lrge/issues/32
[i33]: https://github.com/mbhall88/lrge/issues/33
[i34]: https://github.com/mbhall88/lrge/issues/34
[i35]: https://github.com/mbhall88/lrge/issues/35
[i36]: https://github.com/mbhall88/lrge/issues/36
[i37]: https://github.com/mbhall88/lrge/issues/37
[i38]: https://github.com/mbhall88/lrge/issues/38
[i39]: https://github.com/mbhall88/lrge/issues/39
[pr40]: https://github.com/mbhall88/lrge/pull/40
[pr45]: https://github.com/mbhall88/lrge/pull/45
[i46]: https://github.com/mbhall88/lrge/issues/46
[i47]: https://github.com/mbhall88/lrge/issues/47
[i48]: https://github.com/mbhall88/lrge/issues/48
[c2]: https://github.com/mbhall88/lrge/issues/29#issuecomment-5478078544
