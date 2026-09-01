# Post-publication corrections: working notes

**Status: working notes, not the post.** This file accumulates everything needed to write a
technical update for lrge users once both mechanisms behind [#29][i29] are addressed. Mechanism 2
is fixed and measured. Mechanism 1 is specified but not implemented, so its impact numbers do not
exist yet. Sections marked **TODO** are the gaps.

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

The 25 sub-0.5x runs partition cleanly by platform, and the partition turns out to track the two
mechanisms:

- **10 PacBio** — all tested, all improve under the mechanism-2 fix. See §3.
- **15 ONT** — dominated by *N. gonorrhoeae* (9 runs) and *E. faecium* (2). Mechanism-1 territory.
  See §4.

That platform split is a genuine finding and belongs in the post, but note the caveat in §6: it is
currently an observation about which runs were tested, not a demonstration that ONT failures are
immune to mechanism 2.

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

## 4. Mechanism 1 — depth skew (specified, not implemented)

Issue [#29][i29] proper. Ticket chain: [#33][i33] (detect skew, observation only) → [#34][i34] (the
fix) → [#35][i35] (low-memory fallback) → [#36][i36] (fit constants) → [#37][i37] (release) →
[#38][i38] (re-derive quantiles). [#32][i32] and [#39][i39] are closed.

### 4.1 What is wrong

The median of per-read estimates collapses when a subpopulation of query reads comes from
high-depth sequence. `SRR26465560` (*E. faecium*, the run #29 was reported against) returns 315 kb
against a truth of 3.20 Mb, because six plasmids totalling 81 kb sit at up to 167x chromosomal
depth and supply the majority of the reads. The estimator returns the size of the over-represented
element rather than the genome.

### 4.2 The planned fix

Depth-aware read selection. Each read's depth is estimated as the **median** of its minimizer
counts — a median suppresses reads lying wholly inside a high-depth element while sparing reads that
merely span one, and spanning reads carry the most genome-size information. Reads are retained with
probability `min(1, C/depth)`. Target and query reads are drawn from a single normalized pool.
Normalization and subsampling happen together via weighted reservoir sampling, so the pass count is
unchanged. `--normalize auto|always|never`, defaulting to `auto`; an unskewed run takes the existing
path and produces byte-identical output.

### 4.3 Impact — **TODO**

This is the section that cannot be written yet, and it is the one that matters most for the post,
because **mechanism 1 changes default output and therefore changes published numbers**.

- [ ] Re-run the 15 ONT sub-0.5x runs after #34 lands; record before/after.
- [ ] Re-run a control set of healthy runs to confirm byte-identical output where skew is absent.
- [ ] Quantify how many of the 3,370 two-set runs change at all.
- [ ] Decide which published figures/tables would move, and by how much.
- [ ] Note the interaction with #38 (reported interval quantiles were fitted on the old behaviour).

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

- [ ] **Are the 15 ONT sub-0.5x failures immune to mechanism 2, or merely untested against it?**
      The platform split in §2 is currently an artefact of which runs were chosen. Running `-F` on
      the ONT failures would settle it. (Partially addressed — see §8.)
- [ ] Should `--max-overhang-ratio` adapt to overlap density rather than being a constant? (§3.6)
- [ ] Should `-F` ever default on, for some detectable class of input? Current answer: no. (§3.5)
- [ ] Does the target-depth relationship in §3.4 hold on ONT data?
- [ ] What is the right guidance for users whose read set is too small to fill `-T`? The 338- and
      238-target runs are pathological, and lrge warns but proceeds.

---

## 7. Reproduction

Working directory `/scratch/user/uqmhal11/lrge-issue29` (**volatile — scratch, not backed up**). The
durable artefacts have been copied here:

- [`rerun_estimates.tsv`](./rerun_estimates.tsv) — one row per run: read stats, target bases and
  depth, published lrge estimates, re-run estimates for all three variants, infinite-estimate counts.
- [`perread/`](./perread) — per-read estimates for every variant, gzipped. These are what the
  infinite-estimate counts are computed from.
- `nottested_31.sh` / `xoryzae_31.sh` — download, prep and run three variants for one accession.
  Prep mirrors `paper/workflow/scripts/download.sh`.
- `verify_31_clean.sh` — the default-behaviour-unchanged check (§3.2).
- `ont_F.sh` — `-F` on ONT accessions already prepped locally.

Not copied: raw FASTQs (~17 GB) and lrge logs (~48 MB). The logs hold `Total target bases`, which is
the only field in the tables above not recoverable from the TSV.

Operational note: seven simultaneous ENA transfers tripped ascp auth throttling and two downloads
failed; re-running them serially succeeded against unchanged paths. Serialise if repeating.

---

## 8. Running log

- **2026-08-31** — #31 implemented, reviewed, merged (PR #40). Default behaviour verified unchanged
  on five accessions.
- **2026-08-31** — three *X. oryzae* runs re-run; mechanism 2 confirmed. Posted to #29.
- **2026-08-31** — remaining seven PacBio sub-0.5x runs re-run; all improve, four land within 0.8%
  of truth. Posted to [#29][c2].
- **2026-09-01** — `-F` tested on the two most extreme ONT failures (SRR26465560, SRR26465526) to
  probe the platform split in §2. **TODO: record result here.**

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
[c2]: https://github.com/mbhall88/lrge/issues/29#issuecomment-5478078544
