//! Deciding whether a run's overlaps come from repeats rather than from shared loci.
//!
//! An internal match is two reads that align over a stretch in the middle of both, leaving long
//! unaligned tails hanging off either end. Reads from the same locus do not look like that; reads
//! that happen to carry the same repeat do. `PafRecord::is_internal` is the test, and
//! `--filter-contained` is the switch that throws those mappings away.
//!
//! Throwing them away helps enormously on some inputs and hurts on most. A per-read estimate
//! divides the target count by the read's overlap count, so removing overlaps raises every
//! estimate: over 27 accessions whose genome size is known, filtering moved the estimate up on
//! every single one, by 8% on the mildest and by more than sevenfold on the worst. That is why the
//! filter is off by default, and why turning it on unconditionally trades one failure mode for
//! another.
//!
//! What this module adds is the third option. The same pass that collects overlaps can count what
//! the filter would have discarded without discarding it, so a run can measure how much of its
//! overlap evidence rests on repeats and decide for itself. See [`InternalMatchTally`].

use std::fmt;
use std::sync::atomic::{AtomicU64, Ordering};

use log::{debug, warn};

use crate::InternalFilter;

// How much of the overlap evidence has to be internal matches before a run filters them out.
//
// This was fitted the way the depth-skew constants were, on runs whose genome size is known, but on
// far fewer of them: the 27 accessions held locally rather than the paper's whole benchmark. Those
// 27 were chosen to be rich in the failures depth normalization leaves behind, so they say more
// about where the filter helps than about how often it is asked for.
//
// What those 27 say is that the runs already within 10% of their true size are all quiet: the
// highest of the thirteen sits at 0.607, so any threshold above that leaves every one of them
// untouched. The runs the filter rescues are not so tidy. Six of the seven sit at 0.778 or above,
// well clear of the quiet population, but the seventh sits at 0.424, in among them, and no
// threshold can take it without taking correct runs too.
//
// So the useful range is everything between those two populations, and across the whole of
// (0.607, 0.778) the outcome is the same: 19 of the 27 land within 10% of the truth, against 13 for
// never filtering and 9 for always filtering, and no run that was already correct is disturbed.
// This sits in the middle of that range, so neither edge is close.
//
// Four runs still fail under it, and they fail because they need the filter without looking like
// it. `SRR10259778` is the clearest: filtering multiplies its estimate by 3.79 and brings it into
// range, on an internal-match share of 0.26. Whatever separates it from the runs that sit at the
// same share and are already correct is not in this statistic.
//
// Those figures come from an instrumented pass, in
// `paper/corrections/issue36_internal_match_probe.tsv`, and were then reproduced by running the
// same 27 accessions through this code in all three modes:
// `paper/corrections/issue36_dynamic_filter.tsv`. The argument is in
// `paper/corrections/README_issue36.md`.
pub(crate) const INTERNAL_MATCH_THRESHOLD: f64 = 0.7;

/// Counts what internal-match filtering would discard, on a run that has not decided to discard it.
///
/// Both counts are of *unique overlaps*, because that is what a per-read estimate divides by. A
/// query read that maps to the same target five times is one overlap either way, and a target that
/// a read reaches only through internal matches is an overlap the filter would take away.
///
/// This is shared across the mapping threads and holds two counters, so the cost of carrying it is
/// two relaxed adds per call to [`observe`][Self::observe].
#[derive(Debug, Default)]
pub(crate) struct InternalMatchTally {
    unfiltered: AtomicU64,
    kept: AtomicU64,
}

impl InternalMatchTally {
    /// Fold in a count of overlaps and however many of them would survive the filter.
    ///
    /// Callers that can attribute overlaps to a single read call this once per read; the
    /// all-vs-all strategy counts read pairs, which belong to two reads at once, so it calls this
    /// once for the whole pass instead. Either way the totals are what the verdict is read off.
    pub(crate) fn observe(&self, unfiltered: usize, kept: usize) {
        debug_assert!(kept <= unfiltered);
        self.unfiltered
            .fetch_add(unfiltered as u64, Ordering::Relaxed);
        self.kept.fetch_add(kept as u64, Ordering::Relaxed);
    }

    /// Read the verdict off the counts. Call this once the mapping pass is done.
    pub(crate) fn report(&self) -> InternalMatchReport {
        InternalMatchReport {
            unfiltered: self.unfiltered.load(Ordering::Relaxed),
            kept: self.kept.load(Ordering::Relaxed),
        }
    }
}

/// What an [`InternalMatchTally`] came to, and whether it asks for the filter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct InternalMatchReport {
    unfiltered: u64,
    kept: u64,
}

impl InternalMatchReport {
    /// The share of the run's overlaps that internal matches account for.
    ///
    /// `None` when there were no overlaps at all, which is a run with nothing to decide about.
    pub(crate) fn drop_fraction(&self) -> Option<f64> {
        if self.unfiltered == 0 {
            return None;
        }
        Some(1.0 - (self.kept as f64 / self.unfiltered as f64))
    }

    /// Whether enough of the overlaps are internal matches to filter them out.
    pub(crate) fn repeat_driven(&self) -> bool {
        self.drop_fraction()
            .is_some_and(|fraction| fraction > INTERNAL_MATCH_THRESHOLD)
    }
}

impl fmt::Display for InternalMatchReport {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.drop_fraction() {
            Some(fraction) => {
                let verdict = if self.repeat_driven() {
                    "Repeat-driven overlaps detected"
                } else {
                    "Repeat-driven overlaps not detected"
                };
                write!(
                    formatter,
                    "{verdict} (internal matches account for {:.1}% of {} overlaps)",
                    fraction * 100.0,
                    self.unfiltered
                )
            }
            None => write!(
                formatter,
                "Repeat-driven overlaps not assessed (no overlaps)"
            ),
        }
    }
}

/// Whether a run should exclude its internal matches, given what it was asked for and what the
/// mapping pass measured.
///
/// Only [`InternalFilter::Auto`] consults the tally, and only it says anything about the tally in
/// the log. A run that was told what to do has nothing to report.
pub(crate) fn filter_internal_matches(mode: InternalFilter, tally: &InternalMatchTally) -> bool {
    match mode {
        InternalFilter::Never => false,
        InternalFilter::Always => true,
        InternalFilter::Auto => {
            let report = tally.report();
            // Filtering moves the estimate on every input, so a run that engages it is a run whose
            // answer would have been very different a moment ago. That is worth a warning.
            if report.repeat_driven() {
                warn!("{report}");
            } else {
                debug!("{report}");
            }
            report.repeat_driven()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn report(unfiltered: u64, kept: u64) -> InternalMatchReport {
        let tally = InternalMatchTally::default();
        tally.observe(unfiltered as usize, kept as usize);
        tally.report()
    }

    #[test]
    fn a_run_with_no_overlaps_has_no_fraction_and_no_verdict() {
        let report = report(0, 0);

        assert_eq!(report.drop_fraction(), None);
        assert!(!report.repeat_driven());
        assert!(report.to_string().contains("not assessed"));
    }

    #[test]
    fn the_fraction_is_the_share_the_filter_would_remove() {
        assert_eq!(report(100, 25).drop_fraction(), Some(0.75));
        assert_eq!(report(100, 100).drop_fraction(), Some(0.0));
    }

    #[test]
    fn the_threshold_is_exclusive_so_a_run_exactly_on_it_is_left_alone() {
        let on_it = report(
            1000,
            (1000.0 * (1.0 - INTERNAL_MATCH_THRESHOLD)).round() as u64,
        );

        assert_eq!(on_it.drop_fraction(), Some(INTERNAL_MATCH_THRESHOLD));
        assert!(!on_it.repeat_driven());
    }

    /// The gap the constant was fitted into: the noisiest run that is already correct, and the
    /// quietest run the filter has to rescue. Neither may end up on the wrong side of the
    /// threshold. See [`INTERNAL_MATCH_THRESHOLD`].
    #[test]
    fn the_fitted_populations_fall_either_side_of_the_threshold() {
        let at = |fraction: f64| {
            let total = 1_000_000;
            report(total, (total as f64 * (1.0 - fraction)).round() as u64)
        };

        assert!(
            !at(0.607_020).repeat_driven(),
            "SRR26465560 is already correct"
        );
        assert!(
            at(0.777_531).repeat_driven(),
            "SRR30357566 needs the filter"
        );
    }

    #[test]
    fn counts_accumulate_across_reads() {
        let tally = InternalMatchTally::default();
        tally.observe(10, 1);
        tally.observe(30, 9);

        assert_eq!(tally.report().drop_fraction(), Some(0.75));
    }

    #[test]
    fn a_told_run_ignores_the_tally() {
        let repeat_driven = InternalMatchTally::default();
        repeat_driven.observe(100, 1);

        assert!(!filter_internal_matches(
            InternalFilter::Never,
            &repeat_driven
        ));
        assert!(filter_internal_matches(
            InternalFilter::Always,
            &InternalMatchTally::default()
        ));
    }

    #[test]
    fn an_auto_run_follows_the_tally() {
        let repeat_driven = InternalMatchTally::default();
        repeat_driven.observe(100, 1);
        let ordinary = InternalMatchTally::default();
        ordinary.observe(100, 99);

        assert!(filter_internal_matches(
            InternalFilter::Auto,
            &repeat_driven
        ));
        assert!(!filter_internal_matches(InternalFilter::Auto, &ordinary));
    }

    /// A run with no overlaps at all has nothing to correct, and the correction it would make is
    /// the one that cannot be undone: filtering can only take overlaps away.
    #[test]
    fn an_auto_run_that_saw_nothing_does_not_filter() {
        assert!(!filter_internal_matches(
            InternalFilter::Auto,
            &InternalMatchTally::default()
        ));
    }

    #[test]
    fn the_verdict_names_itself_either_way() {
        assert!(report(100, 10).to_string().contains("detected (internal"));
        assert!(report(100, 90).to_string().contains("not detected"));
    }
}
