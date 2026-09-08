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
//! What this module adds is a middle setting. The same pass that collects overlaps can count what
//! the filter would have discarded without discarding it, so a run can measure how much of its
//! overlap evidence rests on repeats and act only when that share is high. The share it has to
//! clear is the caller's to choose, because the benchmark says the best value depends on what the
//! caller is willing to trade: see [`DEFAULT_INTERNAL_MATCH_THRESHOLD`][crate::DEFAULT_INTERNAL_MATCH_THRESHOLD].
//! See [`InternalMatchTally`].

use std::fmt;
use std::sync::atomic::{AtomicU64, Ordering};

use log::{debug, warn};

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

    /// Read the verdict off the counts, against the share the run had to clear. Call this once the
    /// mapping pass is done.
    pub(crate) fn report(&self, threshold: Option<f64>) -> InternalMatchReport {
        InternalMatchReport {
            unfiltered: self.unfiltered.load(Ordering::Relaxed),
            kept: self.kept.load(Ordering::Relaxed),
            threshold,
        }
    }
}

/// What an [`InternalMatchTally`] came to, and the share a run had to clear, if it asked to.
///
/// Every run measures, because the count rides along with the one the estimate needs and a share a
/// run cannot see is a threshold nobody can choose. Only a run that asked for the filter has a
/// threshold to be measured against, and only that run gets a verdict.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct InternalMatchReport {
    unfiltered: u64,
    kept: u64,
    threshold: Option<f64>,
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
    ///
    /// False for a run that did not ask, however repeat-driven it turns out to be.
    pub(crate) fn repeat_driven(&self) -> bool {
        match (self.drop_fraction(), self.threshold) {
            (Some(fraction), Some(threshold)) => fraction > threshold,
            _ => false,
        }
    }
}

impl fmt::Display for InternalMatchReport {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Some(fraction) = self.drop_fraction() else {
            return write!(
                formatter,
                "Repeat-driven overlaps not assessed (no overlaps)"
            );
        };
        let share = format!(
            "internal matches account for {:.1}% of {} overlaps",
            fraction * 100.0,
            self.unfiltered
        );
        match self.threshold {
            // A run that did not ask for the filter gets the measurement and no verdict, because
            // the verdict would be against a threshold it never chose.
            None => write!(formatter, "Overlap composition: {share}"),
            Some(threshold) => {
                let verdict = if self.repeat_driven() {
                    "Repeat-driven overlaps detected"
                } else {
                    "Repeat-driven overlaps not detected"
                };
                write!(
                    formatter,
                    "{verdict} ({share}, against a threshold of {:.1}%)",
                    threshold * 100.0
                )
            }
        }
    }
}

/// Whether a run should exclude its internal matches, given the share it was asked to clear and
/// what the mapping pass measured.
///
/// `None` is a run that never asked for the filter. It still reports what it measured, because a
/// share nobody can see is a threshold nobody can choose, but it never acts on it. A threshold of
/// zero filters whatever internal matches the run found, which is what asking for the filter
/// unconditionally amounts to.
pub(crate) fn filter_internal_matches(threshold: Option<f64>, tally: &InternalMatchTally) -> bool {
    let report = tally.report(threshold);
    // Filtering moves the estimate on every input, so a run that engages it is a run whose answer
    // would have been very different a moment ago. That is worth a warning.
    if report.repeat_driven() {
        warn!("{report}");
    } else {
        debug!("{report}");
    }
    report.repeat_driven()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DEFAULT_INTERNAL_MATCH_THRESHOLD as DEFAULT;

    fn report(unfiltered: u64, kept: u64) -> InternalMatchReport {
        at(unfiltered, kept, DEFAULT)
    }

    fn at(unfiltered: u64, kept: u64, threshold: f64) -> InternalMatchReport {
        let tally = InternalMatchTally::default();
        tally.observe(unfiltered as usize, kept as usize);
        tally.report(Some(threshold))
    }

    /// A share, expressed as the counts a run would have produced to reach it.
    fn share(fraction: f64) -> InternalMatchReport {
        let total = 1_000_000;
        report(total, (total as f64 * (1.0 - fraction)).round() as u64)
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
        let on_it = share(DEFAULT);

        assert_eq!(on_it.drop_fraction(), Some(DEFAULT));
        assert!(!on_it.repeat_driven());
    }

    /// What the benchmark asks of the default: leave alone every run the estimator already puts
    /// within 10% of the truth. The noisiest of those over 3,370 accessions is `SRR13170267`, at a
    /// share of 0.796, and filtering it would take it from 1.095x to 4.169x. See
    /// [`DEFAULT_INTERNAL_MATCH_THRESHOLD`][crate::DEFAULT_INTERNAL_MATCH_THRESHOLD].
    #[test]
    fn the_noisiest_already_correct_run_stays_below_the_default() {
        assert!(
            !share(0.796).repeat_driven(),
            "SRR13170267 is already correct"
        );
        assert!(share(0.845).repeat_driven(), "SRR30357565 needs the filter");
    }

    #[test]
    fn counts_accumulate_across_reads() {
        let tally = InternalMatchTally::default();
        tally.observe(10, 1);
        tally.observe(30, 9);

        assert_eq!(tally.report(Some(DEFAULT)).drop_fraction(), Some(0.75));
    }

    #[test]
    fn a_run_that_never_asked_does_not_filter() {
        let repeat_driven = InternalMatchTally::default();
        repeat_driven.observe(100, 1);

        assert!(!filter_internal_matches(None, &repeat_driven));
    }

    /// A run that never asked still measures, and says what it found without passing judgement on
    /// it. That is what lets someone choose a threshold from one ordinary run.
    #[test]
    fn a_run_that_never_asked_reports_the_share_without_a_verdict() {
        let tally = InternalMatchTally::default();
        tally.observe(100, 10);
        let line = tally.report(None).to_string();

        assert!(
            line.contains("account for 90.0% of 100 overlaps"),
            "got: {line}"
        );
        assert!(
            !line.contains("detected"),
            "a run that did not ask gets no verdict: {line}"
        );
        assert!(!line.contains("threshold"), "got: {line}");
    }

    /// Asking for the filter unconditionally is a threshold of zero, which is what the flag meant
    /// before it took a share.
    #[test]
    fn a_threshold_of_zero_filters_whatever_the_run_found() {
        let barely = InternalMatchTally::default();
        barely.observe(1_000_000, 999_999);

        assert!(filter_internal_matches(Some(0.0), &barely));
    }

    /// A run with no internal matches at all has none to remove, so even a threshold of zero
    /// leaves its estimate alone. The two are the same answer by different routes.
    #[test]
    fn a_run_with_nothing_to_filter_is_not_filtered_at_any_threshold() {
        let clean = InternalMatchTally::default();
        clean.observe(1_000, 1_000);

        assert!(!filter_internal_matches(Some(0.0), &clean));
    }

    #[test]
    fn the_threshold_the_caller_gives_is_the_one_that_decides() {
        let two_thirds = InternalMatchTally::default();
        two_thirds.observe(300, 100);

        assert!(filter_internal_matches(Some(0.5), &two_thirds));
        assert!(!filter_internal_matches(Some(0.9), &two_thirds));
    }

    #[test]
    fn an_auto_run_that_saw_nothing_does_not_filter() {
        assert!(!filter_internal_matches(
            Some(DEFAULT),
            &InternalMatchTally::default()
        ));
    }

    #[test]
    fn the_verdict_names_itself_either_way_and_says_what_it_measured_against() {
        assert!(report(100, 10).to_string().contains("detected (internal"));
        assert!(report(100, 90).to_string().contains("not detected"));
        assert!(at(100, 10, 0.5).to_string().contains("threshold of 50.0%"));
    }
}
