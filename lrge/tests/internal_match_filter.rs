mod common;

use common::{circular_read, pseudo_random_dna, run};
use std::io::Write;
use tempfile::NamedTempFile;

const CHROMOSOME_SIZE: usize = 60_000;
const FLANK_SIZE: usize = 1_500;
/// Long enough that two reads sharing it align over a stretch several times their overhangs.
const REPEAT_SIZE: usize = 3_000;
const READ_LENGTH: usize = FLANK_SIZE + REPEAT_SIZE + FLANK_SIZE;
const CHROMOSOME_READS: usize = 200;
const REPEAT_READS: usize = 400;

/// What every test here runs under, less the filter mode it appends.
///
/// Depth normalization would also fire on these inputs, because a repeat carried by two thirds of
/// the reads is exactly what it looks for. Turning it off is what leaves the filter as the only
/// thing that can move the estimate.
const ARGUMENTS: [&str; 9] = [
    "-T",
    "200",
    "-Q",
    "100",
    "--seed",
    "42",
    "-vv",
    "--normalize",
    "never",
];

fn write_read(input: &mut NamedTempFile, name: &str, sequence: &[u8]) {
    writeln!(input, ">{name}\n{}", String::from_utf8_lossy(sequence)).unwrap();
}

/// Reads tiling a chromosome, plus `repeat_reads` reads that share a repeat and nothing else.
///
/// A repeat-bearing read is unique sequence, then the repeat, then more unique sequence. Two of
/// them align over the repeat in the middle with a flank hanging off either end, which is the
/// internal match `--filter-contained` exists to remove. Two chromosome reads overlap end to end
/// instead, so they are what the filter must leave alone.
fn reads_sharing_a_repeat(repeat_reads: usize) -> NamedTempFile {
    let chromosome = pseudo_random_dna(CHROMOSOME_SIZE, 1);
    let repeat = pseudo_random_dna(REPEAT_SIZE, 2);
    // Every flank is cut from its own stretch of this, so no two repeat-bearing reads share
    // anything but the repeat.
    let flanks = pseudo_random_dna(repeat_reads.max(1) * 2 * FLANK_SIZE, 3);
    let mut input = NamedTempFile::new().unwrap();

    for index in 0..CHROMOSOME_READS {
        let start = index * 137 % CHROMOSOME_SIZE;
        write_read(
            &mut input,
            &format!("chromosome{index}"),
            &circular_read(&chromosome, start, READ_LENGTH),
        );
    }

    for index in 0..repeat_reads {
        let left = index * 2 * FLANK_SIZE;
        let mut read = flanks[left..left + FLANK_SIZE].to_vec();
        read.extend_from_slice(&repeat);
        read.extend_from_slice(&flanks[left + FLANK_SIZE..left + 2 * FLANK_SIZE]);
        write_read(&mut input, &format!("repeat{index}"), &read);
    }

    input
}

/// Run with the filter off, or at a given share.
///
/// The share is written with an equals sign because that is the only form `-F` takes. A share of
/// zero excludes whatever internal matches the run finds, which is what asking for the filter
/// unconditionally amounts to and what these fixtures are built to trigger.
fn estimate(input: &NamedTempFile, share: Option<&str>, extra: &[&str]) -> (u64, String) {
    let mut arguments = ARGUMENTS.to_vec();
    if let Some(share) = share {
        arguments.push(share);
    }
    arguments.extend_from_slice(extra);
    run(input, &arguments)
}

#[test]
fn an_input_whose_overlaps_are_repeats_filters_them_without_being_told_to() {
    let input = reads_sharing_a_repeat(REPEAT_READS);
    let (never, never_log) = estimate(&input, None, &[]);
    let (always, _) = estimate(&input, Some("-F=0"), &[]);
    let (auto, auto_log) = estimate(&input, Some("-F=0.8"), &[]);

    assert!(
        auto_log.contains("Repeat-driven overlaps detected"),
        "an input built entirely out of internal matches was not detected:\n{auto_log}"
    );
    assert!(auto_log
        .lines()
        .any(|line| line.contains("WARN") && line.contains("Repeat-driven overlaps detected")));
    assert_eq!(
        auto, always,
        "auto detected the repeats but did not estimate as if it had"
    );
    assert!(
        never < auto,
        "keeping the internal matches ({never}) should read lower than dropping them ({auto})"
    );
    // a run that never asked still says what it measured, so a threshold can be chosen from it
    assert!(!never_log.contains("Repeat-driven overlaps"));
    assert!(never_log.contains("Overlap composition: internal matches account for"));
}

#[test]
fn an_input_without_repeats_is_left_alone() {
    let input = reads_sharing_a_repeat(0);
    let (never, _) = estimate(&input, None, &[]);
    let (auto, auto_log) = estimate(&input, Some("-F=0.8"), &[]);

    assert!(
        auto_log.contains("Repeat-driven overlaps not detected"),
        "an input with no repeats was called repeat-driven:\n{auto_log}"
    );
    assert_eq!(
        auto, never,
        "auto found no repeats but still moved the estimate"
    );
}

/// `--use-min-ref` counts overlaps on the other side of the mapping, in its own loop, so it needs
/// its own check that the verdict reaches the estimate.
#[test]
fn the_inverse_mapping_path_reaches_the_same_verdict() {
    let input = reads_sharing_a_repeat(REPEAT_READS);
    let (never, _) = estimate(&input, None, &["--use-min-ref"]);
    let (always, _) = estimate(&input, Some("-F=0"), &["--use-min-ref"]);
    let (auto, auto_log) = estimate(&input, Some("-F=0.8"), &["--use-min-ref"]);

    assert!(auto_log.contains("Repeat-driven overlaps detected"));
    assert_eq!(auto, always);
    assert!(never < auto);
}

/// All-vs-all counts overlaps as read pairs rather than per query read, so it too has its own loop.
#[test]
fn all_vs_all_reaches_the_same_verdict() {
    let input = reads_sharing_a_repeat(REPEAT_READS);
    let arguments = [
        "--num",
        "300",
        "--seed",
        "42",
        "-vv",
        "--normalize",
        "never",
    ];
    let with = |share: Option<&str>| {
        let mut arguments = arguments.to_vec();
        if let Some(share) = share {
            arguments.push(share);
        }
        run(&input, &arguments)
    };

    let (never, _) = with(None);
    let (always, _) = with(Some("-F=0"));
    let (auto, auto_log) = with(Some("-F=0.8"));

    assert!(auto_log.contains("Repeat-driven overlaps detected"));
    assert_eq!(auto, always);
    assert!(never < auto);
}

/// Bare `-F` is the fitted share, so on an input built out of internal matches it filters, and on
/// one built without them it does not.
#[test]
fn a_bare_filter_flag_uses_the_fitted_share() {
    let repeats = reads_sharing_a_repeat(REPEAT_READS);
    let (at_default, _) = estimate(&repeats, Some("-F=0.8"), &[]);
    let mut arguments = ARGUMENTS.to_vec();
    arguments.push("-F");
    let (bare, _) = run(&repeats, &arguments);

    assert_eq!(bare, at_default);

    let plain = reads_sharing_a_repeat(0);
    let (unfiltered, _) = estimate(&plain, None, &[]);
    let mut arguments = ARGUMENTS.to_vec();
    arguments.push("-F");
    let (bare_on_plain, _) = run(&plain, &arguments);

    assert_eq!(bare_on_plain, unfiltered);
}

/// The share the caller gives is the one that decides, so the same input filters at a low threshold
/// and does not at one above its share.
#[test]
fn the_share_the_caller_gives_is_the_one_that_decides() {
    let input = reads_sharing_a_repeat(REPEAT_READS);
    let (unfiltered, _) = estimate(&input, None, &[]);
    let (low, low_log) = estimate(&input, Some("-F=0.1"), &[]);
    let (high, high_log) = estimate(&input, Some("-F=0.99"), &[]);

    assert!(low_log.contains("Repeat-driven overlaps detected"));
    assert!(low_log.contains("against a threshold of 10.0%"));
    assert!(high_log.contains("Repeat-driven overlaps not detected"));
    assert!(high_log.contains("against a threshold of 99.0%"));
    assert!(unfiltered < low);
    assert_eq!(
        high, unfiltered,
        "a share nothing can clear must change nothing"
    );
}
