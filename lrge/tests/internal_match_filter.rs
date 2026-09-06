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

fn estimate(input: &NamedTempFile, mode: &str, extra: &[&str]) -> (u64, String) {
    let mut arguments = ARGUMENTS.to_vec();
    arguments.extend(["--filter-contained", mode]);
    arguments.extend_from_slice(extra);
    run(input, &arguments)
}

#[test]
fn an_input_whose_overlaps_are_repeats_filters_them_without_being_told_to() {
    let input = reads_sharing_a_repeat(REPEAT_READS);
    let (never, never_log) = estimate(&input, "never", &[]);
    let (always, _) = estimate(&input, "always", &[]);
    let (auto, auto_log) = estimate(&input, "auto", &[]);

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
    assert!(!never_log.contains("Repeat-driven overlaps"));
}

#[test]
fn an_input_without_repeats_is_left_alone() {
    let input = reads_sharing_a_repeat(0);
    let (never, _) = estimate(&input, "never", &[]);
    let (auto, auto_log) = estimate(&input, "auto", &[]);

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
    let (never, _) = estimate(&input, "never", &["--use-min-ref"]);
    let (always, _) = estimate(&input, "always", &["--use-min-ref"]);
    let (auto, auto_log) = estimate(&input, "auto", &["--use-min-ref"]);

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
    let with = |mode: &str| {
        let mut arguments = arguments.to_vec();
        arguments.extend(["--filter-contained", mode]);
        run(&input, &arguments)
    };

    let (never, _) = with("never");
    let (always, _) = with("always");
    let (auto, auto_log) = with("auto");

    assert!(auto_log.contains("Repeat-driven overlaps detected"));
    assert_eq!(auto, always);
    assert!(never < auto);
}

/// Bare `-F` meant "filter everything" before it took a mode, and the estimate it gave has to be
/// the estimate it still gives.
#[test]
fn a_bare_filter_flag_still_filters_everything() {
    let input = reads_sharing_a_repeat(REPEAT_READS);
    let (always, _) = estimate(&input, "always", &[]);
    let mut arguments = ARGUMENTS.to_vec();
    arguments.push("-F");
    let (bare, _) = run(&input, &arguments);

    assert_eq!(bare, always);
}
