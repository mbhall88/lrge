mod common;

use common::{next_draw, pseudo_random_dna, run, write_reads};
use tempfile::NamedTempFile;

const CHROMOSOME_SIZE: usize = 20_000;
const ELEMENT_SIZE: usize = 2_000;
const READ_LENGTH: usize = 800;
const CHROMOSOME_READS: usize = 400;
const ELEMENT_READS: usize = 4_000;

/// A chromosome long enough that [`LOW_DEPTH_CHROMOSOME_READS`] of it is eight-fold coverage.
const LOW_DEPTH_CHROMOSOME_SIZE: usize = 120_000;
const LOW_DEPTH_READ_LENGTH: usize = 2_000;
const LOW_DEPTH_CHROMOSOME_READS: usize = 500;
/// Enough copies of the element that it is four fifths of the reads, as a plasmid can be.
const LOW_DEPTH_ELEMENT_READS: usize = 2_000;
/// One base in this many is replaced at random, which is about what a nanopore read carries.
///
/// The errors are what put these inputs at a profile median depth of one. Nearly every minimizer
/// an error creates is seen in the one read that has it, and at this read length and error rate
/// those singletons outnumber the minimizers the genome repeats, so the median count over the
/// detection sample is one however deep the sequencing is. Error-free reads cannot reach that
/// regime at a coverage the estimator can still work at: their minimizer counts are the coverage.
const LOW_DEPTH_ERROR_DENOMINATOR: u64 = 10;
/// Where the error generator starts, so an input is the same one every run.
const LOW_DEPTH_ERROR_SEED: u64 = 99;

fn circular_read(sequence: &[u8], start: usize) -> Vec<u8> {
    common::circular_read(sequence, start, READ_LENGTH)
}

/// A read of `source`, wrapping at its end, with about one base in
/// [`LOW_DEPTH_ERROR_DENOMINATOR`] replaced at random.
fn error_prone_read(source: &[u8], start: usize, state: &mut u64) -> Vec<u8> {
    let mut read: Vec<u8> = (0..LOW_DEPTH_READ_LENGTH)
        .map(|offset| source[(start + offset) % source.len()])
        .collect();
    for base in &mut read {
        if next_draw(state).is_multiple_of(LOW_DEPTH_ERROR_DENOMINATOR) {
            *base = b"ACGT"[(next_draw(state) & 3) as usize];
        }
    }
    read
}

/// What both low-coverage tests run under, less the normalization mode each appends.
///
/// `-vv` is what puts the profile's median depth in the log, which is the number these tests are
/// really about.
const LOW_DEPTH_ARGUMENTS: [&str; 8] = [
    "-T",
    "200",
    "-Q",
    "100",
    "--seed",
    "42",
    "-vv",
    "--normalize",
];

/// Error-prone reads of a chromosome at eight-fold coverage, plus `element_reads` reads of a
/// high-copy element.
///
/// Ask for no element reads to get the even-depth version of the same input, which is what says
/// whether a median depth of one is enough to call an ordinary input skewed by mistake.
fn low_depth_reads(element_reads: usize) -> NamedTempFile {
    let chromosome = pseudo_random_dna(LOW_DEPTH_CHROMOSOME_SIZE, 1);
    let element = pseudo_random_dna(ELEMENT_SIZE, 2);
    let mut state = LOW_DEPTH_ERROR_SEED;
    let mut input = NamedTempFile::new().unwrap();

    write_reads(
        &mut input,
        "chromosome",
        &chromosome,
        LOW_DEPTH_CHROMOSOME_READS,
        |source, start| error_prone_read(source, start, &mut state),
    );
    write_reads(
        &mut input,
        "element",
        &element,
        element_reads,
        |source, start| error_prone_read(source, start, &mut state),
    );

    input
}

fn skewed_reads() -> NamedTempFile {
    let chromosome = pseudo_random_dna(CHROMOSOME_SIZE, 1);
    let element = pseudo_random_dna(ELEMENT_SIZE, 2);
    let mut input = NamedTempFile::new().unwrap();

    write_reads(
        &mut input,
        "chromosome",
        &chromosome,
        CHROMOSOME_READS,
        circular_read,
    );
    write_reads(
        &mut input,
        "element",
        &element,
        ELEMENT_READS,
        circular_read,
    );

    input
}

#[test]
fn two_set_normalization_corrects_a_high_copy_element() {
    let input = skewed_reads();
    let (legacy, legacy_log) = run(
        &input,
        &[
            "-T",
            "200",
            "-Q",
            "100",
            "--seed",
            "42",
            "--normalize",
            "never",
        ],
    );
    let (normalized, normalized_log) = run(
        &input,
        &[
            "-T",
            "200",
            "-Q",
            "100",
            "--seed",
            "42",
            "--normalize",
            "auto",
        ],
    );

    assert!(legacy < 10_000, "legacy estimate was {legacy}");
    assert!(
        (15_000..=35_000).contains(&normalized),
        "normalized estimate was {normalized}"
    );
    assert!(!legacy_log.contains("depth normalization retained"));
    assert!(normalized_log.contains("Depth skew detected"));
    assert!(normalized_log.contains("depth normalization retained"));
    assert!(normalized_log
        .lines()
        .any(|line| line.contains("WARN") && line.contains("Depth skew detected")));
}

#[test]
fn all_vs_all_normalization_corrects_a_high_copy_element() {
    let input = skewed_reads();
    let (legacy, _) = run(
        &input,
        &["--num", "300", "--seed", "42", "--normalize", "never"],
    );
    let (normalized, normalized_log) = run(
        &input,
        &["--num", "300", "--seed", "42", "--normalize", "auto"],
    );

    assert!(legacy < 10_000, "legacy estimate was {legacy}");
    assert!(
        (15_000..=35_000).contains(&normalized),
        "normalized estimate was {normalized}"
    );
    assert!(normalized_log.contains("depth normalization retained"));
    assert!(normalized_log
        .lines()
        .any(|line| line.contains("WARN") && line.contains("Depth skew detected")));
}

#[test]
fn a_small_normalized_pool_warns_and_still_estimates() {
    let input = skewed_reads();
    let (estimate, log) = run(
        &input,
        &[
            "-T",
            "800",
            "-Q",
            "800",
            "--seed",
            "42",
            "--normalize",
            "always",
        ],
    );

    assert!(estimate > 0);
    assert!(log.contains("Depth normalization forced"));
    assert!(log
        .lines()
        .any(|line| line.contains("INFO") && line.contains("Depth normalization forced")));
    assert!(log.contains("Normalized read pool is smaller than requested"));
}

#[test]
fn the_low_memory_path_reproduces_the_buffered_estimate() {
    let input = skewed_reads();
    let arguments = [
        "-T",
        "200",
        "-Q",
        "100",
        "--seed",
        "42",
        "--normalize",
        "auto",
    ];
    let (buffered, buffered_log) = run(&input, &arguments);
    let mut low_memory_arguments = arguments.to_vec();
    low_memory_arguments.extend(["--max-read-buffer", "1"]);
    let (low_memory, low_memory_log) = run(&input, &low_memory_arguments);

    assert_eq!(
        buffered, low_memory,
        "buffered estimate {buffered} differs from low-memory estimate {low_memory}"
    );
    assert!(!buffered_log.contains("Selected reads by position"));
    assert!(low_memory_log.contains("Selected reads by position"));
    assert!(low_memory_log.contains("depth normalization retained"));
}

#[test]
fn a_generous_budget_keeps_the_buffered_path() {
    let input = skewed_reads();
    let (estimate, log) = run(
        &input,
        &[
            "-T",
            "200",
            "-Q",
            "100",
            "--seed",
            "42",
            "--normalize",
            "auto",
            "--max-read-buffer",
            "4G",
        ],
    );

    assert!(estimate > 0);
    assert!(!log.contains("Selected reads by position"));
}

/// Reads are scored for retention in parallel, so a thread count has to be free to change without
/// moving the estimate. Both normalized selection paths score every read, so both are checked.
#[test]
fn the_estimate_does_not_depend_on_the_thread_count() {
    let input = skewed_reads();
    let arguments = [
        "-T",
        "200",
        "-Q",
        "100",
        "--seed",
        "42",
        "--normalize",
        "auto",
    ];

    for buffer in ["4G", "1"] {
        let mut single = arguments.to_vec();
        single.extend(["-t", "1", "--max-read-buffer", buffer]);
        let mut many = arguments.to_vec();
        many.extend(["-t", "4", "--max-read-buffer", buffer]);

        let (single_estimate, single_log) = run(&input, &single);
        let (many_estimate, _) = run(&input, &many);

        assert!(single_log.contains("Depth skew detected"));
        assert_eq!(
            single_log.contains("Selected reads by position"),
            buffer == "1"
        );
        assert_eq!(
            single_estimate, many_estimate,
            "at a {buffer} read buffer, one thread estimated {single_estimate}, four estimated {many_estimate}"
        );
    }
}

#[test]
fn normalization_corrects_a_high_copy_element_at_a_median_depth_of_one() {
    let input = low_depth_reads(LOW_DEPTH_ELEMENT_READS);
    let (legacy, _) = run(
        &input,
        &[LOW_DEPTH_ARGUMENTS.as_slice(), &["never"]].concat(),
    );
    let (normalized, normalized_log) = run(
        &input,
        &[LOW_DEPTH_ARGUMENTS.as_slice(), &["auto"]].concat(),
    );

    // The trailing space keeps this off a median depth of ten or nineteen.
    assert!(
        normalized_log.contains("Median depth is 1 "),
        "this input is meant to normalize against a median depth of one:\n{normalized_log}"
    );
    assert!(normalized_log.contains("Depth skew detected"));
    assert!(legacy < 20_000, "legacy estimate was {legacy}");
    assert!(
        (80_000..=180_000).contains(&normalized),
        "normalized estimate was {normalized}, against a chromosome of {LOW_DEPTH_CHROMOSOME_SIZE}"
    );
}

/// A median depth of one is a small number to divide the high percentile by, so the question is
/// whether an ordinary input gets called skewed down there by arithmetic alone. It does not.
#[test]
fn an_even_depth_input_is_not_normalized_at_a_median_depth_of_one() {
    let input = low_depth_reads(0);
    // Only a run that normalizes builds a profile, so the median this input would normalize
    // against has to be read off a forced one.
    let (_, forced_log) = run(
        &input,
        &[LOW_DEPTH_ARGUMENTS.as_slice(), &["always"]].concat(),
    );
    let (_, log) = run(
        &input,
        &[LOW_DEPTH_ARGUMENTS.as_slice(), &["auto"]].concat(),
    );

    // The trailing space keeps this off a median depth of ten or nineteen.
    assert!(
        forced_log.contains("Median depth is 1 "),
        "this input is meant to normalize against a median depth of one:\n{forced_log}"
    );
    assert!(
        log.contains("Depth skew not detected"),
        "an even-depth input was called skewed:\n{log}"
    );
}
