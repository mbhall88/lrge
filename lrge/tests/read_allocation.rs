mod common;

use common::{pseudo_random_dna, run, write_reads};
use tempfile::NamedTempFile;

const GENOME_SIZE: usize = 20_000;
const READ_LENGTH: usize = 800;
/// Twelve-fold coverage, which is thin but enough that most query reads find a target to overlap.
const READ_COUNT: usize = 300;

/// Reads of one small genome, at a count deliberately too low to fill the requests below.
fn small_input() -> NamedTempFile {
    let genome = pseudo_random_dna(GENOME_SIZE, 1);
    let mut input = NamedTempFile::new().unwrap();

    write_reads(&mut input, "read", &genome, READ_COUNT, |source, start| {
        common::circular_read(source, start, READ_LENGTH)
    });

    input
}

/// Normalization is off throughout so that what reaches each set is the allocation rule alone.
const ARGUMENTS: [&str; 4] = ["--seed", "42", "--normalize", "never"];

#[test]
fn an_input_too_small_for_both_sets_keeps_the_requested_ratio() {
    let input = small_input();

    let (_, log) = run(
        &input,
        &[ARGUMENTS.as_slice(), &["-T", "200", "-Q", "150"]].concat(),
    );

    // 300 reads split in the requested 4:3, rather than 150 and 150.
    assert!(
        log.contains("Using 172 target reads and 128 query reads"),
        "log was:\n{log}"
    );
    // Fifty reads short of 350 is not the regime the second warning is for.
    assert!(
        !log.contains("less than half the reads requested"),
        "log was:\n{log}"
    );
}

#[test]
fn the_target_mode_keeps_the_full_query_request() {
    let input = small_input();

    let (_, log) = run(
        &input,
        &[
            ARGUMENTS.as_slice(),
            &["-T", "200", "-Q", "150", "--shortfall", "target"],
        ]
        .concat(),
    );

    assert!(
        log.contains("Using 150 target reads and 150 query reads"),
        "log was:\n{log}"
    );
}

/// An input smaller than the query request used to be a hard error, so LRGE refused to run on any
/// input of 5,000 reads or fewer at the defaults.
#[test]
fn an_input_below_the_query_request_still_estimates() {
    let input = small_input();

    let (estimate, log) = run(&input, ARGUMENTS.as_slice());

    assert!(
        log.contains("Using 200 target reads and 100 query reads"),
        "log was:\n{log}"
    );
    // Three hundred reads against a request for fifteen thousand is the regime that warrants it.
    assert!(
        log.contains("less than half the reads requested"),
        "log was:\n{log}"
    );
    assert!(estimate > 0, "estimate was {estimate}");
}
