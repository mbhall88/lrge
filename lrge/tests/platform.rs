//! `-P` picks the minimap2 overlap preset, and the log says which one a run used.
//!
//! The flag reached the aligner from the library but not from the command line: `main` parsed it,
//! printed it, and never handed it to a builder, so every release up to 1.0.0 overlapped PacBio
//! reads with the Nanopore preset. These run the whole way through rather than the argument
//! parsing, because argument parsing was never what was broken.

mod common;

use common::{circular_read, pseudo_random_dna, run, write_reads};
use tempfile::NamedTempFile;

const GENOME_SIZE: usize = 20_000;
const READ_LENGTH: usize = 800;
const READS: usize = 400;
const READ_SETS: [&str; 4] = ["-T", "300", "-Q", "100"];

fn reads() -> NamedTempFile {
    let genome = pseudo_random_dna(GENOME_SIZE, 7);
    let mut input = NamedTempFile::new().unwrap();
    write_reads(&mut input, "read", &genome, READS, |source, start| {
        circular_read(source, start, READ_LENGTH)
    });
    input
}

#[test]
fn the_platform_picks_the_overlap_preset() {
    let input = reads();

    for (platform, preset) in [("ont", "ava-ont"), ("pb", "ava-pb")] {
        let mut arguments = vec!["-v", "-P", platform];
        arguments.extend(READ_SETS);
        let (_, log) = run(&input, &arguments);

        assert!(
            log.contains(&format!("with the {preset} preset")),
            "-P {platform} should overlap with {preset}:\n{log}"
        );
    }
}

#[test]
fn nanopore_is_the_platform_a_run_that_says_nothing_gets() {
    let input = reads();

    let mut arguments = vec!["-v"];
    arguments.extend(READ_SETS);
    let (_, log) = run(&input, &arguments);

    assert!(log.contains("with the ava-ont preset"), "{log}");
}

#[test]
fn the_all_vs_all_strategy_takes_the_platform_too() {
    let input = reads();
    let (_, log) = run(&input, &["-v", "-n", "200", "-P", "pb"]);

    assert!(log.contains("with the ava-pb preset"), "{log}");
}
