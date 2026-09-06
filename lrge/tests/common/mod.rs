//! Fixture building and process running shared by the integration tests.
//!
//! Each integration test file is its own binary, so a helper only one of them uses is dead code in
//! the others. The allow is what lets this module hold the whole shared set rather than the subset
//! every test happens to need.
#![allow(dead_code)]

use assert_cmd::Command;
use std::io::Write;
use tempfile::NamedTempFile;

/// One step of the generator the inputs are built from.
pub fn next_draw(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1);
    *state >> 32
}

pub fn pseudo_random_dna(len: usize, seed: u64) -> Vec<u8> {
    let mut state = seed;
    (0..len)
        .map(|_| b"ACGT"[(next_draw(&mut state) & 3) as usize])
        .collect()
}

/// A read of `len` bases from `source`, wrapping at its end.
pub fn circular_read(source: &[u8], start: usize, len: usize) -> Vec<u8> {
    (0..len)
        .map(|offset| source[(start + offset) % source.len()])
        .collect()
}

/// Write `count` reads cut from `source`, named `prefix0` onwards.
///
/// The starts step by a stride coprime with every source length here, so the reads walk the whole
/// of it rather than piling up on one stretch.
pub fn write_reads(
    input: &mut NamedTempFile,
    prefix: &str,
    source: &[u8],
    count: usize,
    mut cut: impl FnMut(&[u8], usize) -> Vec<u8>,
) {
    for index in 0..count {
        let read = cut(source, index * 137 % source.len());
        writeln!(
            input,
            ">{prefix}{index}\n{}",
            String::from_utf8(read).unwrap()
        )
        .unwrap();
    }
}

/// Run `lrge` over `input`, returning its estimate and its log.
pub fn run(input: &NamedTempFile, arguments: &[&str]) -> (u64, String) {
    let output = Command::cargo_bin("lrge")
        .unwrap()
        .arg(input.path())
        .args(arguments)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "lrge failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let estimate = String::from_utf8(output.stdout)
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    (estimate, String::from_utf8(output.stderr).unwrap())
}
