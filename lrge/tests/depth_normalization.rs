use assert_cmd::Command;
use std::io::Write;
use tempfile::NamedTempFile;

const CHROMOSOME_SIZE: usize = 20_000;
const ELEMENT_SIZE: usize = 2_000;
const READ_LENGTH: usize = 800;

fn pseudo_random_dna(len: usize, seed: u64) -> Vec<u8> {
    let mut state = seed;
    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            b"ACGT"[((state >> 32) & 3) as usize]
        })
        .collect()
}

fn circular_read(sequence: &[u8], start: usize) -> Vec<u8> {
    (0..READ_LENGTH)
        .map(|offset| sequence[(start + offset) % sequence.len()])
        .collect()
}

fn skewed_reads() -> NamedTempFile {
    let chromosome = pseudo_random_dna(CHROMOSOME_SIZE, 1);
    let element = pseudo_random_dna(ELEMENT_SIZE, 2);
    let mut input = NamedTempFile::new().unwrap();

    for index in 0..400 {
        let read = circular_read(&chromosome, index * 137 % chromosome.len());
        writeln!(
            input,
            ">chromosome{index}\n{}",
            String::from_utf8(read).unwrap()
        )
        .unwrap();
    }
    for index in 0..4_000 {
        let read = circular_read(&element, index * 137 % element.len());
        writeln!(
            input,
            ">element{index}\n{}",
            String::from_utf8(read).unwrap()
        )
        .unwrap();
    }

    input
}

fn run(input: &NamedTempFile, arguments: &[&str]) -> (u64, String) {
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
