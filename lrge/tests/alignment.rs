use assert_cmd::Command;
use std::io::Write;
use tempfile::NamedTempFile;

#[test]
fn test_sam_input() {
    let mut temp_file = NamedTempFile::new().unwrap();
    writeln!(
        temp_file,
        "@HD\tVN:1.6\tSO:unsorted\nREAD1\t4\t*\t0\t0\t*\t*\t0\t0\tGATTACA\t!!!!!!!\nREAD2\t4\t*\t0\t0\t*\t*\t0\t0\tGATTACA\t!!!!!!!"
    )
    .unwrap();

    let mut cmd = Command::cargo_bin("lrge").unwrap();
    // We expect failure because 2 reads won't generate a finite estimate,
    // but it should NOT fail due to parsing.
    cmd.arg(temp_file.path())
        .arg("-T")
        .arg("1")
        .arg("-Q")
        .arg("1")
        .assert()
        .failure()
        .stderr(predicates::str::contains(
            "No finite estimates were generated",
        ));
}

#[test]
fn test_mapped_sam_fails() {
    let mut temp_file = NamedTempFile::new().unwrap();
    writeln!(
        temp_file,
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\nREAD1\t0\tchr1\t1\t0\t7M\t*\t0\t0\tGATTACA\t!!!!!!!\nREAD2\t4\t*\t0\t0\t*\t*\t0\t0\tGATTACA\t!!!!!!!"
    )
    .unwrap();

    let mut cmd = Command::cargo_bin("lrge").unwrap();
    cmd.arg(temp_file.path())
        .arg("-T")
        .arg("1")
        .arg("-Q")
        .arg("1")
        .assert()
        .failure()
        .stderr(predicates::str::contains(
            "Mapped records are not supported",
        ));
}

#[test]
fn test_toy_bam_input() {
    let mut cmd = Command::cargo_bin("lrge").unwrap();
    let bam_path = std::path::Path::new("tests").join("data").join("toy.bam");

    if bam_path.exists() {
        // Use a fixed seed to ensure deterministic behavior across platforms
        cmd.arg(bam_path)
            .arg("-T")
            .arg("10")
            .arg("-Q")
            .arg("5")
            .arg("--seed")
            .arg("1")
            .assert()
            .success();
    }
}
