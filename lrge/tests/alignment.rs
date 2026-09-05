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
        // What this test is for is that a BAM reaches the estimator. Five hundred reads cannot
        // size a genome: most of the fifty query reads find no overlap among the hundred target
        // reads, so the estimate rests on a handful of them and swings by seed, which is why the
        // seed is pinned. The fixture is also small enough that detection samples every read of it
        // and calls it skewed, so it normalizes against a profile median depth of one. That regime
        // has since been measured on inputs whose genome size is known, and normalizing there
        // lands closer to the truth, so this runs the default path instead of stepping around it.
        // See `paper/corrections/issue56_low_depth_normalization.tsv`.
        cmd.arg(bam_path)
            .arg("-T")
            .arg("100")
            .arg("-Q")
            .arg("50")
            .arg("--seed")
            .arg("6")
            .assert()
            .success();
    }
}
