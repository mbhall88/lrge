use clap::{builder::ArgPredicate, Parser};
use std::ffi::OsStr;
use std::path::PathBuf;

const TARGET_NUM_READS: &str = "5000";
const QUERY_NUM_READS: &str = "10000";

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Opts {
    /// Input FASTQ file
    #[arg(name = "INPUT", value_parser = check_path_exists)]
    input: PathBuf,

    /// Number of reads to use (for all-vs-all strategy)
    #[arg(short, long = "num", value_name = "NUM", conflicts_with_all = &["target_num_reads", "query_num_reads"])]
    num_reads: Option<usize>,

    /// Target number of reads to use (for two-set strategy; default)
    #[arg(short = 'T', long = "target", value_name = "TARGET", default_value_if("num_reads", ArgPredicate::IsPresent, None), default_value = TARGET_NUM_READS)]
    target_num_reads: Option<usize>,

    /// Query number of reads to use (for two-set strategy; default)
    #[arg(short = 'Q', long = "query", value_name = "QUERY", default_value_if("num_reads", ArgPredicate::IsPresent, None), default_value = QUERY_NUM_READS)]
    query_num_reads: Option<usize>,
}

/// A utility function that allows the CLI to error if a path doesn't exist
fn check_path_exists<S: AsRef<OsStr> + ?Sized>(s: &S) -> Result<PathBuf, String> {
    let path = PathBuf::from(s);
    if path.exists() {
        Ok(path)
    } else {
        Err(format!("{:?} does not exist", path))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const BIN: &str = env!("CARGO_BIN_NAME");
    #[test]
    fn check_path_exists_it_doesnt() {
        let result = check_path_exists(OsStr::new("fake.path"));
        assert!(result.is_err())
    }

    #[test]
    fn check_path_it_does() {
        let actual = check_path_exists(OsStr::new("Cargo.toml")).unwrap();
        let expected = PathBuf::from("Cargo.toml");
        assert_eq!(actual, expected)
    }

    #[test]
    fn cli_no_args() {
        let opts = Opts::try_parse_from([BIN]);
        assert!(opts.is_err());
        assert!(opts
            .unwrap_err()
            .to_string()
            .contains("error: the following required arguments were not provided"));
    }

    #[test]
    fn cli_with_input() {
        let opts = Opts::try_parse_from([BIN, "Cargo.toml"]).unwrap();

        assert_eq!(opts.input, PathBuf::from("Cargo.toml"));
        assert_eq!(
            opts.target_num_reads,
            Some(TARGET_NUM_READS.parse().unwrap())
        );
        assert_eq!(opts.query_num_reads, Some(QUERY_NUM_READS.parse().unwrap()));
    }

    #[test]
    fn cli_with_num_reads() {
        let opts = Opts::try_parse_from([BIN, "Cargo.toml", "--num", "100"]).unwrap();

        assert_eq!(opts.input, PathBuf::from("Cargo.toml"));
        assert_eq!(opts.num_reads, Some(100));
        assert_eq!(opts.target_num_reads, None);
        assert_eq!(opts.query_num_reads, None);
    }

    #[test]
    fn cli_with_target_and_query_reads() {
        let opts =
            Opts::try_parse_from([BIN, "Cargo.toml", "--target", "100", "--query", "200"]).unwrap();
        assert_eq!(opts.input, PathBuf::from("Cargo.toml"));
        assert_eq!(opts.num_reads, None);
        assert_eq!(opts.target_num_reads, Some(100));
        assert_eq!(opts.query_num_reads, Some(200));
    }

    #[test]
    fn cli_with_num_reads_and_target_reads_and_query_reads() {
        let opts = Opts::try_parse_from([
            BIN,
            "Cargo.toml",
            "--num",
            "100",
            "--target",
            "200",
            "--query",
            "300",
        ]);
        assert!(opts.is_err());
        assert!(opts
            .unwrap_err()
            .to_string()
            .contains("error: the argument '--num <NUM>' cannot be used with"));
    }

    #[test]
    fn cli_with_num_reads_and_target_reads() {
        let opts = Opts::try_parse_from([BIN, "Cargo.toml", "--num", "100", "--target", "200"]);
        assert!(opts.is_err());
        assert!(opts
            .unwrap_err()
            .to_string()
            .contains("error: the argument '--num <NUM>' cannot be used with"));
    }

    #[test]
    fn cli_with_num_reads_and_query_reads() {
        let opts = Opts::try_parse_from([BIN, "Cargo.toml", "--num", "100", "--query", "200"]);
        assert!(opts.is_err());
        assert!(opts
            .unwrap_err()
            .to_string()
            .contains("error: the argument '--num <NUM>' cannot be used with"));
    }

    #[test]
    fn cli_with_target_reads_no_query_reads() {
        let opts = Opts::try_parse_from([BIN, "Cargo.toml", "--target", "100"]).unwrap();
        assert_eq!(opts.target_num_reads, Some(100));
        assert_eq!(opts.query_num_reads, Some(QUERY_NUM_READS.parse().unwrap()));
    }

    #[test]
    fn cli_with_query_reads_no_target_reads() {
        let opts = Opts::try_parse_from([BIN, "Cargo.toml", "--query", "100"]).unwrap();
        assert_eq!(opts.query_num_reads, Some(100));
        assert_eq!(
            opts.target_num_reads,
            Some(TARGET_NUM_READS.parse().unwrap())
        );
    }
}
