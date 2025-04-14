use clap::{builder::ArgPredicate, Parser};
use std::ffi::OsStr;
use std::path::PathBuf;

const TARGET_NUM_READS: &str = "10000";
const QUERY_NUM_READS: &str = "5000";
const MAX_OVERHANG_RATIO: &str = "0.2";

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    /// Input FASTQ file
    #[arg(name = "INPUT", value_parser = check_path_exists)]
    pub input: PathBuf,

    /// Output file for the estimate
    #[arg(short, long, value_name = "OUTPUT", default_value = "-")]
    pub output: String,

    /// Target number of reads to use (for two-set strategy; default)
    #[arg(short = 'T', long = "target", value_name = "INT", default_value_if("num_reads", ArgPredicate::IsPresent, None), default_value = TARGET_NUM_READS)]
    pub target_num_reads: Option<usize>,

    /// Query number of reads to use (for two-set strategy; default)
    #[arg(short = 'Q', long = "query", value_name = "INT", default_value_if("num_reads", ArgPredicate::IsPresent, None), default_value = QUERY_NUM_READS)]
    pub query_num_reads: Option<usize>,

    /// Number of reads to use (for all-vs-all strategy)
    #[arg(short, long = "num", value_name = "INT", conflicts_with_all = &["target_num_reads", "query_num_reads"])]
    pub num_reads: Option<usize>,

    /// Sequencing platform of the reads
    #[arg(short = 'P', long, value_name = "PLATFORM", value_parser = ["ont", "pb"], default_value = "ont")]
    pub platform: String,

    /// Exclude overlaps for internal matches
    #[arg(short = 'F', long = "filter-contained")]
    pub filter_contained: bool,

    /// Number of threads to use
    #[arg(short, long, value_name = "INT", default_value = "1")]
    pub threads: usize,

    /// Don't clean up temporary files
    #[arg(short = 'C', long)]
    pub keep_temp: bool,

    /// Temporary directory for storing intermediate files
    #[arg(short = 'D', long = "temp", value_name = "DIR")]
    pub temp_dir: Option<PathBuf>,

    /// Random seed to use - making the estimate repeatable
    #[clap(short = 's', long = "seed", value_name = "INT")]
    pub seed: Option<u64>,

    /// Take the estimate as the median of all estimates, *including infinite estimates*
    #[arg(short = '8', long = "inf", hide_short_help = true)]
    pub with_infinity: bool,

    /// I neeeeeed that precision! Output the estimate as a floating point number
    #[arg(short = 'f', long = "float-my-boat", hide_short_help = true)]
    pub precise: bool,

    /// The lower quantile to use for the estimate
    #[arg(long = "q1", value_name = "FLOAT", default_value_t = liblrge::estimate::LOWER_QUANTILE, value_parser = validate_low_quantile, hide_short_help = true)]
    pub lower_q: f32,

    /// The upper quantile to use for the estimate
    #[arg(long = "q3", value_name = "FLOAT", default_value_t = liblrge::estimate::UPPER_QUANTILE, value_parser = validate_high_quantile, hide_short_help = true)]
    pub upper_q: f32,

    /// Maximum overhang size to alignment length ratio for internal overlap filtering
    #[arg(long = "max-overhang-ratio", value_name = "FLOAT", default_value = MAX_OVERHANG_RATIO, value_parser = validate_overhang_ratio, hide_short_help = true)]
    pub max_overhang_ratio: f32,

    /// Use the smaller Q/T dataset as minimap2 reference (for two-set strategy)
    #[arg(long = "use-min-ref", hide_short_help = true)]
    pub use_min_ref: bool,

    /// `-q` only show errors and warnings. `-qq` only show errors. `-qqq` shows nothing.
    #[arg(short, long, action = clap::ArgAction::Count, conflicts_with = "verbose")]
    pub quiet: u8,

    /// `-v` show debug output. `-vv` show trace output.
    #[arg(short, long, action = clap::ArgAction::Count)]
    pub verbose: u8,
}

/// A utility function that allows the CLI to error if a path doesn't exist
fn check_path_exists<S: AsRef<OsStr> + ?Sized>(s: &S) -> Result<PathBuf, String> {
    let path = PathBuf::from(s);
    if path.exists() {
        Ok(path)
    } else {
        Err(format!("{} does not exist", path.to_string_lossy()))
    }
}

/// A generic value parser to ensure the value is within the specified range
fn validate_quantile(s: &str, min: f32, max: f32) -> Result<f32, String> {
    let value: f32 = s
        .parse()
        .map_err(|_| format!("`{}` is not a valid number", s))?;
    if value > min && value < max {
        Ok(value)
    } else {
        Err(format!(
            "Value `{}` must be greater than {} and less than {}",
            s, min, max
        ))
    }
}

/// A value parser for the lower quantile
fn validate_low_quantile(s: &str) -> Result<f32, String> {
    validate_quantile(s, 0.0, 0.5)
}

/// A value parser for the upper quantile
fn validate_high_quantile(s: &str) -> Result<f32, String> {
    validate_quantile(s, 0.5, 1.0)
}

/// A value parser for the maximum overhang ratio
fn validate_overhang_ratio(s: &str) -> Result<f32, String> {
    let value: f32 = s
        .parse()
        .map_err(|_| format!("`{}` is not a valid number", s))?;
    if (0.0..=1.0).contains(&value) {
        Ok(value)
    } else {
        Err(format!("Value `{}` must be between 0.0 and 1.0", s))
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
    fn test_validate_quantile() {
        assert!(validate_quantile("0.1", 0.0, 0.5).is_ok());
        assert!(validate_quantile("0.5", 0.0, 0.5).is_err());
        assert!(validate_quantile("0", 0.0, 0.5).is_err());
        assert!(validate_quantile("-0.1", 0.0, 0.5).is_err());
        assert!(validate_quantile("abc", 0.0, 0.5).is_err());
        assert!(validate_quantile("0.6", 0.5, 1.0).is_ok());
        assert!(validate_quantile("1.0", 0.5, 1.0).is_err());
    }

    #[test]
    fn cli_no_args() {
        let opts = Args::try_parse_from([BIN]);
        assert!(opts.is_err());
        assert!(opts
            .unwrap_err()
            .to_string()
            .contains("error: the following required arguments were not provided"));
    }

    #[test]
    fn cli_with_input() {
        let opts = Args::try_parse_from([BIN, "Cargo.toml"]).unwrap();

        assert_eq!(opts.input, PathBuf::from("Cargo.toml"));
        assert_eq!(
            opts.target_num_reads,
            Some(TARGET_NUM_READS.parse().unwrap())
        );
        assert_eq!(opts.query_num_reads, Some(QUERY_NUM_READS.parse().unwrap()));
    }

    #[test]
    fn cli_with_num_reads() {
        let opts = Args::try_parse_from([BIN, "Cargo.toml", "--num", "100"]).unwrap();

        assert_eq!(opts.input, PathBuf::from("Cargo.toml"));
        assert_eq!(opts.num_reads, Some(100));
        assert_eq!(opts.target_num_reads, None);
        assert_eq!(opts.query_num_reads, None);
    }

    #[test]
    fn cli_with_target_and_query_reads() {
        let opts =
            Args::try_parse_from([BIN, "Cargo.toml", "--target", "100", "--query", "200"]).unwrap();
        assert_eq!(opts.input, PathBuf::from("Cargo.toml"));
        assert_eq!(opts.num_reads, None);
        assert_eq!(opts.target_num_reads, Some(100));
        assert_eq!(opts.query_num_reads, Some(200));
    }

    #[test]
    fn cli_with_num_reads_and_target_reads_and_query_reads() {
        let opts = Args::try_parse_from([
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
            .contains("error: the argument '--num <INT>' cannot be used with"));
    }

    #[test]
    fn cli_with_num_reads_and_target_reads() {
        let opts = Args::try_parse_from([BIN, "Cargo.toml", "--num", "100", "--target", "200"]);
        assert!(opts.is_err());
        assert!(opts
            .unwrap_err()
            .to_string()
            .contains("error: the argument '--num <INT>' cannot be used with"));
    }

    #[test]
    fn cli_with_num_reads_and_query_reads() {
        let opts = Args::try_parse_from([BIN, "Cargo.toml", "--num", "100", "--query", "200"]);
        assert!(opts.is_err());
        assert!(opts
            .unwrap_err()
            .to_string()
            .contains("error: the argument '--num <INT>' cannot be used with"));
    }

    #[test]
    fn cli_with_target_reads_no_query_reads() {
        let opts = Args::try_parse_from([BIN, "Cargo.toml", "--target", "100"]).unwrap();
        assert_eq!(opts.target_num_reads, Some(100));
        assert_eq!(opts.query_num_reads, Some(QUERY_NUM_READS.parse().unwrap()));
    }

    #[test]
    fn cli_with_query_reads_no_target_reads() {
        let opts = Args::try_parse_from([BIN, "Cargo.toml", "--query", "100"]).unwrap();
        assert_eq!(opts.query_num_reads, Some(100));
        assert_eq!(
            opts.target_num_reads,
            Some(TARGET_NUM_READS.parse().unwrap())
        );
    }

    #[test]
    fn cli_with_quiet() {
        let opts = Args::try_parse_from([BIN, "Cargo.toml", "-q"]).unwrap();
        assert_eq!(opts.quiet, 1);
    }

    #[test]
    fn cli_with_verbose() {
        let opts = Args::try_parse_from([BIN, "Cargo.toml", "-v"]).unwrap();
        assert_eq!(opts.verbose, 1);
    }

    #[test]
    fn cli_with_verbose_verbose() {
        let opts = Args::try_parse_from([BIN, "Cargo.toml", "-vv"]).unwrap();
        assert_eq!(opts.verbose, 2);
    }

    #[test]
    fn cli_with_verbose_verbose_verbose() {
        let opts = Args::try_parse_from([BIN, "Cargo.toml", "-vvv"]).unwrap();
        assert_eq!(opts.verbose, 3);
    }

    #[test]
    fn cli_with_quiet_verbose() {
        let opts = Args::try_parse_from([BIN, "Cargo.toml", "-qv"]);
        assert!(opts.is_err());
        assert!(opts
            .unwrap_err()
            .to_string()
            .contains("error: the argument '--quiet...' cannot be used with"));
    }
}
