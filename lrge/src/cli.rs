use clap::{builder::ArgPredicate, Parser};
use std::ffi::OsStr;
use std::path::PathBuf;

const TARGET_NUM_READS: &str = "10000";
const QUERY_NUM_READS: &str = "5000";
const MAX_OVERHANG_RATIO: &str = "0.2";
const MAX_READ_BUFFER: &str = "1G";

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    /// Input FASTQ, FASTA, or unaligned BAM/CRAM/SAM file
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

    /// Control depth-aware read normalization
    #[arg(long, value_name = "MODE", default_value = "auto")]
    pub normalize: liblrge::Normalization,

    /// How to split an input too small to supply both read sets [scale, target]
    #[arg(long, value_name = "MODE", default_value = "scale")]
    pub shortfall: liblrge::Shortfall,

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
    ///
    /// Only meaningful alongside -F/--filter-contained, which this option requires.
    #[arg(long = "max-overhang-ratio", value_name = "FLOAT", default_value = MAX_OVERHANG_RATIO, value_parser = validate_overhang_ratio, requires = "filter_contained", hide_short_help = true)]
    pub max_overhang_ratio: f32,

    /// Use the smaller Q/T dataset as minimap2 reference (for two-set strategy)
    #[arg(long = "use-min-ref", hide_short_help = true)]
    pub use_min_ref: bool,

    /// Cap on the memory used to buffer selected reads when normalizing (e.g. 512M, 1.5G)
    ///
    /// Above this, lrge buffers read positions and reads the input one extra time. The reads
    /// selected for a given seed are the same either way.
    #[arg(long = "max-read-buffer", value_name = "SIZE", default_value = MAX_READ_BUFFER, value_parser = parse_byte_size, hide_short_help = true)]
    pub max_read_buffer: u64,

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
        .map_err(|_| format!("`{s}` is not a valid number",))?;
    if value > min && value < max {
        Ok(value)
    } else {
        Err(format!(
            "Value `{s}` must be greater than {min} and less than {max}",
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

/// A value parser for a memory size, in bytes or with a binary suffix such as `512M` or `1.5G`
///
/// A fractional size is scaled by its suffix and truncated to whole bytes, so `1.5G` is
/// 1610612736 and `0.5K` is 512.
fn parse_byte_size(s: &str) -> Result<u64, String> {
    let text = s.trim();
    let split = text
        .find(|c: char| !c.is_ascii_digit() && c != '.')
        .unwrap_or(text.len());
    let (number, suffix) = text.split_at(split);
    let invalid = || format!("`{s}` is not a valid size");
    let too_large = || format!("Size `{s}` is too large");

    if number.is_empty() {
        return Err(invalid());
    }
    let multiplier = 1_u128
        << parse_size_suffix(suffix)
            .ok_or_else(|| format!("`{s}` has an unknown size suffix; use K, M, G, or T"))?;

    let (whole, fraction) = number.split_once('.').unwrap_or((number, ""));
    if whole.is_empty() && fraction.is_empty() {
        return Err(invalid());
    }

    let whole: u128 = match whole {
        "" => 0,
        digits => digits.parse().map_err(|_| too_large())?,
    };
    let mut bytes = whole.checked_mul(multiplier).ok_or_else(too_large)?;

    if !fraction.is_empty() {
        // Scale the fraction with integers so a size like 1.5G stays exact.
        let digits: u128 = fraction.parse().map_err(|_| invalid())?;
        let divisor = 10_u128
            .checked_pow(u32::try_from(fraction.len()).map_err(|_| invalid())?)
            .ok_or_else(invalid)?;
        let scaled = digits.checked_mul(multiplier).ok_or_else(too_large)? / divisor;
        bytes = bytes.checked_add(scaled).ok_or_else(too_large)?;
    }

    u64::try_from(bytes).map_err(|_| too_large())
}

/// The power of two a size suffix stands for, ignoring an optional `B`, `iB`, and case
fn parse_size_suffix(suffix: &str) -> Option<u32> {
    let suffix = suffix.trim_start();
    let suffix = suffix
        .strip_suffix("iB")
        .or_else(|| suffix.strip_suffix("IB"))
        .or_else(|| suffix.strip_suffix('B'))
        .or_else(|| suffix.strip_suffix('b'))
        .unwrap_or(suffix);
    match suffix.to_ascii_uppercase().as_str() {
        "" => Some(0),
        "K" => Some(10),
        "M" => Some(20),
        "G" => Some(30),
        "T" => Some(40),
        _ => None,
    }
}

/// A value parser for the maximum overhang ratio
fn validate_overhang_ratio(s: &str) -> Result<f32, String> {
    let value: f32 = s
        .parse()
        .map_err(|_| format!("`{s}` is not a valid number",))?;
    if (0.0..=1.0).contains(&value) {
        Ok(value)
    } else {
        Err(format!("Value `{s}` must be between 0.0 and 1.0",))
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
    fn byte_size_accepts_plain_bytes() {
        assert_eq!(parse_byte_size("1024").unwrap(), 1024);
    }

    #[test]
    fn byte_size_accepts_binary_suffixes() {
        for (text, expected) in [
            ("2K", 2 << 10),
            ("2k", 2 << 10),
            ("2KB", 2 << 10),
            ("2KiB", 2 << 10),
            ("512M", 512 << 20),
            ("4G", 4u64 << 30),
            ("1T", 1u64 << 40),
        ] {
            assert_eq!(parse_byte_size(text).unwrap(), expected, "parsing {text}");
        }
    }

    #[test]
    fn byte_size_allows_a_space_before_the_suffix() {
        assert_eq!(parse_byte_size("4 G").unwrap(), 4u64 << 30);
    }

    #[test]
    fn byte_size_accepts_fractional_sizes() {
        for (text, expected) in [
            ("1.5G", 1_610_612_736),
            ("0.5K", 512),
            (".5K", 512),
            ("1.25G", 1_342_177_280),
            ("2.5MiB", 2_621_440),
            ("1.0G", 1 << 30),
        ] {
            assert_eq!(parse_byte_size(text).unwrap(), expected, "parsing {text}");
        }
    }

    #[test]
    fn byte_size_truncates_a_fraction_of_a_byte() {
        assert_eq!(parse_byte_size("1.5").unwrap(), 1);
    }

    #[test]
    fn byte_size_rejects_a_value_that_overflows_its_suffix() {
        // a shift that fits in u64 can still shift every bit of the value out
        assert!(parse_byte_size("16777216T").is_err());
        assert!(parse_byte_size("18446744073709551615G").is_err());
    }

    #[test]
    fn byte_size_rejects_nonsense() {
        for text in ["", "G", "-1", ".", "1.2.3", "4X", "four", "1.5X"] {
            assert!(parse_byte_size(text).is_err(), "accepted {text}");
        }
    }

    #[test]
    fn cli_defaults_the_read_buffer_cap() {
        let opts = Args::try_parse_from([BIN, "Cargo.toml"]).unwrap();

        assert_eq!(opts.max_read_buffer, liblrge::DEFAULT_MAX_READ_BUFFER);
    }

    #[test]
    fn cli_accepts_a_read_buffer_cap() {
        let opts = Args::try_parse_from([BIN, "Cargo.toml", "--max-read-buffer", "512M"]).unwrap();

        assert_eq!(opts.max_read_buffer, 512 << 20);
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
        assert_eq!(opts.normalize, liblrge::Normalization::Auto);
        assert_eq!(
            opts.target_num_reads,
            Some(TARGET_NUM_READS.parse().unwrap())
        );
        assert_eq!(opts.query_num_reads, Some(QUERY_NUM_READS.parse().unwrap()));
    }

    #[test]
    fn cli_accepts_normalization_modes() {
        for (mode, expected) in [
            ("auto", liblrge::Normalization::Auto),
            ("always", liblrge::Normalization::Always),
            ("never", liblrge::Normalization::Never),
        ] {
            let opts = Args::try_parse_from([BIN, "Cargo.toml", "--normalize", mode]).unwrap();
            assert_eq!(opts.normalize, expected);
        }
    }

    #[test]
    fn cli_rejects_invalid_normalization_mode() {
        let error = Args::try_parse_from([BIN, "Cargo.toml", "--normalize", "sometimes"])
            .unwrap_err()
            .to_string();

        assert!(error.contains("expected auto, always, or never"));
    }

    #[test]
    fn cli_accepts_shortfall_modes() {
        for (mode, expected) in [
            ("scale", liblrge::Shortfall::Scale),
            ("target", liblrge::Shortfall::Target),
        ] {
            let opts = Args::try_parse_from([BIN, "Cargo.toml", "--shortfall", mode]).unwrap();
            assert_eq!(opts.shortfall, expected);
        }
    }

    #[test]
    fn cli_defaults_to_scaling_the_shortfall() {
        let opts = Args::try_parse_from([BIN, "Cargo.toml"]).unwrap();
        assert_eq!(opts.shortfall, liblrge::Shortfall::Scale);
    }

    #[test]
    fn cli_rejects_invalid_shortfall_mode() {
        let error = Args::try_parse_from([BIN, "Cargo.toml", "--shortfall", "query"])
            .unwrap_err()
            .to_string();

        assert!(error.contains("expected scale or target"));
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

    #[test]
    fn cli_max_overhang_ratio_without_filter_contained_is_an_error() {
        let opts = Args::try_parse_from([BIN, "Cargo.toml", "--max-overhang-ratio", "0.05"]);

        let err = opts.unwrap_err().to_string();
        assert!(
            err.contains("--filter-contained"),
            "error should name the flag that is missing, got: {err}"
        );
    }

    #[test]
    fn cli_max_overhang_ratio_with_filter_contained_is_accepted() {
        let opts = Args::try_parse_from([BIN, "Cargo.toml", "-F", "--max-overhang-ratio", "0.05"])
            .unwrap();

        assert!(opts.filter_contained);
        assert_eq!(opts.max_overhang_ratio, 0.05);
    }

    #[test]
    fn cli_default_overhang_ratio_matches_the_library_default() {
        // the CLI has to spell the default as a string for clap; keep it in step with the
        // single definition in liblrge rather than letting the two drift
        assert_eq!(
            MAX_OVERHANG_RATIO.parse::<f32>().unwrap(),
            liblrge::DEFAULT_MAX_OVERHANG_RATIO
        );
    }

    #[test]
    fn cli_max_overhang_ratio_default_does_not_require_filter_contained() {
        // the default value must not trip the requirement - a plain run has to keep working
        let opts = Args::try_parse_from([BIN, "Cargo.toml"]).unwrap();

        assert!(!opts.filter_contained);
        assert_eq!(opts.max_overhang_ratio, MAX_OVERHANG_RATIO.parse().unwrap());
    }
}
