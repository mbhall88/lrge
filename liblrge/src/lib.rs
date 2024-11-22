//! # liblrge
//!
//! `liblrge` is a Rust library that provides utilities for estimating genome size for a given set
//! of reads.
//!
//! You can find a command-line interface (CLI) tool that uses this library in the [`lrge`][lrge] crate.
//!
//! [lrge]: https://crates.io/crates/lrge
//!
//! ## Usage
//!
//! The library provides two strategies for estimating genome size:
//!
//! ### [`TwoSetStrategy`]
//!
//! The two-set strategy uses two (random) sets of reads to estimate the genome size. The query set, which is
//! generally smaller, is overlapped against a target set of reads. A genome size estimate is generated
//! for each read in the query set, based on the number of overlaps and the average read length.
//! The median of these estimates is taken as the final genome size estimate.
//!
//! ```no_run
//! use liblrge::{Estimate, TwoSetStrategy};
//! use liblrge::twoset::{Builder, DEFAULT_TARGET_NUM_READS, DEFAULT_QUERY_NUM_READS};
//!
//! let input = "path/to/reads.fastq";
//! let mut strategy = Builder::new()
//!    .target_num_reads(DEFAULT_TARGET_NUM_READS)
//!    .query_num_reads(DEFAULT_QUERY_NUM_READS)
//!    .threads(4)
//!    .build(input);
//!
//! let est_result = strategy.estimate(false, None, None).expect("Failed to generate estimate");
//! let estimate = est_result.estimate;
//! // do something with the estimate
//! ```
//!
//! ### [`AvaStrategy`]
//!
//! The all-vs-all (ava) strategy takes a (random) set of reads and overlaps it against itself to
//! estimate the genome size. The genome size estimate is generated for each read in the set, based on the
//! number of overlaps and the average read length - minus the read being assessed. The median of these
//! estimates is taken as the final genome size estimate.
//!
//! ```no_run
//! use liblrge::{Estimate, AvaStrategy};
//! use liblrge::ava::{Builder, DEFAULT_AVA_NUM_READS};
//!
//! let input = "path/to/reads.fastq";
//! let mut strategy = Builder::new()
//!    .num_reads(DEFAULT_AVA_NUM_READS)
//!   .threads(4)
//!   .build(input);
//!
//! let est_result = strategy.estimate(false, None, None).expect("Failed to generate estimate");
//! let estimate = est_result.estimate;
//! // do something with the estimate
//! ```
//!
//! ## Features
//!
//! This library includes optional support for compressed file formats, controlled by feature flags.
//! By default, the `compression` feature is enabled, which activates support for all included
//! compression formats.
//!
//! ### Available Features
//!
//! - **compression** (default): Enables all available compression formats (`gzip`, `zstd`, `bzip2`, `xz`).
//! - **gzip**: Enables support for gzip-compressed files (`.gz`) using the [`flate2`][flate2] crate.
//! - **zstd**: Enables support for zstd-compressed files (`.zst`) using the [`zstd`][zstd] crate.
//! - **bzip2**: Enables support for bzip2-compressed files (`.bz2`) using the [`bzip2`][bzip2] crate.
//! - **xz**: Enables support for xz-compressed files (`.xz`) using the [`liblzma`][xz] crate.
//!
//! ### Enabling and Disabling Features
//!
//! By default, all compression features are enabled. However, you can selectively enable or disable them
//! in your `Cargo.toml` to reduce dependencies or target specific compression formats:
//!
//! To **disable all compression features**:
//!
//! ```toml
//! liblrge = { version = "0.1.0", default-features = false }
//! ```
//!
//! To enable only specific compression formats, list the desired features in `Cargo.toml`:
//!
//! ```toml
//! liblrge = { version = "0.1.0", default-features = false, features = ["gzip", "zstd"] }
//! ```
//!
//! In this example, only `gzip` (`flate2`) and `zstd` are enabled, so `liblrge` will support `.gz`
//! and `.zst` files.
//!
//! ## Compression Detection
//!
//! The library uses [**magic bytes**][magic] at the start of the file to detect its compression
//! format before deciding how to read it. Supported formats include gzip, zstd, bzip2, and xz, with
//! automatic decompression if the [appropriate feature](#features) is enabled.
//!
//! [flate2]: https://crates.io/crates/flate2
//! [zstd]: https://crates.io/crates/zstd
//! [xz]: https://crates.io/liblzma
//! [bzip2]: https://crates.io/crates/bzip2
//! [magic]: https://en.wikipedia.org/wiki/Magic_number_(programming)#In_files
//!
//! ## Disabling logging
//!
//! `liblrge` will output some logging information via the [`log`][log] crate. If you wish to
//! suppress this logging you can configure the logging level in your application. For example, using
//! the [`env_logger`][env_logger] crate you can do the following:
//!
//! ```
//! use log::LevelFilter;
//!
//! let mut log_builder = env_logger::Builder::new();
//! log_builder
//!     .filter(None, LevelFilter::Info)
//!     .filter_module("liblrge", LevelFilter::Off);
//! log_builder.init();
//!
//! // Your application code here
//! ```
//!
//! This will set the global logging level to `Info` and disable all logging from the `liblrge` library.
//!
//! [log]: https://crates.io/crates/log
//! [env_logger]: https://crates.io/crates/env_logger
// todo add link to paper
#[deny(missing_docs)]
pub mod ava;
pub mod error;
pub mod estimate;
pub(crate) mod io;
pub(crate) mod minimap2;
pub mod twoset;

use rand::rngs::StdRng;
use rand::SeedableRng;

pub use self::ava::AvaStrategy;
pub use self::estimate::Estimate;
pub use self::twoset::TwoSetStrategy;
use std::str::FromStr;

/// A type alias for `Result` with [`LrgeError`][crate::error::LrgeError] as the error type.
pub type Result<T> = std::result::Result<T, error::LrgeError>;

/// The sequencing platform used to generate the reads.
///
/// # Examples
///
/// ```
/// use std::str::FromStr;
/// use liblrge::Platform;
///
/// for platform in ["pacbio", "pb"] {
///     assert_eq!(Platform::from_str(platform).unwrap(), Platform::PacBio);
/// }
///
/// for platform in ["nanopore", "ont"] {
///     assert_eq!(Platform::from_str(platform).unwrap(), Platform::Nanopore);
/// }
/// ```
#[derive(Debug, Default, Eq, PartialEq)]
pub enum Platform {
    PacBio,
    #[default]
    Nanopore,
}

impl FromStr for Platform {
    type Err = error::LrgeError;

    fn from_str(s: &str) -> Result<Self> {
        match s.to_lowercase().as_str() {
            "pacbio" | "pb" => Ok(Platform::PacBio),
            "nanopore" | "ont" => Ok(Platform::Nanopore),
            _ => Err(error::LrgeError::InvalidPlatform(s.to_string())),
        }
    }
}

/// Generate a shuffled list of `k` indices from 0 to `n`.
///
/// # Arguments
///
/// * `k`: The number of indices to generate.
/// * `n`: The maximum value for the range (exclusive).
/// * `seed`: An optional seed for the random number generator.
pub(crate) fn unique_random_set(k: usize, n: u32, seed: Option<u64>) -> Vec<u32> {
    // Initialize RNG, using the seed if provided
    let mut rng = match seed {
        Some(seed_value) => StdRng::seed_from_u64(seed_value),
        None => StdRng::from_entropy(),
    };

    if k > n as usize {
        panic!(
            "Cannot generate {} unique values from a range of 0 to {}",
            k, n
        );
    }

    rand::seq::index::sample(&mut rng, n as usize, k)
        .into_iter()
        .map(|x| x as u32)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn test_unique_random_set_basic_functionality() {
        let k = 5;
        let n = 100;

        for _ in 0..1000 {
            let result = unique_random_set(k, n, None);

            // Check that result has exactly k elements
            assert_eq!(result.len(), k);

            // Check that all elements are within the range 0 to n-1
            assert!(result.iter().all(|&x| x < n));

            // check that all elements are unique
            assert_eq!(result.len(), result.iter().collect::<HashSet<_>>().len());
        }
    }

    #[test]
    fn test_unique_random_set_with_seed() {
        let k = 5;
        let n = 1000000;
        let seed = Some(42);

        // Generate two sets with the same seed
        let result1 = unique_random_set(k, n, seed);
        let result2 = unique_random_set(k, n, seed);

        // They should be the same due to the same seed
        assert_eq!(result1, result2);
    }

    #[test]
    fn test_unique_random_set_without_seed() {
        let k = 5;
        let n = 10000000;

        // Generate two sets without a seed
        let result1 = unique_random_set(k, n, None);
        let result2 = unique_random_set(k, n, None);

        // They should generally be different
        assert_ne!(result1, result2);
    }

    #[test]
    #[should_panic(expected = "Cannot generate")]
    fn test_unique_random_set_k_greater_than_n() {
        let k = 10;
        let n = 5;

        // This should panic as k > n is impossible for unique values
        unique_random_set(k, n, None);
    }
}
