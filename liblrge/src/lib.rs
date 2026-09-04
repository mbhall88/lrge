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
//! let input = "path/to/reads.fastq"; // or .fasta, .bam, .cram, .sam
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
//! let input = "path/to/reads.fastq"; // or .fasta, .bam, .cram, .sam
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
//! This library includes optional support for compressed file formats and alignment formats, controlled by feature flags.
//! By default, the `compression` and `alignment` features are enabled.
//!
//! ### Available Features
//!
//! - **compression** (default): Enables all available compression formats (`gzip`, `zstd`, `bzip2`, `xz`).
//! - **alignment** (default): Enables support for unaligned BAM, CRAM, and SAM formats using the [`noodles`][noodles] crate.
//! - **gzip**: Enables support for gzip-compressed files (`.gz`) using the [`flate2`][flate2] crate.
//! - **zstd**: Enables support for zstd-compressed files (`.zst`) using the [`zstd`][zstd] crate.
//! - **bzip2**: Enables support for bzip2-compressed files (`.bz2`) using the [`bzip2`][bzip2] crate.
//! - **xz**: Enables support for xz-compressed files (`.xz`) using the [`liblzma`][xz] crate.
//!
//! ### Enabling and Disabling Features
//!
//! By default, all features are enabled. However, you can selectively enable or disable them
//! in your `Cargo.toml` to reduce dependencies or target specific formats:
//!
//! To **disable all optional features**:
//!
//! ```toml
//! liblrge = { version = "0.2.2", default-features = false }
//! ```
//!
//! To enable only specific features, list them in `Cargo.toml`:
//!
//! ```toml
//! liblrge = { version = "0.2.2", default-features = false, features = ["gzip", "alignment"] }
//! ```
//!
//! ## Format Detection
//!
//! The library uses [**magic bytes**][magic] at the start of the file to detect its compression
//! format and content type before deciding how to read it. Supported formats include:
//! - **FASTX**: FASTA and FASTQ (via `needletail`).
//! - **Alignment**: BAM, CRAM, and SAM (via `noodles`). Alignment files must be **unaligned**.
//! - **Compression**: gzip, zstd, bzip2, and xz (automatic decompression if the [appropriate feature](#features) is enabled).
//!
//! [flate2]: https://crates.io/crates/flate2
//! [zstd]: https://crates.io/crates/zstd
//! [xz]: https://crates.io/liblzma
//! [bzip2]: https://crates.io/crates/bzip2
//! [noodles]: https://crates.io/crates/noodles
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
//! [doi]: https://doi.org/10.1101/2024.11.27.625777
#[deny(missing_docs)]
pub mod ava;
pub(crate) mod depth_skew;
pub mod error;
pub mod estimate;
pub(crate) mod io;
pub(crate) mod minimap2;
pub(crate) mod read_selection;
pub mod twoset;

/// The default maximum ratio of overhang to alignment length above which an alignment is
/// treated as an internal match. Shared by both strategies.
pub const DEFAULT_MAX_OVERHANG_RATIO: f32 = 0.2;

/// The default cap on the bytes of selected reads that depth normalization holds in memory.
///
/// Normalization buffers the reads it selects so it can write them once sampling is done. A
/// request large enough to project past this cap is served instead by a low-memory path that
/// buffers read positions and re-reads the input, trading one extra pass for a bounded buffer.
/// The two paths select the same reads for a given seed, so the cap governs how much memory a
/// run takes while the estimate stays the same.
///
/// The cap covers the read buffer alone: the depth profile, the minimap2 index, and the overlap
/// stage are all outside it. It is also applied to a projection made from the mean read length,
/// so an input whose long reads survive normalization in unusual numbers can still buffer past
/// it. Such a run says by how much once it is done.
pub const DEFAULT_MAX_READ_BUFFER: u64 = 1 << 30;

use rand::rngs::StdRng;
use rand::SeedableRng;

pub use self::ava::AvaStrategy;
pub use self::estimate::Estimate;
pub use self::twoset::TwoSetStrategy;
use std::str::FromStr;

/// A type alias for `Result` with [`LrgeError`][crate::error::LrgeError] as the error type.
pub type Result<T> = std::result::Result<T, error::LrgeError>;

/// Controls whether read-depth normalization is applied before estimation.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum Normalization {
    /// Detect depth skew and normalize only when skew is present.
    #[default]
    Auto,
    /// Normalize every input without consulting the skew verdict.
    ///
    /// Detection still runs, because what normalization measures depth against is read off the
    /// minimizers detection samples. The verdict is what this mode ignores, not the work.
    Always,
    /// Disable depth-skew detection and normalization.
    Never,
}

impl FromStr for Normalization {
    type Err = String;

    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "auto" => Ok(Self::Auto),
            "always" => Ok(Self::Always),
            "never" => Ok(Self::Never),
            _ => Err(format!(
                "invalid normalization mode '{value}'; expected auto, always, or never"
            )),
        }
    }
}

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

/// Sample `k` distinct indices from 0 to `n`, optionally in proportion to weights.
///
/// # Arguments
///
/// * `k`: The number of indices to generate.
/// * `n`: The maximum value for the range (exclusive).
/// * `seed`: An optional seed for the random number generator.
/// * `weights`: An optional weight for each index.
pub(crate) fn sample_unique_indices(
    k: usize,
    n: u32,
    seed: Option<u64>,
    weights: Option<&[f64]>,
) -> Vec<u32> {
    // Initialize RNG, using the seed if provided
    let mut rng = match seed {
        Some(seed_value) => StdRng::seed_from_u64(seed_value),
        None => StdRng::from_rng(&mut rand::rng()),
    };

    if k > n as usize {
        panic!("Cannot generate {k} unique values from a range of 0 to {n}",);
    }

    if let Some(weights) = weights {
        assert_eq!(weights.len(), n as usize, "expected one weight per read");
        assert!(
            weights
                .iter()
                .all(|weight| weight.is_finite() && *weight >= 0.0),
            "read weights must be finite and non-negative"
        );

        let uniform = weights
            .first()
            .is_none_or(|first| *first > 0.0 && weights.iter().all(|weight| weight == first));
        if !uniform {
            let positive_weights = weights.iter().filter(|weight| **weight > 0.0).count();
            if k > positive_weights {
                panic!("Cannot generate {k} unique values from a range of 0 to {positive_weights}");
            }

            return rand::seq::index::sample_weighted(
                &mut rng,
                n as usize,
                |index| weights[index],
                k,
            )
            .expect("read weights were validated")
            .into_iter()
            .map(|index| index as u32)
            .collect();
        }
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
    fn test_sample_unique_indices_basic_functionality() {
        let k = 5;
        let n = 100;

        for _ in 0..1000 {
            let result = sample_unique_indices(k, n, None, None);

            // Check that result has exactly k elements
            assert_eq!(result.len(), k);

            // Check that all elements are within the range 0 to n-1
            assert!(result.iter().all(|&x| x < n));

            // check that all elements are unique
            assert_eq!(result.len(), result.iter().collect::<HashSet<_>>().len());
        }
    }

    #[test]
    fn test_sample_unique_indices_with_seed() {
        let k = 5;
        let n = 1000000;
        let seed = Some(42);

        // Generate two sets with the same seed
        let result1 = sample_unique_indices(k, n, seed, None);
        let result2 = sample_unique_indices(k, n, seed, None);

        // They should be the same due to the same seed
        assert_eq!(result1, result2);
    }

    #[test]
    fn test_sample_unique_indices_without_seed() {
        let k = 5;
        let n = 10000000;

        // Generate two sets without a seed
        let result1 = sample_unique_indices(k, n, None, None);
        let result2 = sample_unique_indices(k, n, None, None);

        // They should generally be different
        assert_ne!(result1, result2);
    }

    #[test]
    #[should_panic(expected = "Cannot generate")]
    fn test_sample_unique_indices_k_greater_than_n() {
        let k = 10;
        let n = 5;

        // This should panic as k > n is impossible for unique values
        sample_unique_indices(k, n, None, None);
    }

    #[test]
    fn test_weighted_random_set_draws_in_proportion_to_weights() {
        let weights = [1.0, 3.0];
        let mut selections = [0_u32; 2];

        for seed in 0..4096 {
            let selected = sample_unique_indices(1, 2, Some(seed), Some(&weights));
            selections[selected[0] as usize] += 1;
        }

        let ratio = selections[1] as f64 / selections[0] as f64;
        assert!((2.5..3.5).contains(&ratio), "selection ratio was {ratio}");
    }

    #[test]
    fn test_weighted_random_set_with_seed() {
        let weights = [1.0, 2.0, 3.0, 4.0, 5.0];
        let seed = Some(42);

        let result1 = sample_unique_indices(2, 5, seed, Some(&weights));
        let result2 = sample_unique_indices(2, 5, seed, Some(&weights));

        assert_eq!(result1, result2);
    }

    #[test]
    #[should_panic(expected = "Cannot generate 4 unique values from a range of 0 to 3")]
    fn test_weighted_random_set_k_greater_than_n() {
        sample_unique_indices(4, 3, Some(42), Some(&[1.0, 2.0, 3.0]));
    }

    #[test]
    #[should_panic(expected = "Cannot generate")]
    fn test_weighted_random_set_k_greater_than_positive_weights() {
        sample_unique_indices(2, 3, Some(42), Some(&[1.0, 0.0, 0.0]));
    }
}
