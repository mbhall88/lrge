//! # liblrge
//!
//! `liblrge` is a Rust library that provides utilities for estimating genome size for a given set
//! of reads.
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
// #![deny(missing_docs)]
pub mod estimate;
pub(crate) mod io;
pub mod twoset;

pub use self::estimate::Estimate;
pub use self::twoset::TwoSetStrategy;
use log::debug;
use rand::prelude::SliceRandom;
use rand::{random, SeedableRng};

/// Returns a vector of indices for the number of elements `n`, but shuffled.
///
/// If a `seed` is provided, the shuffle will be deterministic.
pub(crate) fn shuffled_indices(n: u32, seed: Option<u64>) -> Vec<u32> {
    let mut indices: Vec<u32> = (0..n).collect();

    let mut rng = match seed {
        Some(s) => rand_pcg::Pcg64::seed_from_u64(s),
        None => {
            let seed = random();
            debug!("Using seed: {}", seed);
            rand_pcg::Pcg64::seed_from_u64(seed)
        }
    };

    indices.shuffle(&mut rng);
    indices
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shuffled_indices() {
        let n = 3;
        let mut num_times_shuffled = 0;
        let iterations = 100;
        for _ in 0..iterations {
            let idxs = shuffled_indices(n, None);
            if idxs != vec![0, 1, 2] {
                num_times_shuffled += 1;
            }
        }
        // chances of shuffling the same way - i.e., [0, 1, 2] - 100 times in a row is 3.054936363499605e-151
        assert!(num_times_shuffled > 0 && num_times_shuffled < iterations)
    }

    #[test]
    fn test_shuffled_indices_with_seed() {
        let n = 4;
        let seed = 42;
        let indices = shuffled_indices(n, Some(seed));

        for _ in 0..100 {
            let new_indices = shuffled_indices(n, Some(seed));
            assert_eq!(indices, new_indices);
        }
    }
}
