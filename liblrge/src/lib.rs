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
