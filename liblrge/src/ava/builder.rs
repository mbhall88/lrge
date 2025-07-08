use std::path::{Path, PathBuf};

use super::{AvaStrategy, DEFAULT_AVA_NUM_READS};
use crate::Platform;

/// A builder for [`AvaStrategy`].
pub struct Builder {
    num_reads: usize,
    num_bases: usize,
    remove_internal: bool,
    max_overhang_size: i32,
    max_overhang_ratio: f32,
    tmpdir: PathBuf,
    threads: usize,
    seed: Option<u64>,
    platform: Platform,
}

impl Default for Builder {
    fn default() -> Self {
        let tmpdir = std::env::temp_dir();
        Self {
            num_reads: DEFAULT_AVA_NUM_READS,
            num_bases: 0,
            remove_internal: false,
            max_overhang_size: 1000,
            max_overhang_ratio: 0.8,
            tmpdir,
            threads: 1,
            seed: None,
            platform: Platform::default(),
        }
    }
}

impl Builder {
    /// Create a new builder with the default settings.
    ///
    /// # Examples
    ///
    /// ```
    /// use liblrge::ava::Builder;
    ///
    /// let builder = Builder::new();
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the number of reads for the strategy. By default, this is [`DEFAULT_AVA_NUM_READS`].
    ///
    /// # Examples
    ///
    /// ```
    /// use liblrge::ava::Builder;
    ///
    /// let builder = Builder::new().num_reads(1000);
    /// ```
    pub fn num_reads(mut self, num_reads: usize) -> Self {
        self.num_reads = num_reads;
        self
    }

    /// Set option for removing the overlaps representing internal matches
    pub fn remove_internal(mut self, do_filt: bool, size: i32, ratio: f32) -> Self {
        self.remove_internal = do_filt;
        if do_filt {
            self.max_overhang_size = size;
            self.max_overhang_ratio = ratio;
        }
        self
    }

    /// Set the temporary directory for the strategy. By default, this is the value of the `TMPDIR`
    /// environment variable.
    ///
    /// # Examples
    ///
    /// ```
    /// use liblrge::ava::Builder;
    /// use std::path::PathBuf;
    ///
    /// let builder = Builder::new().tmpdir(PathBuf::from("/tmp"));
    /// ```
    pub fn tmpdir<P: AsRef<Path>>(mut self, tmpdir: P) -> Self {
        self.tmpdir = tmpdir.as_ref().to_path_buf();
        self
    }

    /// Set the number of threads to use with minimap2. By default, this is `1`.
    ///
    /// # Examples
    ///
    /// ```
    /// use liblrge::ava::Builder;
    ///
    /// let builder = Builder::new().threads(4);
    /// ```
    pub fn threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    /// Set the seed for the strategy. By default (`None`), the seed will be
    /// [randomly generated](https://docs.rs/rand/latest/rand/fn.random.html).
    ///
    /// # Examples
    ///
    /// ```
    /// use liblrge::ava::Builder;
    ///
    /// let builder = Builder::new().seed(Some(42));
    /// ```
    pub fn seed(mut self, seed: Option<u64>) -> Self {
        self.seed = seed;
        self
    }

    /// Set the sequencing platform for the reads. By default, this is [`Platform::default()`].
    ///
    /// # Examples
    ///
    /// ```
    /// use liblrge::{ava::Builder, Platform};
    ///
    /// let builder = Builder::new().platform(Platform::PacBio);
    /// ```
    pub fn platform(mut self, platform: Platform) -> Self {
        self.platform = platform;
        self
    }

    /// Build the [`AvaStrategy`], using the reads from the given `input` file.
    ///
    /// # Examples
    ///
    /// ```
    /// use liblrge::ava::Builder;
    ///
    /// let builder = Builder::new().num_reads(1000);
    /// let strategy = builder.build("reads.fq");
    /// ```
    pub fn build<P: AsRef<Path>>(self, input: P) -> AvaStrategy {
        AvaStrategy {
            input: input.as_ref().to_path_buf(),
            num_reads: self.num_reads,
            num_bases: self.num_bases,
            remove_internal: self.remove_internal,
            max_overhang_size: self.max_overhang_size,
            max_overhang_ratio: self.max_overhang_ratio,
            tmpdir: self.tmpdir,
            threads: self.threads,
            seed: self.seed,
            platform: self.platform,
        }
    }
}
