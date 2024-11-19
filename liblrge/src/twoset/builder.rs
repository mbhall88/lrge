use crate::Platform;
use std::path::Path;
use std::path::PathBuf;

use super::{TwoSetStrategy, DEFAULT_QUERY_NUM_READS, DEFAULT_TARGET_NUM_READS};

/// A builder for [`TwoSetStrategy`].
pub struct Builder {
    target_num_reads: usize,
    query_num_reads: usize,
    tmpdir: PathBuf,
    threads: usize,
    seed: Option<u64>,
    platform: Platform,
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            target_num_reads: DEFAULT_TARGET_NUM_READS,
            query_num_reads: DEFAULT_QUERY_NUM_READS,
            tmpdir: env!("TMPDIR").into(),
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
    /// use liblrge::twoset::Builder;
    ///
    /// let builder = Builder::new();
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the number of target reads for the strategy. By default, this is [`DEFAULT_TARGET_NUM_READS`].
    ///
    /// The target reads are the (generally) smaller set of reads that the query reads are
    /// overlapped against.
    ///
    /// # Examples
    ///
    /// ```
    /// use liblrge::twoset::Builder;
    ///
    /// let builder = Builder::new().target_num_reads(1000);
    /// ```
    pub fn target_num_reads(mut self, target_num_reads: usize) -> Self {
        self.target_num_reads = target_num_reads;
        self
    }

    /// Set the number of query reads for the strategy. By default, this is [`DEFAULT_QUERY_NUM_READS`].
    ///
    /// The query reads are the (generally) larger set of reads that are overlapped against the
    /// target reads.
    ///
    /// # Examples
    ///
    /// ```
    /// use liblrge::twoset::Builder;
    ///
    /// let builder = Builder::new().query_num_reads(1000);
    /// ```
    pub fn query_num_reads(mut self, query_num_reads: usize) -> Self {
        self.query_num_reads = query_num_reads;
        self
    }

    /// Set the number of threads to use with minimap2. By default, this is 1.
    pub fn threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    /// Set the temporary directory for the strategy. By default, this is the `TMPDIR` environment
    /// variable.
    ///
    /// The directory will not be created, nor will it be cleaned up after the strategy is run.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::path::Path;
    /// use liblrge::twoset::Builder;
    ///
    /// let builder = Builder::new().tmpdir(Path::new("/my-temp-dir"));
    /// ```
    ///
    /// If you want the temporary directory to clean up after the strategy is run, you can use the
    /// [`tempfile`](https://crates.io/crates/tempfile) crate to create a temporary directory.
    pub fn tmpdir<P: AsRef<Path>>(mut self, tmpdir: P) -> Self {
        self.tmpdir = tmpdir.as_ref().to_path_buf();
        self
    }

    /// Set the seed for the strategy. By default, the seed will be
    /// [randomly generated](https://docs.rs/rand/latest/rand/fn.random.html).
    ///
    /// # Examples
    ///
    /// ```
    /// use liblrge::twoset::Builder;
    ///
    /// let builder = Builder::new().seed(42);
    /// ```
    pub fn seed(mut self, seed: u64) -> Self {
        self.seed = Some(seed);
        self
    }

    /// Set the sequencing platform for the strategy. By default, this is [`Platform::Nanopore`].
    ///
    /// # Examples
    ///
    /// ```
    /// use liblrge::{twoset::Builder, Platform};
    ///
    /// let builder = Builder::new().platform(Platform::PacBio);
    /// ```
    pub fn platform(mut self, platform: Platform) -> Self {
        self.platform = platform;
        self
    }

    /// Build the [`TwoSetStrategy`], using the given number of target and query reads.
    ///
    /// # Examples
    ///
    /// ```
    /// use liblrge::twoset::Builder;
    ///
    /// let strategy = Builder::new().target_num_reads(1000).build("input.fastq");
    /// ```
    pub fn build<P: AsRef<Path>>(self, input: P) -> TwoSetStrategy {
        TwoSetStrategy {
            input: input.as_ref().to_path_buf(),
            target_num_reads: self.target_num_reads,
            query_num_reads: self.query_num_reads,
            tmpdir: self.tmpdir,
            threads: self.threads,
            seed: self.seed,
            platform: self.platform,
        }
    }
}
