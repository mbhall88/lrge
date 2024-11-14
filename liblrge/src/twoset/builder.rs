use std::path::Path;
use std::path::PathBuf;

use super::TwoSetStrategy;

/// A builder for [`TwoSetStrategy`].
pub struct Builder {
    tmpdir: PathBuf,
    seed: Option<u64>,
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            tmpdir: env!("TMPDIR").into(),
            seed: None,
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

    /// Set the temporary directory for the strategy. By default, this is the `TMPDIR` environment
    /// variable.
    ///
    /// The directory must exist and be writable. In addition, it will not be cleaned up after the
    /// strategy is run.
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

    /// Build the [`TwoSetStrategy`], using the given number of target and query reads.
    ///
    /// # Examples
    ///
    /// ```
    /// use liblrge::twoset::Builder;
    ///
    /// let strategy = Builder::new().seed(42).build(1000, 100);
    /// ```
    pub fn build(self, target_num_reads: usize, query_num_reads: usize) -> TwoSetStrategy {
        TwoSetStrategy {
            target_num_reads,
            query_num_reads,
            tmpdir: self.tmpdir,
            seed: self.seed,
        }
    }
}
