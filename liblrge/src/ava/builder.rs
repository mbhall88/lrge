use std::path::{Path, PathBuf};

use super::{AvaStrategy, DEFAULT_AVA_NUM_READS};
use crate::{Normalization, Platform};
use crate::{DEFAULT_MAX_OVERHANG_RATIO, DEFAULT_MAX_READ_BUFFER};

/// A builder for [`AvaStrategy`].
pub struct Builder {
    num_reads: usize,
    num_bases: usize,
    internal_filter: Option<f64>,
    max_overhang_ratio: f32,
    tmpdir: PathBuf,
    threads: usize,
    seed: Option<u64>,
    platform: Platform,
    normalization: Normalization,
    max_read_buffer: u64,
}

impl Default for Builder {
    fn default() -> Self {
        let tmpdir = std::env::temp_dir();
        Self {
            num_reads: DEFAULT_AVA_NUM_READS,
            num_bases: 0,
            internal_filter: None,
            max_overhang_ratio: DEFAULT_MAX_OVERHANG_RATIO,
            tmpdir,
            threads: 1,
            seed: None,
            platform: Platform::default(),
            normalization: Normalization::default(),
            max_read_buffer: DEFAULT_MAX_READ_BUFFER,
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

    /// Set the share of a run's overlaps that has to be internal matches before they are excluded,
    /// and the maximum ratio of overhang to alignment length above which a mapping counts as one.
    ///
    /// `None` keeps every overlap, which is the default. `Some(0.0)` excludes whatever internal
    /// matches the run finds. In between, the run measures its own share during the overlap pass
    /// and excludes them only if it clears the threshold; see
    /// [`DEFAULT_INTERNAL_MATCH_THRESHOLD`][crate::DEFAULT_INTERNAL_MATCH_THRESHOLD] for what the
    /// benchmark says about choosing one.
    ///
    /// The ratio is stored either way, so it survives being set before the filter is turned on.
    pub fn internal_filter(mut self, threshold: Option<f64>, ratio: f32) -> Self {
        self.internal_filter = threshold;
        self.max_overhang_ratio = ratio;
        self
    }

    /// Set option for removing the overlaps representing internal matches, and the maximum
    /// ratio of overhang to alignment length above which a mapping counts as one.
    #[deprecated(
        since = "0.4.0",
        note = "use `internal_filter`, which takes the share of a run's overlaps that internal \
                matches have to exceed rather than a flag. `remove_internal(true, ratio)` is \
                `internal_filter(Some(0.0), ratio)`, and `remove_internal(false, ratio)` is \
                `internal_filter(None, ratio)`."
    )]
    pub fn remove_internal(self, filter_contained: bool, ratio: f32) -> Self {
        // Excluding every internal match is a threshold of zero: any run holding one at all has a
        // share above it, and a run holding none has nothing the filter could have taken.
        self.internal_filter(filter_contained.then_some(0.0), ratio)
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

    /// Set when depth-aware read normalization is applied. The default is
    /// [`Normalization::Auto`].
    pub fn normalization(mut self, normalization: Normalization) -> Self {
        self.normalization = normalization;
        self
    }

    /// Set the cap on the bytes of selected reads that depth normalization may hold in memory
    /// at once. By default, this is [`DEFAULT_MAX_READ_BUFFER`].
    ///
    /// A request projected to exceed the cap falls back to a low-memory path that pays one extra
    /// pass over the input. Both paths select the same reads for a given seed, so the cap governs
    /// how much memory a run takes while the estimate stays the same.
    ///
    /// # Examples
    ///
    /// ```
    /// use liblrge::ava::Builder;
    ///
    /// let builder = Builder::new().max_read_buffer(4 << 30);
    /// ```
    pub fn max_read_buffer(mut self, bytes: u64) -> Self {
        self.max_read_buffer = bytes;
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
            internal_filter: self.internal_filter,
            max_overhang_ratio: self.max_overhang_ratio,
            tmpdir: self.tmpdir,
            threads: self.threads,
            seed: self.seed,
            platform: self.platform,
            normalization: self.normalization,
            max_read_buffer: self.max_read_buffer,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DEFAULT_INTERNAL_MATCH_THRESHOLD;

    #[test]
    fn internal_filter_keeps_the_ratio_when_the_filter_is_off() {
        // the ratio is an independent setting - storing it only when the filter is enabled
        // silently discards a caller's choice
        let strategy = Builder::new().internal_filter(None, 0.05).build("reads.fq");

        assert_eq!(strategy.internal_filter, None);
        assert_eq!(strategy.max_overhang_ratio, 0.05);
    }

    #[test]
    fn internal_filter_keeps_the_ratio_at_every_threshold() {
        for threshold in [
            None,
            Some(0.0),
            Some(0.5),
            Some(DEFAULT_INTERNAL_MATCH_THRESHOLD),
        ] {
            let strategy = Builder::new()
                .internal_filter(threshold, 0.05)
                .build("reads.fq");

            assert_eq!(strategy.internal_filter, threshold);
            assert_eq!(strategy.max_overhang_ratio, 0.05);
        }
    }

    /// The old flag is kept so a caller has a release to migrate in, and it has to mean what it
    /// meant: excluding every internal match, which is a threshold of zero.
    #[test]
    #[allow(deprecated)]
    fn remove_internal_still_maps_onto_a_threshold() {
        let on = Builder::new().remove_internal(true, 0.05).build("reads.fq");
        let off = Builder::new()
            .remove_internal(false, 0.05)
            .build("reads.fq");

        assert_eq!(on.internal_filter, Some(0.0));
        assert_eq!(on.max_overhang_ratio, 0.05);
        assert_eq!(off.internal_filter, None);
        assert_eq!(off.max_overhang_ratio, 0.05);
    }

    #[test]
    fn internal_filter_defaults_to_off_with_the_default_ratio() {
        let strategy = Builder::new().build("reads.fq");

        assert_eq!(strategy.internal_filter, None);
        assert_eq!(strategy.max_overhang_ratio, DEFAULT_MAX_OVERHANG_RATIO);
    }
}
