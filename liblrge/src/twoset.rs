//! A strategy that compares overlaps between two sets of reads.
mod builder;

pub use self::builder::Builder;
use crate::{io, Estimate, shuffled_indices};
use log::debug;
use std::path::{Path, PathBuf};

pub const DEFAULT_TARGET_NUM_READS: usize = 5000;
pub const DEFAULT_QUERY_NUM_READS: usize = 10000;

/// A strategy that compares overlaps between two sets of reads.
///
/// The convention is to use a smaller set of target reads and a larger set of query reads. The
/// query reads are overlapped with the target reads and an estimated genome size is calculated
/// for **each target read** based on the number of overlaps it has with the query set.
pub struct TwoSetStrategy {
    /// Path to the FASTQ file.
    input: PathBuf,
    /// The number of target reads to use in the strategy.
    target_num_reads: usize,
    /// The number of query reads to use in the strategy.
    query_num_reads: usize,
    /// The directory to which all intermediate files will be written.
    tmpdir: PathBuf,
    /// The (optional) seed to use for randomly selecting reads.
    seed: Option<u64>,
}

impl TwoSetStrategy {
    /// Create a new `TwoSetStrategy` with the default settings, using the given input file.
    ///
    /// To customize the strategy, use the [`Builder`] interface.
    pub fn new<P: AsRef<Path>>(input: P) -> Self {
        let builder = Builder::default();

        builder.build(input)
    }

    fn split_fastq(&self) -> std::io::Result<(PathBuf, PathBuf)> {
        debug!("Counting records in FASTQ file...");
        let num_reads = {
            let mut reader = io::open_file(&self.input)?;
            io::count_fastq_records(&mut reader)?
        };
        debug!("Found {} reads in FASTQ file", num_reads);

        if num_reads > u32::MAX as usize {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "Number of reads exceeds maximum supported value of {}",
                    u32::MAX
                ),
            ));
        }
        let mut indices = shuffled_indices(num_reads as u32, self.seed);
        todo!("Split the indices into target and query sizes")
    }
}

impl Estimate for TwoSetStrategy {
    fn generate_estimates(&self) -> Vec<f32> {
        todo!()
    }
}
