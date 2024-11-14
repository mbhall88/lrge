//! A strategy that compares overlaps between two sets of reads.
mod builder;

pub use self::builder::Builder;
use crate::Estimate;
use std::path::PathBuf;

/// A strategy that compares overlaps between two sets of reads.
///
/// The convention is to use a smaller set of target reads and a larger set of query reads. The
/// query reads are overlapped with the target reads and an estimated genome size is calculated
/// for **each target read** based on the number of overlaps it has with the query set.
pub struct TwoSetStrategy {
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
    /// Create a new `TwoSetStrategy` with the default settings and the given number of target and
    /// query reads.
    ///
    /// To customize the strategy, use the [`Builder`] interface.
    pub fn new(target_num_reads: usize, query_num_reads: usize) -> Self {
        let builder = Builder::default();

        builder.build(target_num_reads, query_num_reads)
    }
}

impl Estimate for TwoSetStrategy {
    fn generate_estimates(&self) -> Vec<(&[u8], f32)> {
        todo!()
    }
}
