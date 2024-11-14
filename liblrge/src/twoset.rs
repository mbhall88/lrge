//! Code for the two-set strategy
mod builder;

pub use self::builder::Builder;
use crate::Estimate;
use std::path::PathBuf;
pub struct TwoSetStrategy {
    target_num_reads: usize,
    query_num_reads: usize,
    tmpdir: PathBuf,
    seed: Option<u64>,
}

impl TwoSetStrategy {
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
