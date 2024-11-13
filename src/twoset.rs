//! Code for the two-set strategy
use std::path::Path;
use anyhow::Result;
pub fn twoset_strategy<P: AsRef<Path>>(
    _input: P,
    _target_num_reads: usize,
    _query_num_reads: usize,
) -> Result<Vec<(&[u8], f32)>> {
    todo!()
}