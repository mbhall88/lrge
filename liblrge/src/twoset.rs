//! A strategy that compares overlaps between two sets of reads.
mod builder;

pub use self::builder::Builder;
use crate::{io, unique_random_set, Estimate};
use log::{debug, warn};
use std::collections::HashSet;
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

    fn split_fastq(&mut self) -> std::io::Result<(PathBuf, PathBuf)> {
        debug!("Counting records in FASTQ file...");
        let n_fq_reads = {
            let mut reader = io::open_file(&self.input)?;
            io::count_fastq_records(&mut reader)?
        };
        debug!("Found {} reads in FASTQ file", n_fq_reads);

        if n_fq_reads > u32::MAX as usize {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "Number of reads exceeds maximum supported value of {}",
                    u32::MAX
                ),
            ));
        }

        let mut n_req_reads = self.target_num_reads + self.query_num_reads;

        if n_fq_reads <= self.target_num_reads {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "Number of reads in FASTQ file ({}) is <= target number of reads ({})",
                    n_fq_reads, self.target_num_reads
                ),
            ));
        } else if n_fq_reads < n_req_reads {
            warn!(
                "Number of reads in FASTQ file ({}) is less than the sum of target and query reads ({})",
                n_fq_reads, n_req_reads
            );
            self.query_num_reads = n_fq_reads - self.target_num_reads;
            n_req_reads = n_fq_reads;
            warn!("Using {} query reads instead", self.query_num_reads);
        }

        let indices = unique_random_set(n_req_reads, n_fq_reads as u32, self.seed);
        let (target_indices, query_indices) =
            split_hashset_into_two(indices, self.target_num_reads);
        todo!("Split the indices into target and query sizes")
    }
}

/// Splits a `HashSet` into two separate sets with potentially different sizes.
///
/// This function consumes the original `HashSet` and divides its elements into
/// two new sets, `set1` and `set2`. The size of `set1` is specified by `size_first`,
/// while `set2` will contain the remaining elements. If `size_first` is larger than
/// the number of elements in `original`, all elements are placed in `set1`, and `set2`
/// will be empty.
///
/// # Arguments
///
/// * `original` - The `HashSet` to be split. This set will be consumed by the function,
///                so it will no longer be accessible after the function call.
/// * `size_first` - The number of elements to place in the first set, `set1`.
///
/// # Returns
///
/// A tuple containing:
/// * `HashSet<T>` - The first set (`set1`), containing up to `size_first` elements.
/// * `HashSet<T>` - The second set (`set2`), containing the remaining elements.
///
/// # Panics
///
/// This function will panic if `size_first` is larger than `original.len()`.
///
pub(crate) fn split_hashset_into_two<T: std::hash::Hash + Eq>(
    mut original: HashSet<T>,
    size_first: usize,
) -> (HashSet<T>, HashSet<T>) {
    let mut set1 = HashSet::with_capacity(size_first);
    let mut set2 = HashSet::with_capacity(original.len().saturating_sub(size_first));

    // Drain items from `original`, moving items into `set1` until it reaches the desired size
    for item in original.drain() {
        if set1.len() < size_first {
            set1.insert(item);
        } else {
            set2.insert(item);
        }
    }

    (set1, set2)
}

impl Estimate for TwoSetStrategy {
    fn generate_estimates(&self) -> Vec<f32> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn test_basic_split() {
        let original: HashSet<_> = vec![1, 2, 3, 4, 5].into_iter().collect();

        let (set1, set2) = split_hashset_into_two(original, 3);

        assert_eq!(set1.len(), 3);
        assert_eq!(set2.len(), 2);
    }

    #[test]
    fn test_all_elements_in_set1() {
        let original: HashSet<_> = vec![1, 2, 3].into_iter().collect();

        let (set1, set2) = split_hashset_into_two(original, 5);

        assert_eq!(set1.len(), 3);
        assert_eq!(set2.len(), 0);
    }

    #[test]
    fn test_all_elements_in_set2() {
        let original: HashSet<_> = vec![1, 2, 3].into_iter().collect();

        let (set1, set2) = split_hashset_into_two(original, 0);

        assert_eq!(set1.len(), 0);
        assert_eq!(set2.len(), 3);
    }

    #[test]
    fn test_no_elements_lost() {
        let original: HashSet<_> = vec![1, 2, 3, 4].into_iter().collect();

        let (set1, set2) = split_hashset_into_two(original.clone(), 2);

        // Verify no elements were lost
        let combined: HashSet<_> = set1.union(&set2).collect();
        assert_eq!(combined.len(), original.len());
        for elem in &original {
            assert!(combined.contains(elem));
        }
    }
}
