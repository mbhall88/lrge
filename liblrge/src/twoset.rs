//! A strategy that compares overlaps between two different sets of reads.
mod builder;

use std::collections::HashSet;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use crossbeam_channel as channel;
use log::{debug, trace, warn};
use needletail::{parse_fastx_file, parse_fastx_reader};
use rayon::prelude::*;

pub use self::builder::Builder;
use crate::estimate::per_read_estimate;
use crate::io::FastqRecordExt;
use crate::minimap2::{AlignerWrapper, Preset};
use crate::{error::LrgeError, io, unique_random_set, Estimate, Platform};

pub const DEFAULT_TARGET_NUM_READS: usize = 10_000;
pub const DEFAULT_QUERY_NUM_READS: usize = 5_000;

/// A strategy that compares overlaps between two sets of reads.
///
/// The convention is to use a smaller set of query reads and a larger set of target reads. The
/// query reads are overlapped with the target reads and an estimated genome size is calculated
/// for **each query read** based on the number of overlaps it has with the target set.
pub struct TwoSetStrategy {
    /// Path to the FASTQ file.
    input: PathBuf,
    /// The number of target reads to use in the strategy.
    target_num_reads: usize,
    /// The number of query reads to use in the strategy.
    query_num_reads: usize,
    /// The directory to which all intermediate files will be written.
    tmpdir: PathBuf,
    /// Number of threads to use with minimap2.
    threads: usize,
    /// The (optional) seed to use for randomly selecting reads.
    seed: Option<u64>,
    /// Sequencing platform of the reads.
    platform: Platform,
}

impl TwoSetStrategy {
    /// Create a new `TwoSetStrategy` with the default settings, using the given input file.
    ///
    /// To customise the strategy, use the [`Builder`] interface.
    pub fn new<P: AsRef<Path>>(input: P) -> Self {
        let builder = Builder::default();

        builder.build(input)
    }

    fn split_fastq(&mut self) -> crate::Result<(PathBuf, PathBuf, f32)> {
        debug!("Counting records in FASTQ file...");
        let n_fq_reads = {
            let mut reader = io::open_file(&self.input)?;
            io::count_fastq_records(&mut reader)?
        };
        debug!("Found {} reads in FASTQ file", n_fq_reads);

        if n_fq_reads > u32::MAX as usize {
            let msg = format!(
                "Number of reads in FASTQ file ({}) exceeds maximum allowed value ({})",
                n_fq_reads,
                u32::MAX
            );
            return Err(LrgeError::TooManyReadsError(msg));
        }

        let mut n_req_reads = self.target_num_reads + self.query_num_reads;

        if n_fq_reads <= self.query_num_reads {
            let msg = format!(
                "Number of reads in FASTQ file ({}) is <= query number of reads ({})",
                n_fq_reads, self.query_num_reads
            );
            return Err(LrgeError::TooFewReadsError(msg));
        } else if n_fq_reads < n_req_reads {
            warn!(
                "Number of reads in FASTQ file ({}) is less than the sum of target and query reads ({})",
                n_fq_reads, n_req_reads
            );
            self.target_num_reads = n_fq_reads - self.query_num_reads;
            n_req_reads = n_fq_reads;
            warn!("Using {} target reads", self.target_num_reads);
        }

        let indices = unique_random_set(n_req_reads, n_fq_reads as u32, self.seed);
        let (mut target_indices, mut query_indices) =
            split_into_hashsets(indices, self.target_num_reads);

        let target_file = self.tmpdir.join("target.fq");
        let query_file = self.tmpdir.join("query.fq");

        let reader = io::open_file(&self.input)?;
        let mut fastx_reader = parse_fastx_reader(reader).map_err(|e| {
            LrgeError::FastqParseError(format!("Error parsing input FASTQ file: {}", e))
        })?;

        debug!("Writing target and query reads to temporary files...");
        let mut target_writer = File::create(&target_file).map(BufWriter::new)?;
        let mut query_writer = File::create(&query_file).map(BufWriter::new)?;
        let mut sum_target_len = 0;
        let mut idx: u32 = 0;
        while let Some(r) = fastx_reader.next() {
            // we can unwrap here because we know the file is valid from when we counted the records
            let record = r.unwrap();

            if target_indices.remove(&idx) {
                record.write(&mut target_writer, None).map_err(|e| {
                    LrgeError::IoError(std::io::Error::new(std::io::ErrorKind::Other, e))
                })?;
                sum_target_len += record.num_bases();
            } else if query_indices.remove(&idx) {
                record.write(&mut query_writer, None).map_err(|e| {
                    LrgeError::IoError(std::io::Error::new(std::io::ErrorKind::Other, e))
                })?;
            }

            if target_indices.is_empty() && query_indices.is_empty() {
                break;
            }

            idx += 1;
        }

        let avg_target_len = sum_target_len as f32 / self.target_num_reads as f32;
        debug!("Target reads written to: {}", target_file.display());
        debug!("Query reads written to: {}", query_file.display());
        debug!("Average target read length: {}", avg_target_len);

        Ok((target_file, query_file, avg_target_len))
    }

    /// Align the query reads to the target reads and write the overlaps to a PAF file.
    fn align_reads(
        &self,
        aln_wrapper: AlignerWrapper,
        query_file: PathBuf,
        avg_target_len: f32,
    ) -> Result<Vec<f32>, LrgeError> {
        // Bounded channel to control memory usage - i.e., 10000 records in the channel at a time
        let (sender, receiver) = channel::bounded(10000);
        let aligner = Arc::clone(&aln_wrapper.aligner); // Shared reference for the producer thread
        let overlap_threshold = aln_wrapper.aligner.mapopt.min_chain_score as u32;

        // Producer: Read FASTQ records and send them to the channel
        let producer = std::thread::spawn(move || -> Result<(), LrgeError> {
            let mut fastx_reader = parse_fastx_file(query_file).map_err(|e| {
                LrgeError::FastqParseError(format!("Error parsing query FASTQ file: {}", e))
            })?;

            while let Some(record) = fastx_reader.next() {
                match record {
                    Ok(rec) => {
                        let msg =
                            io::Message::Data((rec.read_id().to_owned(), rec.seq().into_owned()));
                        if sender.send(msg).is_err() {
                            break; // Exit if the receiver is dropped
                        }
                    }
                    Err(e) => {
                        return Err(LrgeError::FastqParseError(format!(
                            "Error parsing query FASTQ file: {}",
                            e
                        )));
                    }
                }
            }

            // Close the channel to signal that no more records will be sent
            drop(sender);
            Ok(())
        });

        // Open the output PAF file for writing
        let paf_path = self.tmpdir.join("overlaps.paf");
        let mut buf = File::create(&paf_path).map(BufWriter::new)?;
        let writer = csv::WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_writer(&mut buf);
        let writer = Arc::new(Mutex::new(writer)); // thread-safe writer

        // set the number of threads to use with rayon in the following mapping code
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .map_err(|e| {
                LrgeError::ThreadError(format!("Error setting number of threads: {}", e))
            })?;

        let estimates = Vec::with_capacity(self.query_num_reads);
        let estimates = Arc::new(Mutex::new(estimates));

        debug!("Aligning reads and writing overlaps to PAF file...");
        // Consumer: Process records from the channel in parallel
        pool.install(|| -> Result<(), LrgeError> {
            receiver
                .into_iter()
                .par_bridge() // Parallelize the processing
                .try_for_each(|record| -> Result<(), LrgeError> {
                    let io::Message::Data((rid, seq)) = record;
                    trace!("Processing read: {:?}", String::from_utf8_lossy(&rid));

                    let mut qname = rid.to_owned();
                    if qname.last() != Some(&0) {
                        // Ensure the qname is null-terminated
                        qname.push(0);
                    }

                    // Use the shared aligner to perform alignment
                    let mappings = aligner.map(&seq, Some(&qname)).map_err(|e| {
                        LrgeError::MapError(format!(
                            "Error mapping read {}: {}",
                            String::from_utf8_lossy(&rid),
                            e
                        ))
                    })?;

                    let mut unique_overlaps = HashSet::new();

                    if !mappings.is_empty() {
                        {
                            let mut writer_lock = writer.lock().unwrap();
                            for mapping in &mappings {
                                // write the PafRecord to the PAF file
                                writer_lock.serialize(mapping)?;
                                unique_overlaps.insert(mapping.target_name.clone());
                            }
                        }
                    } else {
                        debug!(
                            "No mappings found for read: {:?}",
                            String::from_utf8_lossy(&rid)
                        );
                    }

                    let est = per_read_estimate(
                        seq.len(),
                        avg_target_len,
                        self.target_num_reads,
                        unique_overlaps.len(),
                        overlap_threshold,
                    );

                    trace!("Estimate for {}: {}", String::from_utf8_lossy(&rid), est);

                    {
                        // Lock the estimates vector and push the estimate
                        let mut estimates_lock = estimates.lock().unwrap();
                        estimates_lock.push(est);
                    }

                    Ok(())
                })?;
            Ok(())
        })?;

        // Wait for the producer to finish
        producer.join().map_err(|e| {
            LrgeError::ThreadError(format!("Thread panicked when joining: {:?}", e))
        })??;

        debug!("Overlaps written to: {:?}", paf_path);

        // we extract the estimates from the Arc and Mutex
        let estimates = Arc::try_unwrap(estimates)
            .map_err(|_| {
                LrgeError::ThreadError(
                    "Error unwrapping estimates Arc<Mutex<Vec<f32>>>".to_string(),
                )
            })?
            .into_inner()
            .map_err(|_| {
                LrgeError::ThreadError("Error unwrapping estimates Mutex<Vec<f32>>".to_string())
            })?;

        Ok(estimates)
    }
}

impl Estimate for TwoSetStrategy {
    fn generate_estimates(&mut self) -> crate::Result<Vec<f32>> {
        let (target_file, query_file, avg_target_len) = self.split_fastq()?;

        let preset = match self.platform {
            Platform::PacBio => Preset::AvaPb,
            Platform::Nanopore => Preset::AvaOnt,
        };

        let aligner = AlignerWrapper::new(&target_file, self.threads, preset, true)?;
        let estimates = self.align_reads(aligner, query_file, avg_target_len)?;

        Ok(estimates)
    }
}

/// Splits a `Vec` into two separate sets with potentially different sizes.
///
/// This function consumes the original `Vec` and divides its elements into
/// two new sets, `set1` and `set2`. The size of `set1` is specified by `size_first`,
/// while `set2` will contain the remaining elements. If `size_first` is larger than
/// the number of elements in `original`, all elements are placed in `set1`, and `set2`
/// will be empty.
///
/// # Arguments
///
/// * `original` - The `Vec` to be split. This set will be consumed by the function,
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
pub(crate) fn split_into_hashsets<T: std::hash::Hash + Eq>(
    mut original: Vec<T>,
    size_first: usize,
) -> (HashSet<T>, HashSet<T>) {
    let mut first_set = HashSet::with_capacity(size_first);
    let mut second_set = HashSet::with_capacity(original.len().saturating_sub(size_first));

    // Fill the first set
    for _ in 0..size_first.min(original.len()) {
        if let Some(element) = original.pop() {
            first_set.insert(element);
        }
    }

    // Fill the second set with the remaining elements
    while let Some(element) = original.pop() {
        second_set.insert(element);
    }

    (first_set, second_set)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn test_basic_split() {
        let original = vec![1, 2, 3, 4, 5];

        let (set1, set2) = split_into_hashsets(original, 3);

        assert_eq!(set1.len(), 3);
        assert_eq!(set2.len(), 2);
    }

    #[test]
    fn test_all_elements_in_set1() {
        let original = vec![1, 2, 3];

        let (set1, set2) = split_into_hashsets(original, 5);

        assert_eq!(set1.len(), 3);
        assert_eq!(set2.len(), 0);
    }

    #[test]
    fn test_all_elements_in_set2() {
        let original = vec![1, 2, 3];

        let (set1, set2) = split_into_hashsets(original, 0);

        assert_eq!(set1.len(), 0);
        assert_eq!(set2.len(), 3);
    }

    #[test]
    fn test_no_elements_lost() {
        let original = vec![1, 2, 3, 4];

        let (set1, set2) = split_into_hashsets(original.clone(), 2);

        // Verify no elements were lost
        let combined: HashSet<_> = set1.union(&set2).collect();
        assert_eq!(combined.len(), original.len());
        for elem in &original {
            assert!(combined.contains(elem));
        }
    }
}
