//! A strategy that compares overlaps between two different sets of reads.
//!
//! The convention is to use a smaller set of query reads and a larger set of target reads. The query
//! reads are overlapped with the target reads and an estimated genome size is calculated for
//! **each query read** based on the number of overlaps it has with the target set.
//!
//! This strategy is generally faster than the [`AvaStrategy`][crate::AvaStrategy] and uses less memory.
//!
//! # Examples
//!
//! You probably want to use the [`Builder`] interface to customise the strategy.
//!
//! ```no_run
//! use liblrge::{Estimate, TwoSetStrategy};
//! use liblrge::estimate::{LOWER_QUANTILE, UPPER_QUANTILE};
//! use liblrge::twoset::{Builder, DEFAULT_TARGET_NUM_READS, DEFAULT_QUERY_NUM_READS};
//!
//! let input = "path/to/reads.fastq";
//! let mut strategy = Builder::new()
//!    .target_num_reads(DEFAULT_TARGET_NUM_READS)
//!    .query_num_reads(DEFAULT_QUERY_NUM_READS)
//!    .threads(4)
//!    .seed(Some(42))  // makes the estimate reproducible
//!    .build(input);
//!
//! let finite = true;  // estimate the genome size based on the finite estimates (recommended)
//! let low_q = Some(LOWER_QUANTILE);   // lower quantile for the confidence interval
//! let upper_q = Some(UPPER_QUANTILE); // upper quantile for the confidence interval
//! let est_result = strategy.estimate(finite, low_q, upper_q).expect("Failed to generate estimate");
//! let estimate = est_result.estimate;
//!
//! let no_mapping_count = est_result.no_mapping_count;
//! // you might want to handle cases where some proportion of query reads did not overlap with target reads
//! ```
//!
//! By default, the intermediate target and query reads and overlap files are written to a temporary
//! directory and cleaned up after the strategy object is dropped. This is done via the use of the
//! [`tempfile`](https://crates.io/crates/tempfile) crate. The intermediate read files are placed in
//! the temporary directory and named `target.fq` and `query.fq`, while the overlap file is named
//! `overlaps.paf`.
//!
//! You can set your own temporary directory by using the [`Builder::tmpdir`] method.
mod builder;
use std::cmp;
use std::collections::{HashMap, HashSet};
use std::ffi::CString;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::sync::atomic::AtomicU32;
use std::sync::{Arc, Mutex};

use crossbeam_channel as channel;
use log::{debug, info, trace, warn};
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
///
/// See the [module-level documentation](crate::twoset) for more information and examples.
pub struct TwoSetStrategy {
    /// Path to the FASTQ file.
    input: PathBuf,
    /// The number of target reads to use in the strategy.
    target_num_reads: usize,
    /// The number of target bases to use in the strategy.
    target_num_bases: usize,
    /// The number of query reads to use in the strategy.
    query_num_reads: usize,
    /// The number of query bases to use in the strategy.
    query_num_bases: usize,
    /// Remove overlaps for internal matches.
    remove_internal: bool,
    /// Maximum overhang ratio
    max_overhang_ratio: f32,
    /// Use the smaller Q/T dataset as minimap2 reference
    use_min_ref: bool,
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

    /// The number of target reads
    pub fn target_num_reads(&self) -> usize {
        self.target_num_reads
    }

    /// The number of query reads
    pub fn query_num_reads(&self) -> usize {
        self.query_num_reads
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
                "Number of reads in FASTQ file ({n_fq_reads}) exceeds maximum allowed value ({})",
                u32::MAX
            );
            return Err(LrgeError::TooManyReadsError(msg));
        }

        let mut n_req_reads = self.target_num_reads + self.query_num_reads;

        if n_fq_reads <= self.query_num_reads {
            let msg = format!(
                "Number of reads in FASTQ file ({n_fq_reads}) is <= query number of reads ({})",
                self.query_num_reads
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
            LrgeError::FastqParseError(format!("Error parsing input FASTQ file: {e}",))
        })?;

        debug!("Writing target and query reads to temporary files...");
        let mut target_writer = File::create(&target_file).map(BufWriter::new)?;
        let mut query_writer = File::create(&query_file).map(BufWriter::new)?;
        let mut sum_target_len = 0;
        let mut sum_query_len: usize = 0;
        let mut idx: u32 = 0;
        while let Some(r) = fastx_reader.next() {
            // we can unwrap here because we know the file is valid from when we counted the records
            let record = r.unwrap();

            if target_indices.remove(&idx) {
                record
                    .write(&mut target_writer, None)
                    .map_err(|e| LrgeError::IoError(std::io::Error::other(e)))?;
                sum_target_len += record.num_bases();
            } else if query_indices.remove(&idx) {
                record
                    .write(&mut query_writer, None)
                    .map_err(|e| LrgeError::IoError(std::io::Error::other(e)))?;
                sum_query_len += record.num_bases();
            }

            if target_indices.is_empty() && query_indices.is_empty() {
                break;
            }

            idx += 1;
        }

        self.target_num_bases = sum_target_len;
        self.query_num_bases = sum_query_len;

        let avg_target_len = sum_target_len as f32 / self.target_num_reads as f32;
        let avg_query_len: f32 = sum_query_len as f32 / self.query_num_reads as f32;
        debug!("Target reads written to: {}", target_file.display());
        debug!("Query reads written to: {}", query_file.display());
        debug!("Total target bases: {}", sum_target_len);
        debug!("Total query bases: {}", sum_query_len);
        debug!("Average target read length: {}", avg_target_len);
        debug!("Average query read length: {}", avg_query_len);

        Ok((target_file, query_file, avg_target_len))
    }

    /// Align the query reads to the target reads and write the overlaps to a PAF file.
    fn align_reads(
        &self,
        aln_wrapper: AlignerWrapper,
        query_file: PathBuf,
        avg_target_len: f32,
    ) -> Result<(Vec<f32>, u32), LrgeError> {
        // Bounded channel to control memory usage - i.e., 10000 records in the channel at a time
        let (sender, receiver) = channel::bounded(10000);
        let aligner = Arc::clone(&aln_wrapper.aligner); // Shared reference for the producer thread
        let overlap_threshold = aln_wrapper.aligner.mapopt.min_chain_score as u32;

        // Producer: Read FASTQ records and send them to the channel
        let producer = std::thread::spawn(move || -> Result<(), LrgeError> {
            let mut fastx_reader = parse_fastx_file(query_file).map_err(|e| {
                LrgeError::FastqParseError(format!("Error parsing query FASTQ file: {e}",))
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
                            "Error parsing query FASTQ file: {e}",
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
                LrgeError::ThreadError(format!("Error setting number of threads: {e}",))
            })?;

        let estimates = Vec::with_capacity(self.query_num_reads);
        let estimates = Arc::new(Mutex::new(estimates));
        let no_mapping_count = AtomicU32::new(0);

        debug!("Aligning reads and writing overlaps to PAF file...");
        // Consumer: Process records from the channel in parallel
        pool.install(|| -> Result<(), LrgeError> {
            receiver
                .into_iter()
                .par_bridge() // Parallelize the processing
                .try_for_each(|record| -> Result<(), LrgeError> {
                    let io::Message::Data((rid, seq)) = record;
                    trace!("Processing read: {}", String::from_utf8_lossy(&rid));

                    let qname = CString::new(rid).map_err(|e| {
                        LrgeError::MapError(format!("Error converting read ID to CString: {e}",))
                    })?;

                    // Use the shared aligner to perform alignment
                    let mappings = aligner.map(&seq, Some(&qname)).map_err(|e| {
                        LrgeError::MapError(format!(
                            "Error mapping read {}: {e}",
                            String::from_utf8_lossy(qname.as_bytes()),
                        ))
                    })?;

                    let mut unique_overlaps = HashSet::new();

                    if !mappings.is_empty() {
                        {
                            let mut writer_lock = writer.lock().unwrap();
                            for mapping in &mappings {
                                // write the PafRecord to the PAF file
                                writer_lock.serialize(mapping)?;

                                if self.remove_internal
                                    && mapping.is_internal(self.max_overhang_ratio)
                                {
                                    continue;
                                }
                                unique_overlaps.insert(mapping.target_name.clone());
                            }
                        }
                    } else {
                        trace!(
                            "No overlaps found for read: {}",
                            String::from_utf8_lossy(qname.as_bytes())
                        );
                        no_mapping_count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    }

                    let est = per_read_estimate(
                        seq.len(),
                        avg_target_len,
                        self.target_num_reads,
                        unique_overlaps.len(),
                        overlap_threshold,
                    );

                    trace!(
                        "Estimate for {}: {}",
                        String::from_utf8_lossy(qname.as_bytes()),
                        est
                    );

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
            LrgeError::ThreadError(format!("Thread panicked when joining: {e:?}",))
        })??;

        debug!("Overlaps written to: {}", paf_path.to_string_lossy());

        let no_mapping_count = no_mapping_count.load(std::sync::atomic::Ordering::Relaxed);
        if no_mapping_count > 0 {
            let percent = (no_mapping_count as f32 / self.query_num_reads as f32) * 100.0;
            info!(
                "{} ({:.2}%) query read(s) did not overlap any target reads",
                no_mapping_count, percent
            );
        } else {
            debug!("All query reads overlapped with target reads");
        }

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

        Ok((estimates, no_mapping_count))
    }

    /// Align the target reads to the query reads and write the overlaps to a PAF file.
    fn align_reads_inverse(
        &self,
        aln_wrapper: AlignerWrapper,
        target_file: PathBuf,
        avg_target_len: f32,
    ) -> Result<(Vec<f32>, u32), LrgeError> {
        // Bounded channel to control memory usage - i.e., 10000 records in the channel at a time
        let (sender, receiver) = channel::bounded(10000);
        let aligner = Arc::clone(&aln_wrapper.aligner); // Shared reference for the producer thread
        let overlap_threshold = aln_wrapper.aligner.mapopt.min_chain_score as u32;

        // Producer: Read FASTQ records and send them to the channel
        let producer = std::thread::spawn(move || -> Result<(), LrgeError> {
            let mut fastx_reader = parse_fastx_file(target_file).map_err(|e| {
                LrgeError::FastqParseError(format!("Error parsing query FASTQ file: {e}",))
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
                            "Error parsing query FASTQ file: {e}",
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
                LrgeError::ThreadError(format!("Error setting number of threads: {e}",))
            })?;

        let mut read_lengths: HashMap<Vec<u8>, usize> =
            HashMap::with_capacity(self.query_num_reads);
        let mut ovlap_counter: HashMap<Vec<u8>, usize> =
            HashMap::with_capacity(self.query_num_reads);

        for i in 0..self.query_num_reads {
            unsafe {
                let qname: *mut ::std::os::raw::c_char =
                    (*((*(aln_wrapper.aligner.idx.unwrap())).seq.add(i))).name;
                let qname = std::ffi::CStr::from_ptr(qname).to_bytes().to_vec();
                let qlens: usize =
                    (*((*(aln_wrapper.aligner.idx.unwrap())).seq.add(i))).len as usize;
                // add to read_lengths
                if read_lengths.insert(qname.clone(), qlens).is_some() {
                    return Err(LrgeError::DuplicateReadIdentifier(
                        String::from_utf8_lossy(&qname).to_string(),
                    ));
                }
                // add to ovlap_counter, we insert it with 0 overlaps
                if ovlap_counter.insert(qname.clone(), 0).is_some() {
                    return Err(LrgeError::DuplicateReadIdentifier(
                        String::from_utf8_lossy(&qname).to_string(),
                    ));
                }
            }
        }

        let ovlap_counter = Arc::new(Mutex::new(ovlap_counter));

        debug!("Aligning reads and writing overlaps to PAF file...");
        // Consumer: Process records from the channel in parallel
        pool.install(|| -> Result<(), LrgeError> {
            receiver
                .into_iter()
                .par_bridge() // Parallelize the processing
                .try_for_each(|record| -> Result<(), LrgeError> {
                    let io::Message::Data((rid, seq)) = record;
                    trace!("Processing read: {}", String::from_utf8_lossy(&rid));

                    let tname: CString = CString::new(rid.clone()).map_err(|e| {
                        LrgeError::MapError(format!("Error converting read name to CString: {e}",))
                    })?;

                    // Use the shared aligner to perform alignment
                    let mappings = aligner.map(&seq, Some(&tname)).map_err(|e| {
                        LrgeError::MapError(format!(
                            "Error mapping read {}: {e}",
                            String::from_utf8_lossy(&rid),
                        ))
                    })?;

                    {
                        if !mappings.is_empty() {
                            let mut writer_lock = writer.lock().unwrap();
                            let mut ovlap_counter_lock = ovlap_counter.lock().unwrap();
                            let mut unique_overlaps: HashSet<Vec<u8>> = HashSet::new();
                            let mut overhang: i32;
                            let mut maplen: i32;

                            for mapping in &mappings {
                                // write the PafRecord to the PAF file
                                writer_lock.serialize(mapping)?;

                                if unique_overlaps.contains(&mapping.target_name) {
                                    continue;
                                }

                                if self.remove_internal {
                                    if mapping.strand == '+' {
                                        overhang =
                                            cmp::min(mapping.query_start, mapping.target_start)
                                                + cmp::min(
                                                    mapping.query_len - mapping.query_end,
                                                    mapping.target_len - mapping.target_end,
                                                );
                                    } else {
                                        overhang = cmp::min(
                                            mapping.query_start,
                                            mapping.target_len - mapping.target_end,
                                        ) + cmp::min(
                                            mapping.query_len - mapping.query_end,
                                            mapping.target_start,
                                        );
                                    }
                                    maplen = cmp::max(
                                        mapping.query_end - mapping.query_start,
                                        mapping.target_end - mapping.target_start,
                                    );
                                    if overhang > ((maplen as f32) * self.max_overhang_ratio) as i32
                                    {
                                        continue;
                                    }
                                }

                                *ovlap_counter_lock
                                    .entry(mapping.target_name.clone())
                                    .or_insert(0) += 1;
                                unique_overlaps.insert(mapping.target_name.clone());
                            }
                        }
                    }

                    Ok(())
                })?;
            Ok(())
        })?;

        // Wait for the producer to finish
        producer.join().map_err(|e| {
            LrgeError::ThreadError(format!("Thread panicked when joining: {e:?}",))
        })??;

        debug!("Overlaps written to: {}", paf_path.to_string_lossy());

        let ovlap_counter = Arc::try_unwrap(ovlap_counter)
            .unwrap()
            .into_inner()
            .unwrap();
        let no_mapping_count = AtomicU32::new(0);
        let estimates = ovlap_counter
            .par_iter()
            .map(|(rid, n_ovlaps)| {
                let est = if *n_ovlaps == 0 {
                    no_mapping_count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    trace!(
                        "No overlaps found for read: {}",
                        String::from_utf8_lossy(rid)
                    );
                    f32::INFINITY
                } else {
                    // safe to unwrap the Option here because we know the key exists
                    let read_len = read_lengths.get(rid).unwrap();
                    per_read_estimate(
                        *read_len,
                        avg_target_len,
                        self.target_num_reads,
                        *n_ovlaps,
                        overlap_threshold,
                    )
                };
                trace!("Estimate for {}: {}", String::from_utf8_lossy(rid), est);
                est
            })
            .collect();

        let no_mapping_count = no_mapping_count.load(std::sync::atomic::Ordering::Relaxed);

        if no_mapping_count > 0 {
            let percent = (no_mapping_count as f32 / self.query_num_reads as f32) * 100.0;
            info!(
                "{} ({:.2}%) read(s) did not overlap any other reads",
                no_mapping_count, percent
            );
        } else {
            debug!("All reads had at least one overlap");
        }

        Ok((estimates, no_mapping_count))
    }
}

impl Estimate for TwoSetStrategy {
    fn generate_estimates(&mut self) -> crate::Result<(Vec<f32>, u32)> {
        let (target_file, query_file, avg_target_len) = self.split_fastq()?;

        let preset = match self.platform {
            Platform::PacBio => Preset::AvaPb,
            Platform::Nanopore => Preset::AvaOnt,
        };

        if self.use_min_ref && self.target_num_bases > self.query_num_bases {
            // align target to query
            let aligner = AlignerWrapper::new(&query_file, self.threads, preset, true)?;
            self.align_reads_inverse(aligner, target_file, avg_target_len)
        } else {
            // align query to target
            let aligner = AlignerWrapper::new(&target_file, self.threads, preset, true)?;
            self.align_reads(aligner, query_file, avg_target_len)
        }
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
/// * `original` - The `Vec` to be split. This set will be consumed by the function, so it will no
///   longer be accessible after the function call.
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
