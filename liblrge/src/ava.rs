//! A strategy that compares overlaps between the same set of reads - i.e., all-vs-all.
//!
//! In general, this strategy is less computationally efficient than [`TwoSetStrategy`], but it
//! is slightly more accurate - though that difference in accuracy is not statistically significant.
mod builder;

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use crossbeam_channel as channel;
use log::{debug, trace, warn};
use needletail::{parse_fastx_file, parse_fastx_reader};
use rayon::prelude::*;

pub use self::builder::Builder;
use crate::error::LrgeError;
use crate::estimate::per_read_estimate;
use crate::io::FastqRecordExt;
use crate::minimap2::{AlignerWrapper, Preset};
use crate::{io, unique_random_set, Estimate, Platform};

pub const DEFAULT_AVA_NUM_READS: usize = 25_000;

/// A strategy that compares overlaps between two sets of reads.
///
/// The convention is to use a smaller set of query reads and a larger set of target reads. The
/// query reads are overlapped with the target reads and an estimated genome size is calculated
/// for **each query read** based on the number of overlaps it has with the target set.
pub struct AvaStrategy {
    /// Path to the FASTQ file.
    input: PathBuf,
    /// The number of reads to use in the strategy.
    num_reads: usize,
    /// The directory to which all intermediate files will be written.
    tmpdir: PathBuf,
    /// Number of threads to use with minimap2.
    threads: usize,
    /// The (optional) seed to use for randomly selecting reads.
    seed: Option<u64>,
    /// Sequencing platform of the reads.
    platform: Platform,
}

impl AvaStrategy {
    /// Create a new `AvaStrategy` with the default settings, using the given input file.
    ///
    /// To customise the strategy, use the [`Builder`] interface.
    pub fn new<P: AsRef<Path>>(input: P) -> Self {
        let builder = Builder::default();

        builder.build(input)
    }

    /// Subsample the reads in the input file to `num_reads`.
    fn subsample_reads(&mut self) -> crate::Result<(PathBuf, usize)> {
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

        if n_fq_reads < self.num_reads {
            warn!(
                "Number of reads in FASTQ file ({}) is less than the number requested ({})",
                n_fq_reads, self.num_reads
            );
            self.num_reads = n_fq_reads;
        }

        let mut indices: HashSet<u32> =
            unique_random_set(self.num_reads, n_fq_reads as u32, self.seed)
                .iter()
                .cloned()
                .collect();

        let out_file = self.tmpdir.join("reads.fq");
        let reader = io::open_file(&self.input)?;
        let mut fastx_reader = parse_fastx_reader(reader).map_err(|e| {
            LrgeError::FastqParseError(format!("Error parsing input FASTQ file: {}", e))
        })?;

        debug!("Writing subsampled reads to temporary files...");
        let mut writer = File::create(&out_file).map(BufWriter::new)?;
        let mut sum_len = 0;
        let mut idx: u32 = 0;
        while let Some(r) = fastx_reader.next() {
            // we can unwrap here because we know the file is valid from when we counted the records
            let record = r.unwrap();

            if indices.remove(&idx) {
                record.write(&mut writer, None).map_err(|e| {
                    LrgeError::IoError(std::io::Error::new(std::io::ErrorKind::Other, e))
                })?;
                sum_len += record.num_bases();
            }

            if indices.is_empty() {
                break;
            }

            idx += 1;
        }

        debug!("Reads written to: {}", out_file.display());
        debug!("Total bases written: {}", sum_len);

        Ok((out_file, sum_len))
    }

    /// Align the reads to each other, write them to a PAF file, and calculate the genome size
    /// estimate for each read
    fn align_reads(
        &self,
        aln_wrapper: AlignerWrapper,
        reads_file: PathBuf,
        sum_len: usize,
    ) -> crate::Result<Vec<f32>> {
        // Bounded channel to control memory usage - i.e., 10000 records in the channel at a time
        let (sender, receiver) = channel::bounded(25_000);
        let aligner = Arc::clone(&aln_wrapper.aligner); // Shared reference for the producer thread
        let overlap_threshold = aln_wrapper.aligner.mapopt.min_chain_score as u32;
        let read_lengths: HashMap<Vec<u8>, usize> = HashMap::with_capacity(self.num_reads);
        let read_lengths = Arc::new(Mutex::new(read_lengths));
        let read_lengths_for_producer = Arc::clone(&read_lengths);

        // Producer: Read FASTQ records and send them to the channel
        let producer = std::thread::spawn(move || -> Result<(), LrgeError> {
            let mut fastx_reader = parse_fastx_file(&reads_file).map_err(|e| {
                LrgeError::FastqParseError(format!("Error parsing FASTQ file: {}", e))
            })?;
            let read_lengths = read_lengths_for_producer;

            while let Some(record) = fastx_reader.next() {
                match record {
                    Ok(rec) => {
                        let rid = rec.read_id().to_owned();
                        let msg = io::Message::Data((rid.to_owned(), rec.seq().into_owned()));

                        {
                            // Lock the read_lengths map and insert the read length
                            let mut read_lengths_lock = read_lengths.lock().unwrap();
                            read_lengths_lock.insert(rid, rec.num_bases());
                        }

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

        let ovlap_counter: HashMap<Vec<u8>, usize> = HashMap::with_capacity(self.num_reads);
        let ovlap_counter = Arc::new(Mutex::new(ovlap_counter));
        let seen_pairs: HashSet<(Vec<u8>, Vec<u8>)> = HashSet::with_capacity(self.num_reads);
        let seen_pairs = Arc::new(Mutex::new(seen_pairs));

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

                    {
                        let mut ovlap_counter_lock = ovlap_counter.lock().unwrap();
                        if !mappings.is_empty() {
                            let mut writer_lock = writer.lock().unwrap();
                            let mut seen_pairs_lock = seen_pairs.lock().unwrap();

                            for mapping in &mappings {
                                // write the PafRecord to the PAF file
                                writer_lock.serialize(mapping)?;

                                let tname = &mapping.target_name;

                                if &rid == tname {
                                    // Skip self-overlaps
                                    continue;
                                }

                                let pair = if &rid < tname {
                                    (rid.clone(), tname.clone())
                                } else {
                                    (tname.clone(), rid.clone())
                                };
                                if seen_pairs_lock.contains(&pair) {
                                    continue;
                                } else {
                                    seen_pairs_lock.insert(pair);
                                }

                                *ovlap_counter_lock.entry(tname.clone()).or_insert(0) += 1;
                                *ovlap_counter_lock.entry(rid.clone()).or_insert(0) += 1;
                            }
                        } else {
                            debug!(
                                "No mappings found for read: {:?}",
                                String::from_utf8_lossy(&rid)
                            );
                            // if the qname is not in the ovlap_counter, we insert it with 0 overlaps
                            ovlap_counter_lock.entry(rid.clone()).or_insert(0);
                        }
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

        let ovlap_counter = Arc::try_unwrap(ovlap_counter)
            .unwrap()
            .into_inner()
            .unwrap();
        let read_lengths = Arc::try_unwrap(read_lengths).unwrap().into_inner().unwrap();
        let estimates = ovlap_counter
            .par_iter()
            .map(|(rid, n_ovlaps)| {
                // safe to unwrap the Option here because we know the key exists
                let read_len = read_lengths.get(rid).unwrap();
                let avg_read_len = sum_len as f32 / (self.num_reads - 1) as f32;
                let est = per_read_estimate(
                    *read_len,
                    avg_read_len,
                    self.num_reads - 1,
                    *n_ovlaps,
                    overlap_threshold,
                );
                trace!("Estimate for {}: {}", String::from_utf8_lossy(rid), est);
                est
            })
            .collect();

        Ok(estimates)
    }
}

impl Estimate for AvaStrategy {
    fn generate_estimates(&mut self) -> crate::Result<Vec<f32>> {
        let (reads_file, sum_len) = self.subsample_reads()?;

        let preset = match self.platform {
            Platform::PacBio => Preset::AvaPb,
            Platform::Nanopore => Preset::AvaOnt,
        };

        let aligner = AlignerWrapper::new(&reads_file, self.threads, preset, false)?;
        let estimates = self.align_reads(aligner, reads_file, sum_len)?;

        Ok(estimates)
    }
}
