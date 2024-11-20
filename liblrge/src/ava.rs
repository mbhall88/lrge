//! A strategy that compares overlaps between the same set of reads - i.e., all-vs-all.
//!
//! In general, this strategy is less computationally efficient than [`TwoSetStrategy`], but it
//! is slightly more accurate - though that difference in accuracy is not statistically significant.
mod builder;

use std::collections::HashSet;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use log::{debug, warn};
use needletail::parse_fastx_reader;

use crate::{Estimate, io, Platform, unique_random_set};
use crate::error::LrgeError;
pub use self::builder::Builder;

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
                n_fq_reads,
                self.num_reads
            );
            self.num_reads = n_fq_reads;
        }

        let mut indices: HashSet<u32> = unique_random_set(self.num_reads, n_fq_reads as u32, self.seed).iter().cloned().collect();

        let out_file = self.tmpdir.join("reads.fq");
        let reader = io::open_file(&self.input)?;
        let mut fastx_reader = parse_fastx_reader(reader).map_err(|e| {
            LrgeError::FastqParseError(format!("Error parsing input FASTQ file: {}", e))
        })?;

        debug!("Writing target and query reads to temporary files...");
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
}

impl Estimate for AvaStrategy {
    fn generate_estimates(&mut self) -> crate::Result<Vec<f32>> {
        let (reads_file, sum_len) = self.subsample_reads()?;
        todo!("Implement generate_estimates for AvaStrategy");
    }
}