//! A strategy that compares overlaps between two sets of reads.
mod builder;

pub use self::builder::Builder;
use crate::io::FastqRecordExt;
use crate::minimap2::{Aligner, Preset};
use crate::{error::LrgeError, io, unique_random_set, Estimate, Platform};
use crossbeam_channel as channel;
use log::{debug, warn};
use needletail::{parse_fastx_file, parse_fastx_reader};
use rayon::prelude::*;
use std::collections::HashSet;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::sync::Arc;

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
    /// To customize the strategy, use the [`Builder`] interface.
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

        if n_fq_reads <= self.target_num_reads {
            let msg = format!(
                "Number of reads in FASTQ file ({}) is <= target number of reads ({})",
                n_fq_reads, self.target_num_reads
            );
            return Err(LrgeError::TooFewReadsError(msg));
        } else if n_fq_reads < n_req_reads {
            warn!(
                "Number of reads in FASTQ file ({}) is less than the sum of target and query reads ({})",
                n_fq_reads, n_req_reads
            );
            self.query_num_reads = n_fq_reads - self.target_num_reads;
            n_req_reads = n_fq_reads;
            warn!("Using {} query reads", self.query_num_reads);
        }

        let indices = unique_random_set(n_req_reads, n_fq_reads as u32, self.seed);
        let (mut target_indices, mut query_indices) =
            split_into_hashsets(indices, self.target_num_reads);

        let target_file = self.tmpdir.join("target.fastq");
        let query_file = self.tmpdir.join("query.fastq");

        let reader = io::open_file(&self.input)?;
        let mut fastx_reader = parse_fastx_reader(reader).map_err(|e| {
            LrgeError::FastqParseError(format!("Error parsing input FASTQ file: {}", e))
        })?;

        debug!("Writing target and query reads to temporary files...");
        let mut target_writer = File::create(&target_file).map(BufWriter::new)?;
        let mut query_writer = File::create(&query_file).map(BufWriter::new)?;
        let mut sum_query_len = 0;
        let mut idx: u32 = 0;
        while let Some(r) = fastx_reader.next() {
            // we can unwrap here because we know the file is valid from when we counted the records
            let record = r.unwrap();

            if target_indices.remove(&idx) {
                record.write(&mut target_writer, None).map_err(|e| {
                    LrgeError::IoError(std::io::Error::new(std::io::ErrorKind::Other, e))
                })?;
            } else if query_indices.remove(&idx) {
                record.write(&mut query_writer, None).map_err(|e| {
                    LrgeError::IoError(std::io::Error::new(std::io::ErrorKind::Other, e))
                })?;
                sum_query_len += record.num_bases();
            }

            if target_indices.is_empty() && query_indices.is_empty() {
                break;
            }

            idx += 1;
        }

        let avg_query_len = sum_query_len as f32 / self.query_num_reads as f32;
        debug!("Target reads written to: {}", target_file.display());
        debug!("Query reads written to: {}", query_file.display());
        debug!("Average query read length: {}", avg_query_len);

        Ok((target_file, query_file, avg_query_len))
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

impl Estimate for TwoSetStrategy {
    fn generate_estimates(&mut self) -> crate::Result<Vec<f32>> {
        let (target_file, query_file, _avg_query_len) = self.split_fastq()?;

        let preset = match self.platform {
            Platform::PacBio => Preset::AvaPb,
            Platform::Nanopore => Preset::AvaOnt,
        };

        let aligner = AlignerWrapper::new(&target_file, self.threads, preset, true)?;
        let _alignments = aligner.align_reads(query_file, &self.tmpdir)?;

        // for mapping in alignments {
        //     let mut fields = Vec::new();
        //     fields.push(mapping.query_name.unwrap_or(b"*".to_vec()));
        //     fields.push(
        //         mapping
        //             .query_len
        //             .map(|x| x)
        //             .unwrap_or(b"*".to_vec()),
        //     );
        //     fields.push(mapping.query_start.to_string());
        //     fields.push(mapping.query_end.to_string());
        //     fields.push(mapping.strand.to_string());
        //     fields.push(mapping.target_name.unwrap_or("*".to_string()));
        //     fields.push(mapping.target_len.to_string());
        //     fields.push(mapping.target_start.to_string());
        //     fields.push(mapping.target_end.to_string());
        //     fields.push(mapping.match_len.to_string());
        //     fields.push(mapping.block_len.to_string());
        //     fields.push(mapping.mapq.to_string());
        //     let row = fields.join("\t");
        //     println!("{}", row);
        // }

        let estimates = vec![0.0; self.target_num_reads];
        Ok(estimates)
    }
}

struct AlignerWrapper {
    aligner: Arc<Aligner>, // Shared aligner across threads
}

impl AlignerWrapper {
    fn new(
        target_file: &Path,
        threads: usize,
        preset: Preset,
        dual: bool,
    ) -> Result<Self, LrgeError> {
        let aligner = Aligner::builder()
            .preset(preset.as_bytes())
            .dual(dual)
            .with_index_threads(threads)
            .with_index(target_file, None)
            .unwrap();

        Ok(Self {
            aligner: Arc::new(aligner),
        })
    }

    fn align_reads(&self, query_file: PathBuf, tmpdir: &PathBuf) -> Result<PathBuf, LrgeError> {
        let (sender, receiver) = channel::bounded(1000); // Bounded channel to control memory usage
        let aligner = Arc::clone(&self.aligner); // Shared reference for the producer thread

        // Producer: Read FASTQ records and send them to the channel
        let producer = std::thread::spawn(move || -> Result<(), LrgeError> {
            let mut fastx_reader = parse_fastx_file(query_file).map_err(|e| {
                LrgeError::FastqParseError(format!("Error parsing query FASTQ file: {}", e))
            })?;

            while let Some(record) = fastx_reader.next() {
                match record {
                    Ok(rec) => {
                        let msg = io::Message::Data((
                            rec.read_id().to_owned(),
                            Vec::from(rec.seq().to_owned()),
                        ));
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
            // sender.send(io::Message::End).unwrap();
            // Close the channel to signal that no more records will be sent
            drop(sender);
            Ok(())
        });

        // Consumer: Process records from the channel in parallel
        // let alignments: Vec<Mapping> = receiver
        //     .into_iter()
        //     .par_bridge() // Parallelize the processing
        //     .flat_map(|record| {
        //         let (rid, seq) = match record {
        //             io::Message::Data(data) => data,
        //             io::Message::End => unimplemented!("End message not expected here"),
        //         };
        //         debug!("Processing read: {:?}", String::from_utf8_lossy(&rid));
        //         let mut qname = rid.to_owned();
        //         if !(qname.last() == Some(&0)) {
        //             qname.push(0);
        //         }
        //         // Use the shared aligner to perform alignment
        //         aligner.map(&seq, Some(&rid)).unwrap()
        //     })
        //     .collect();

        // Open the output file for writing
        let paf_path = tmpdir.join("overlaps.paf");
        let mut _file = File::create(&paf_path).map(BufWriter::new)?;

        // Consumer: Process records from the channel in parallel
        receiver
            .into_iter()
            .par_bridge() // Parallelize the processing
            .try_for_each(|record| -> Result<(), LrgeError> {
                let io::Message::Data((rid, seq)) = record;
                debug!("Processing read: {:?}", String::from_utf8_lossy(&rid));
                let mut qname = rid.to_owned();
                if qname.last() != Some(&0) {
                    qname.push(0);
                }
                // Use the shared aligner to perform alignment
                let _mappings = aligner.map(&seq, Some(&rid)).unwrap();
                // for mapping in mappings {
                //     file.write(mapping)?;
                // }
                Ok(())
            })?;

        // Wait for the producer to finish
        producer.join().unwrap();

        Ok(paf_path)
    }
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
