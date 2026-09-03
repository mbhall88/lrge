use std::collections::{BinaryHeap, HashSet};
use std::fmt;
use std::path::Path;
use std::sync::atomic::{AtomicU32, Ordering};

use crossbeam_channel as channel;
use log::debug;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use simd_minimizers::packed_seq::AsciiSeq;

use crate::io;

const KMER_LENGTH: usize = 15;
const WINDOW_WIDTH: usize = 9;
const READ_SAMPLE_DENOMINATOR: u32 = 100;
const SKETCH_ROWS: usize = 4;
// Use 2^20 columns so sketch_index can replace modulo with a bit mask.
const SKETCH_WIDTH: usize = 1 << 20;
// Retain at most 2^16 distinct minimizers to bound the quantile calculation.
const DISTINCT_SAMPLE_SIZE: usize = 1 << 16;
const MIN_DISTINCT_MINIMIZERS: usize = 128;
// Express the 99.9th percentile as 999 parts per thousand to avoid float rounding.
const HIGH_QUANTILE_PER_MILLE: usize = 999;
const SKEW_THRESHOLD: f64 = 16.0;
const RETENTION_TARGET_MULTIPLIER: u32 = 2;
// Bases of sequence a profiling batch gathers before it is handed to a worker. Batching keeps the
// channel out of the way of the sketch work, which is what the profiling pass is really doing.
const PROFILE_BATCH_BASES: usize = 1 << 20;

#[derive(Clone)]
pub(crate) struct DepthSkewReport {
    /// Ratio of the 99.9th percentile minimizer count to the median count.
    pub(crate) score: Option<f64>,
    pub(crate) skewed: bool,
    pub(crate) sampled_records: usize,
}

impl DepthSkewReport {
    /// The report for a run that never asked whether the input is skewed.
    pub(crate) fn not_assessed() -> Self {
        Self {
            score: None,
            skewed: false,
            sampled_records: 0,
        }
    }
}

impl fmt::Display for DepthSkewReport {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.score {
            Some(score) => {
                let verdict = if self.skewed {
                    "Depth skew detected"
                } else {
                    "Depth skew not detected"
                };
                write!(
                    formatter,
                    "{verdict} (99.9th percentile minimizer count is {score:.2}x the median; sampled reads: {})",
                    self.sampled_records
                )
            }
            None => write!(
                formatter,
                "Depth skew not assessed (too few minimizers; sampled reads: {})",
                self.sampled_records
            ),
        }
    }
}

/// Decides whether an input is depth skewed, from a small sample of its reads.
///
/// This rides along with the pass that counts the records, so it is the only depth work an
/// unskewed run pays for. It draws minimizers from roughly one read in
/// [`READ_SAMPLE_DENOMINATOR`] and skips the rest entirely: the reads it does not sample cost it
/// nothing beyond one random number.
///
/// The sample comes from a seeded stream in record order, so this pass has to stay single
/// threaded. Parallelising it would change which reads are drawn and so what a seed means. The
/// full profile has no such constraint, which is why [`profile`] is a separate, parallel pass.
pub(crate) struct DepthSkewDetector {
    rng: StdRng,
    /// One read in `sample_denominator` is drawn. Only tests lower it.
    sample_denominator: u32,
    sketch: CountMinSketch,
    sample: DistinctSample,
    positions: Vec<u32>,
    sampled_records: usize,
    observed_minimizers: usize,
}

impl DepthSkewDetector {
    pub(crate) fn new(seed: Option<u64>) -> Self {
        Self {
            rng: seeded_rng(seed),
            sample_denominator: READ_SAMPLE_DENOMINATOR,
            sketch: CountMinSketch::new(),
            sample: DistinctSample::new(DISTINCT_SAMPLE_SIZE),
            positions: Vec::new(),
            sampled_records: 0,
            observed_minimizers: 0,
        }
    }

    /// Offer a read to the detection sample.
    ///
    /// A read that is not drawn costs one random number and nothing else. That is the whole point:
    /// the minimizers of the other ninety-nine reads in a hundred are never computed.
    pub(crate) fn observe(&mut self, sequence: &[u8]) {
        if !self.rng.random_ratio(1, self.sample_denominator) {
            return;
        }

        self.sampled_records += 1;
        let sketch = &self.sketch;
        let sample = &mut self.sample;
        let observed = &mut self.observed_minimizers;
        for_each_minimizer(sequence, &mut self.positions, |value| {
            sketch.increment(value);
            sample.observe(value);
            *observed += 1;
        });
    }

    pub(crate) fn finish(self) -> DepthDetection {
        debug!(
            "Depth detection drew {} minimizers from {} sampled reads",
            self.observed_minimizers, self.sampled_records
        );
        DepthDetection {
            report: depth_skew_report(&self.sketch, &self.sample.keys, self.sampled_records),
            sample: self.sample,
        }
    }

    #[cfg(test)]
    fn sampling_all(seed: Option<u64>) -> Self {
        Self {
            sample_denominator: 1,
            ..Self::new(seed)
        }
    }
}

/// What the detection pass leaves behind.
pub(crate) struct DepthDetection {
    pub(crate) report: DepthSkewReport,
    /// The minimizers the detector sampled. [`profile`] scores these against the full counts to
    /// find the input's median depth, so keeping them means the profiling pass does not have to
    /// draw a sample of its own.
    pub(crate) sample: DistinctSample,
}

/// Count every minimizer of every read, in parallel, and reduce that to a depth profile.
///
/// This is the expensive pass, so it only runs once normalization is going to happen. Nothing
/// here draws a random number and sketch increments commute, so the profile does not depend on
/// how the work was divided: the same input gives the same profile at any thread count.
///
/// `sample` is the minimizer sample the detector already drew. When it is `None` — a forced
/// normalization, which never runs the detector — the pass draws its own from every read.
pub(crate) fn profile(
    input: &Path,
    threads: usize,
    sample: Option<DistinctSample>,
) -> crate::Result<DepthProfile> {
    let threads = threads.max(1);
    debug!("Profiling read depth across the input on {threads} thread(s)...");
    let sketch = CountMinSketch::new();
    let draw_sample = sample.is_none();

    let (read_result, drawn) = std::thread::scope(|scope| {
        let (sender, receiver) = channel::bounded::<Vec<Vec<u8>>>(2 * threads);
        let workers: Vec<_> = (0..threads)
            .map(|_| {
                let receiver = receiver.clone();
                let sketch = &sketch;
                scope.spawn(move || {
                    let mut positions = Vec::new();
                    let mut sample = draw_sample.then(|| DistinctSample::new(DISTINCT_SAMPLE_SIZE));
                    for batch in receiver {
                        for sequence in &batch {
                            for_each_minimizer(sequence, &mut positions, |value| {
                                sketch.increment(value);
                                if let Some(sample) = &mut sample {
                                    sample.observe(value);
                                }
                            });
                        }
                    }
                    sample
                })
            })
            .collect();
        drop(receiver);

        let mut batch = Vec::new();
        let mut batch_bases = 0_usize;
        let read_result = io::count_records(input, |sequence| {
            batch_bases += sequence.len();
            batch.push(sequence.to_vec());
            if batch_bases >= PROFILE_BATCH_BASES {
                batch_bases = 0;
                // A send only fails once every worker has gone, which the join below reports.
                let _ = sender.send(std::mem::take(&mut batch));
            }
            Ok(())
        });
        if !batch.is_empty() {
            let _ = sender.send(batch);
        }
        drop(sender);

        let drawn: Vec<_> = workers
            .into_iter()
            .map(|worker| worker.join().expect("depth profiling thread panicked"))
            .collect();
        (read_result, drawn)
    });
    read_result?;

    let sample = match sample {
        Some(sample) => sample,
        // The smallest hashes of a union are the smallest hashes of its parts, so merging the
        // workers' samples gives the same minimizers however the reads were divided up.
        None => {
            let mut merged = DistinctSample::new(DISTINCT_SAMPLE_SIZE);
            for worker_sample in drawn.into_iter().flatten() {
                for value in worker_sample.keys {
                    merged.observe(value);
                }
            }
            merged
        }
    };

    let mut counts: Vec<u32> = sample
        .keys
        .iter()
        .map(|value| sketch.estimate(*value))
        .collect();
    counts.sort_unstable();
    let median_count = counts
        .get(counts.len().saturating_sub(1) / 2)
        .copied()
        .unwrap_or(1)
        .max(1);
    debug!(
        "Median depth is {median_count} across {} sampled minimizers",
        counts.len()
    );

    Ok(DepthProfile {
        sketch,
        median_count,
        positions: Vec::new(),
        counts: Vec::new(),
    })
}

/// Feed every canonical minimizer of `sequence` to `observe`.
///
/// A run of non-ACGT bases splits the read; a piece too short to hold a window has no minimizers.
fn for_each_minimizer(sequence: &[u8], positions: &mut Vec<u32>, mut observe: impl FnMut(u64)) {
    for segment in sequence
        .split(|base| !matches!(base, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
    {
        if segment.len() < KMER_LENGTH + WINDOW_WIDTH - 1 {
            continue;
        }

        positions.clear();
        let minimizers = simd_minimizers::canonical_minimizers(KMER_LENGTH, WINDOW_WIDTH)
            .run(AsciiSeq(segment), positions);
        for value in minimizers.values_u64() {
            observe(value);
        }
    }
}

fn seeded_rng(seed: Option<u64>) -> StdRng {
    match seed {
        Some(seed) => StdRng::seed_from_u64(seed),
        None => StdRng::from_rng(&mut rand::rng()),
    }
}

fn depth_skew_report(
    sketch: &CountMinSketch,
    minimizers: &HashSet<u64>,
    sampled_records: usize,
) -> DepthSkewReport {
    let mut counts: Vec<u32> = minimizers
        .iter()
        .map(|value| sketch.estimate(*value))
        .collect();
    if counts.len() < MIN_DISTINCT_MINIMIZERS {
        return DepthSkewReport {
            score: None,
            skewed: false,
            sampled_records,
        };
    }

    counts.sort_unstable();
    let median = counts[(counts.len() - 1) / 2].max(1);
    let high_index = (counts.len() * HIGH_QUANTILE_PER_MILLE)
        .div_ceil(1000)
        .saturating_sub(1);
    let score = counts[high_index] as f64 / median as f64;
    DepthSkewReport {
        score: Some(score),
        skewed: score >= SKEW_THRESHOLD,
        sampled_records,
    }
}

pub(crate) struct DepthProfile {
    sketch: CountMinSketch,
    median_count: u32,
    positions: Vec<u32>,
    counts: Vec<u32>,
}

impl DepthProfile {
    pub(crate) fn retention_probability(&mut self, sequence: &[u8]) -> f64 {
        let sketch = &self.sketch;
        let counts = &mut self.counts;
        counts.clear();
        for_each_minimizer(sequence, &mut self.positions, |value| {
            counts.push(sketch.estimate(value));
        });

        if counts.is_empty() {
            return 1.0;
        }

        let middle = (counts.len() - 1) / 2;
        let (_, read_depth, _) = counts.select_nth_unstable(middle);
        let read_depth = (*read_depth).max(1);
        let target = self
            .median_count
            .saturating_mul(RETENTION_TARGET_MULTIPLIER);
        (target as f64 / read_depth as f64).min(1.0)
    }
}

struct CountMinSketch {
    counters: Box<[AtomicU32]>,
}

impl CountMinSketch {
    fn new() -> Self {
        Self {
            counters: (0..SKETCH_ROWS * SKETCH_WIDTH)
                .map(|_| AtomicU32::new(0))
                .collect(),
        }
    }

    /// Add one to each of the counters this minimizer maps to.
    ///
    /// Addition commutes, so threads sharing a sketch reach the same counts as one thread would.
    /// A counter holds at its maximum rather than wrapping, which takes a read-modify-write rather
    /// than a plain add. Saturating needs about four billion occurrences of one minimizer, so this
    /// is defensive, but leaving the ceiling to a race would leave the counts thread dependent.
    fn increment(&self, value: u64) {
        for row in 0..SKETCH_ROWS {
            let counter = &self.counters[row * SKETCH_WIDTH + sketch_index(value, row)];
            let _ = counter.fetch_update(Ordering::Relaxed, Ordering::Relaxed, |count| {
                count.checked_add(1)
            });
        }
    }

    fn estimate(&self, value: u64) -> u32 {
        (0..SKETCH_ROWS)
            .map(|row| {
                self.counters[row * SKETCH_WIDTH + sketch_index(value, row)].load(Ordering::Relaxed)
            })
            .min()
            .unwrap_or_default()
    }
}

fn sketch_index(value: u64, row: usize) -> usize {
    let row_seed = (row as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15);
    (mix64(value ^ row_seed) as usize) & (SKETCH_WIDTH - 1)
}

fn mix64(mut value: u64) -> u64 {
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

/// The minimizers with the smallest hashes seen so far, at most `capacity` of them.
///
/// Which minimizers this keeps depends only on the set it was shown, not the order, so two of
/// these can be merged by showing one the other's keys.
pub(crate) struct DistinctSample {
    capacity: usize,
    entries: BinaryHeap<(u64, u64)>,
    keys: HashSet<u64>,
}

impl DistinctSample {
    fn new(capacity: usize) -> Self {
        Self {
            capacity,
            entries: BinaryHeap::with_capacity(capacity),
            keys: HashSet::with_capacity(capacity),
        }
    }

    fn observe(&mut self, value: u64) {
        if self.keys.contains(&value) {
            return;
        }

        let entry = (mix64(value ^ 0x243f_6a88_85a3_08d3), value);
        if self.entries.len() < self.capacity {
            self.entries.push(entry);
            self.keys.insert(value);
        } else if self.entries.peek().is_some_and(|largest| entry < *largest) {
            let (_, removed) = self.entries.pop().expect("sample is not empty");
            self.keys.remove(&removed);
            self.entries.push(entry);
            self.keys.insert(value);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn pseudo_random_dna(len: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        (0..len)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1);
                b"ACGT"[((state >> 32) & 3) as usize]
            })
            .collect()
    }

    /// A file of `records` reads of `length` bases, each read distinct from the others.
    fn fasta(records: usize, length: usize) -> tempfile::NamedTempFile {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        for index in 0..records {
            writeln!(
                input,
                ">read{index}\n{}",
                String::from_utf8(pseudo_random_dna(length, index as u64 + 1)).unwrap()
            )
            .unwrap();
        }
        input.flush().unwrap();
        input
    }

    #[test]
    fn distinguishes_high_copy_sequence_from_even_coverage() {
        let chromosome = pseudo_random_dna(20_000, 1);
        let plasmid = pseudo_random_dna(2_000, 2);

        let mut even = DepthSkewDetector::sampling_all(Some(0));
        for _ in 0..40 {
            even.observe(&chromosome);
        }
        let even = even.finish().report;

        let mut skewed = DepthSkewDetector::sampling_all(Some(0));
        for _ in 0..40 {
            skewed.observe(&chromosome);
        }
        for _ in 0..800 {
            skewed.observe(&plasmid);
        }
        let skewed = skewed.finish().report;

        assert!(!even.skewed, "even-coverage score was {:?}", even.score);
        assert!(skewed.skewed, "skewed score was {:?}", skewed.score);
        assert!(skewed.score.unwrap() > even.score.unwrap());
    }

    #[test]
    fn report_states_the_verdict_and_score() {
        let report = DepthSkewReport {
            score: Some(21.25),
            skewed: true,
            sampled_records: 87,
        };

        assert_eq!(
            report.to_string(),
            "Depth skew detected (99.9th percentile minimizer count is 21.25x the median; sampled reads: 87)"
        );
    }

    /// Detection cost has to stay proportional to the reads it samples, not to the whole input.
    /// If the full profile ever creeps back into the counting pass, this is what catches it.
    #[test]
    fn detection_draws_minimizers_only_from_the_sampled_reads() {
        let read = pseudo_random_dna(2_000, 7);
        let mut detector = DepthSkewDetector::new(Some(42));
        for _ in 0..2_000 {
            detector.observe(&read);
        }

        let sampled = detector.sampled_records;
        let drawn = detector.observed_minimizers;
        let per_read = {
            let mut one = DepthSkewDetector::sampling_all(Some(42));
            one.observe(&read);
            one.observed_minimizers
        };

        assert!(per_read > 0, "the test read has no minimizers");
        assert!(
            (5..=60).contains(&sampled),
            "sampled {sampled} of 2000 reads"
        );
        assert_eq!(
            drawn,
            sampled * per_read,
            "minimizers were drawn from reads outside the detection sample"
        );
    }

    #[test]
    fn profile_does_not_depend_on_the_thread_count() {
        // Large enough to fill several batches and to push the minimizer sample past its capacity,
        // so the workers' samples actually have to be merged.
        let input = fasta(800, 5_000);
        let reads: Vec<Vec<u8>> = (0..800)
            .map(|index| pseudo_random_dna(5_000, index + 1))
            .collect();

        let mut single = profile(input.path(), 1, None).unwrap();
        let mut many = profile(input.path(), 4, None).unwrap();

        assert_eq!(single.median_count, many.median_count);
        for read in &reads {
            assert_eq!(
                single.retention_probability(read),
                many.retention_probability(read)
            );
        }
    }

    #[test]
    fn profiling_scores_reads_against_the_detector_sample() {
        let input = fasta(200, 2_000);
        let mut detector = DepthSkewDetector::sampling_all(Some(3));
        io::count_records(input.path(), |sequence| {
            detector.observe(sequence);
            Ok(())
        })
        .unwrap();
        let detection = detector.finish();

        let profile = profile(input.path(), 2, Some(detection.sample)).unwrap();

        assert!(profile.median_count >= 1);
    }
}
