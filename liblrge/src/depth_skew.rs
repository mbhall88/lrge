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
// One minimizer in 2^MINIMIZER_SAMPLE_SHIFT reaches the sketch. The gate is a hash of the
// minimizer, so a minimizer is either always counted or never counted: a sampled minimizer ends up
// with the count it would have had anyway, and every pass over the input agrees on which
// minimizers those are.
//
// What the sketch holds is dominated by sequencing error. A gigabase of long reads carries on the
// order of 200 million minimizer occurrences, and the great majority of the distinct values among
// them are singletons that exist in one read only. Spread over 4,194,304 counters, four to a
// value, that puts about 190 in every counter before a single genuine count arrives, and a
// minimizer whose true depth is 40 reads back as 200. Counting a quarter of them drops that floor
// to about 48 and the same minimizer reads back as 77. The median depth the profile is built on
// falls with it, which is why this change moves estimates. A quarter of a read is still hundreds
// of minimizers to take a median over, so the per-read depth loses little.
//
// A smaller fraction would go further, at the cost of short reads: a 500 base read has around a
// hundred minimizers and an eighth of that is thin. #36 refits this against the full benchmark.
const MINIMIZER_SAMPLE_SHIFT: u32 = 2;
/// Salt for the hash that both draws the minimizer sample and orders [`DistinctSample`].
///
/// Sharing one hash is what keeps the two consistent. The sample is the minimizers whose hash falls
/// below a threshold and [`DistinctSample`] keeps the smallest hashes, so as long as the gate keeps
/// at least [`DISTINCT_SAMPLE_SIZE`] of them the bottom-k sample is the same set it would have been
/// had the gate never been there. Every input a genome is worth estimating from clears that by
/// orders of magnitude. An input that does not hands the quantile a sample the gate has shrunk, and
/// one small enough to fall under [`MIN_DISTINCT_MINIMIZERS`] goes unassessed where it once got a
/// verdict.
const SAMPLE_SALT: u64 = 0x243f_6a88_85a3_08d3;
/// Salt for the hash that places a minimizer in the sketch, independent of [`SAMPLE_SALT`].
const SKETCH_SALT: u64 = 0x9e37_79b9_7f4a_7c15;
// Every counter a minimizer touches sits in one 64-byte block, so a lookup that misses cache misses
// once rather than SKETCH_LANES times. The block is divided into lanes and one counter is drawn
// from each, which keeps the counters distinct; four independent draws over sixteen shared counters
// would repeat one about a third of the time and waste the probe.
const SKETCH_LANES: usize = 4;
const SKETCH_LANE_WIDTH: usize = 4;
const SKETCH_BLOCK_COUNTERS: usize = SKETCH_LANES * SKETCH_LANE_WIDTH;
// 2^18 blocks of sixteen u32 counters is the same 16 MB table, and the same 4,194,304 counters,
// that four rows of 2^20 gave. Only where the counters sit has changed.
//
// Sharing one block does cost a little accuracy, because two minimizers that collide there collide
// in a correlated way rather than independently per row. Simulated at the load a gigabase of long
// reads puts on the table, a minimizer of true count 40 reads back with mean 201 from four
// independent rows and 205 from blocks, and the 99.9th percentile moves from 258 to 289. The gate
// above is what pays for that and more: at a quarter of the load the same minimizer reads back at
// 82. Shrinking the table would make it cache resident and hand the rest of that back.
const SKETCH_BLOCKS: usize = 1 << 18;
const SKETCH_BLOCK_BITS: u32 = SKETCH_BLOCKS.trailing_zeros();
const SKETCH_LANE_BITS: u32 = SKETCH_LANE_WIDTH.trailing_zeros();
// Both the bit masks and the shifts above assume powers of two, and one hash has to supply the
// block index and every lane's counter index.
const _: () = assert!(SKETCH_BLOCKS.is_power_of_two() && SKETCH_LANE_WIDTH.is_power_of_two());
const _: () = assert!(SKETCH_BLOCK_BITS + SKETCH_LANES as u32 * SKETCH_LANE_BITS <= 64);
// A shift of zero would ask for a shift by the full width of the hash.
const _: () = assert!(MINIMIZER_SAMPLE_SHIFT > 0 && MINIMIZER_SAMPLE_SHIFT < u64::BITS);
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
    sketched_minimizers: usize,
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
            sketched_minimizers: 0,
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
        let sketched = &mut self.sketched_minimizers;
        for_each_sampled_minimizer(sequence, &mut self.positions, |value| {
            sketch.increment(value);
            sample.observe(value);
            *sketched += 1;
        });
    }

    pub(crate) fn finish(self) -> DepthDetection {
        debug!(
            "Depth detection sketched {} minimizers from {} sampled reads",
            self.sketched_minimizers, self.sampled_records
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

/// Count every sampled minimizer of every read, in parallel, and reduce that to a depth profile.
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
                            for_each_sampled_minimizer(sequence, &mut positions, |value| {
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

/// Feed the sampled canonical minimizers of `sequence` to `observe`.
///
/// A run of non-ACGT bases splits the read; a piece too short to hold a window has no minimizers.
/// Of the minimizers that remain, only the sampled fraction is passed on, so the pass that builds
/// the sketch and the pass that reads it agree on which minimizers exist.
fn for_each_sampled_minimizer(
    sequence: &[u8],
    positions: &mut Vec<u32>,
    mut observe: impl FnMut(u64),
) {
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
            if is_sampled(value) {
                observe(value);
            }
        }
    }
}

/// The hash that both draws the minimizer sample and orders [`DistinctSample`].
fn sample_hash(value: u64) -> u64 {
    mix64(value ^ SAMPLE_SALT)
}

/// Whether this minimizer is one of the one in `2^`[`MINIMIZER_SAMPLE_SHIFT`] the sketch counts.
fn is_sampled(value: u64) -> bool {
    sample_hash(value) >> (u64::BITS - MINIMIZER_SAMPLE_SHIFT) == 0
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
        for_each_sampled_minimizer(sequence, &mut self.positions, |value| {
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

/// One cache line of counters.
#[repr(C, align(64))]
struct SketchBlock([AtomicU32; SKETCH_BLOCK_COUNTERS]);

struct CountMinSketch {
    blocks: Box<[SketchBlock]>,
}

impl CountMinSketch {
    fn new() -> Self {
        Self {
            blocks: (0..SKETCH_BLOCKS)
                .map(|_| SketchBlock(std::array::from_fn(|_| AtomicU32::new(0))))
                .collect(),
        }
    }

    /// The block this minimizer's counters live in, and which counter of it each lane holds.
    fn counters(&self, value: u64) -> (&SketchBlock, [usize; SKETCH_LANES]) {
        let (block, lanes) = sketch_slots(value);
        (&self.blocks[block], lanes)
    }

    /// Add one to each of the counters this minimizer maps to.
    ///
    /// Addition commutes, so threads sharing a sketch reach the same counts as one thread would.
    /// A counter holds at its maximum rather than wrapping, which takes a compare-and-swap rather
    /// than a plain add. Saturating needs about four billion occurrences of one minimizer, so this
    /// is defensive, but leaving the ceiling to a race would leave the counts thread dependent.
    /// The loop is written out because `fetch_update` is deprecated on newer toolchains and its
    /// replacement does not exist on the one this crate supports.
    fn increment(&self, value: u64) {
        let (block, lanes) = self.counters(value);
        for slot in lanes {
            let counter = &block.0[slot];
            let mut count = counter.load(Ordering::Relaxed);
            while count < u32::MAX {
                match counter.compare_exchange_weak(
                    count,
                    count + 1,
                    Ordering::Relaxed,
                    Ordering::Relaxed,
                ) {
                    Ok(_) => break,
                    Err(observed) => count = observed,
                }
            }
        }
    }

    fn estimate(&self, value: u64) -> u32 {
        let (block, lanes) = self.counters(value);
        lanes
            .iter()
            .map(|slot| block.0[*slot].load(Ordering::Relaxed))
            .min()
            .expect("every minimizer has a counter in every lane")
    }
}

/// Where this minimizer's counters live: which block, and which counter within each of its lanes.
fn sketch_slots(value: u64) -> (usize, [usize; SKETCH_LANES]) {
    let hash = mix64(value ^ SKETCH_SALT);
    let block = (hash as usize) & (SKETCH_BLOCKS - 1);
    let mut lanes = [0; SKETCH_LANES];
    for (lane, slot) in lanes.iter_mut().enumerate() {
        let bits = SKETCH_BLOCK_BITS + lane as u32 * SKETCH_LANE_BITS;
        *slot = lane * SKETCH_LANE_WIDTH + ((hash >> bits) as usize & (SKETCH_LANE_WIDTH - 1));
    }
    (block, lanes)
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

        let entry = (sample_hash(value), value);
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

    /// The minimizer sample a detection pass over `input` leaves behind.
    fn detect(input: &std::path::Path) -> DistinctSample {
        let mut detector = DepthSkewDetector::sampling_all(Some(3));
        io::count_records(input, |sequence| {
            detector.observe(sequence);
            Ok(())
        })
        .unwrap();
        detector.finish().sample
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
        let drawn = detector.sketched_minimizers;
        let per_read = {
            let mut one = DepthSkewDetector::sampling_all(Some(42));
            one.observe(&read);
            one.sketched_minimizers
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

    /// The sketch may overestimate a minimizer, never underestimate it. Drawing every counter from
    /// one block correlates the collisions, so this is worth pinning rather than assuming.
    #[test]
    fn the_sketch_never_undercounts() {
        let sketch = CountMinSketch::new();
        let values: Vec<u64> = (0..50_000u64).map(mix64).collect();
        for (index, value) in values.iter().enumerate() {
            for _ in 0..=index % 7 {
                sketch.increment(*value);
            }
        }

        for (index, value) in values.iter().enumerate() {
            let truth = (index % 7) as u32 + 1;
            assert!(
                sketch.estimate(*value) >= truth,
                "estimate for value {index} fell below its true count of {truth}"
            );
        }
    }

    /// One block is one cache line, and every counter a minimizer touches has to be inside it and
    /// distinct from the others. A repeated counter would spend a probe on an answer already held.
    #[test]
    fn a_minimizer_touches_one_cache_line_and_no_counter_twice() {
        assert_eq!(std::mem::size_of::<SketchBlock>(), 64);
        assert_eq!(std::mem::align_of::<SketchBlock>(), 64);

        for index in 0..20_000u64 {
            let (block, lanes) = sketch_slots(mix64(index));
            assert!(block < SKETCH_BLOCKS);
            assert!(lanes.iter().all(|slot| *slot < SKETCH_BLOCK_COUNTERS));

            let mut distinct = lanes.to_vec();
            distinct.sort_unstable();
            distinct.dedup();
            assert_eq!(
                distinct.len(),
                SKETCH_LANES,
                "lanes {lanes:?} repeat a counter"
            );
        }
    }

    #[test]
    fn the_gate_keeps_about_one_minimizer_in_four() {
        let offered = 100_000;
        let kept = (0..offered as u64)
            .filter(|index| is_sampled(mix64(*index)))
            .count();
        let expected = offered / (1 << MINIMIZER_SAMPLE_SHIFT);

        assert!(
            (expected * 9 / 10..=expected * 11 / 10).contains(&kept),
            "kept {kept} of {offered}, expected about {expected}"
        );
    }

    /// The detector's score and the profile's median depth are both read off the bottom-k sample,
    /// so sampling minimizers must not disturb it. Sharing one hash with the gate is what buys
    /// that: the smallest hashes are exactly the ones the gate keeps.
    #[test]
    fn sampling_minimizers_leaves_the_bottom_k_sample_alone() {
        let capacity = 1 << 10;
        let mut from_sample = DistinctSample::new(capacity);
        let mut from_everything = DistinctSample::new(capacity);
        let mut kept = 0;

        for index in 0..200_000u64 {
            let value = mix64(index);
            from_everything.observe(value);
            if is_sampled(value) {
                from_sample.observe(value);
                kept += 1;
            }
        }

        assert!(
            kept > capacity,
            "the gate kept {kept}, too few to fill the sample"
        );
        assert_eq!(from_sample.keys, from_everything.keys);
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

        let profile = profile(input.path(), 2, Some(detect(input.path()))).unwrap();

        assert!(profile.median_count >= 1);
    }
}
