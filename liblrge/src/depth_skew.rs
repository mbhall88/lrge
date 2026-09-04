use std::collections::{BinaryHeap, HashSet};
use std::fmt;
use std::path::Path;
use std::sync::atomic::{AtomicU32, Ordering};
use std::time::Instant;

use crossbeam_channel as channel;
use log::debug;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use simd_minimizers::packed_seq::{PackedSeqVec, SeqVec};

use crate::io;

const KMER_LENGTH: usize = 15;
const WINDOW_WIDTH: usize = 9;
const READ_SAMPLE_DENOMINATOR: u32 = 100;
/// The fewest reads detection will settle for, when the input holds that many.
///
/// Below this the 99.9th percentile of minimizer counts moves with the draw rather than with the
/// input. `SRR26715166` holds 7,473 reads, so one in a hundred is about 70 of them, and over
/// twelve seeds its score ranged 12x to 22x around a mean of 18.0 with a standard deviation of
/// 2.9. The 16x threshold sat 0.7 standard deviations below that mean and two of the twelve seeds
/// missed a genuinely skewed input. 500 reads puts the threshold about 1.75 standard deviations
/// out. #36 refits this against the benchmark.
pub(crate) const DETECTION_READ_FLOOR: usize = 500;
/// The most the reserve will hold, whatever the reads turn out to be.
///
/// [`DETECTION_READ_FLOOR`] reads is 3.5 MB at a 7 kb read length, but read length has no upper
/// bound and five hundred ultra-long reads would be hundreds of megabytes. Past this the reserve
/// drops its largest hashes, so it holds the reads that fit rather than all of them. That leaves
/// it short of the floor on an input averaging more than about 128 kb a read, which is the trade:
/// what is scarce on an input like that is reads, and there is no bound on how much memory
/// gathering more of them could take.
const DETECTION_RESERVE_BYTES: usize = 64 << 20;
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
// Bases of sequence a scoring batch gathers before it is handed to a worker. Smaller than a
// profiling batch because a scoring batch makes a round trip: it holds the reads until the
// consumer has seen them, so its size sets what the pipeline holds in memory. A batch is at least
// one read, so it runs over this by up to the length of the read that closed it.
const SCORING_BATCH_BASES: usize = 1 << 18;
// Scoring batches a worker may have queued on either side of it. Batches hold the same number of
// bases and the workers are alike, so one is rarely more than a batch behind the others; two of
// slack absorbs the jitter of sharing a machine without holding reads that buy nothing.
const SCORING_QUEUE_BATCHES: usize = 2;
// What the pipeline holds at once, per thread: a queue on each side of a worker, and the batch it
// is scoring. See [`score_reads`].
const SCORING_PIPELINE_BATCHES: usize = 2 * SCORING_QUEUE_BATCHES + 1;

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
/// nothing beyond one random number and, while the reserve is open, a copy.
///
/// One in a hundred is too few reads on a small input, and the record count is not known until the
/// pass that feeds this one has finished. So a [`ReadReserve`] holds on to reads the draw passed
/// over, and [`finish`](Self::finish) tops the sample up from it when the draw came up short of
/// [`DETECTION_READ_FLOOR`]. The reserve closes the moment the draw reaches the floor on its own,
/// which is why a large input pays for it only over its first reads.
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
    /// Detection settles for no fewer reads than this. Only tests lower it.
    read_floor: usize,
    /// Reads the draw passed over, to top the sample up from if it comes up short.
    reserve: ReadReserve,
    scratch: MinimizerScratch,
    observed_records: usize,
    sampled_records: usize,
    sketched_minimizers: usize,
}

impl DepthSkewDetector {
    pub(crate) fn new(seed: Option<u64>) -> Self {
        Self::sampling(seed, READ_SAMPLE_DENOMINATOR, DETECTION_READ_FLOOR)
    }

    /// A detector that draws one read in `one_in` and settles for no fewer than `at_least`.
    ///
    /// The reserve is sized here and nowhere else, so it cannot end up holding too few reads to
    /// reach the floor it is meant to serve.
    fn sampling(seed: Option<u64>, one_in: u32, at_least: usize) -> Self {
        Self {
            rng: seeded_rng(seed),
            sample_denominator: one_in,
            sketch: CountMinSketch::new(),
            sample: DistinctSample::new(DISTINCT_SAMPLE_SIZE),
            read_floor: at_least,
            reserve: ReadReserve::new(at_least, reserve_salt(seed)),
            scratch: MinimizerScratch::default(),
            observed_records: 0,
            sampled_records: 0,
            sketched_minimizers: 0,
        }
    }

    /// Offer a read to the detection sample.
    ///
    /// A read that is not drawn has its minimizers left uncomputed, which is the whole point. It
    /// is copied into the reserve while that is open, and once the draw has reached the floor it
    /// costs one random number and nothing else.
    pub(crate) fn observe(&mut self, sequence: &[u8]) {
        let index = self.observed_records;
        self.observed_records += 1;

        if !self.rng.random_ratio(1, self.sample_denominator) {
            self.reserve.offer(index, sequence);
            return;
        }

        self.sketch_read(sequence);
        if self.sampled_records >= self.read_floor {
            self.reserve.close();
        }
    }

    /// Count a read's minimizers into the sketch and the bottom-k sample.
    fn sketch_read(&mut self, sequence: &[u8]) {
        self.sampled_records += 1;
        let sketch = &self.sketch;
        let sample = &mut self.sample;
        let sketched = &mut self.sketched_minimizers;
        self.scratch.for_each_sampled_minimizer(sequence, |value| {
            sketch.increment(value);
            sample.observe(value);
            *sketched += 1;
        });
    }

    pub(crate) fn finish(mut self) -> DepthDetection {
        // The reserve holds only reads the draw passed over, so nothing here is counted twice.
        let short = self.read_floor.saturating_sub(self.sampled_records);
        if short > 0 {
            let topped_up = self.reserve.take(short);
            debug!(
                "Depth detection drew {} reads, short of {}, so it took {} more from reserve",
                self.sampled_records,
                self.read_floor,
                topped_up.len()
            );
            for sequence in topped_up {
                self.sketch_read(&sequence);
            }
        }

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
        Self::sampling(seed, 1, DETECTION_READ_FLOOR)
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
    let started = Instant::now();
    let sketch = CountMinSketch::new();
    let draw_sample = sample.is_none();

    let (read_result, drawn) = std::thread::scope(|scope| {
        let (sender, receiver) = channel::bounded::<Vec<Vec<u8>>>(2 * threads);
        let workers: Vec<_> = (0..threads)
            .map(|_| {
                let receiver = receiver.clone();
                let sketch = &sketch;
                scope.spawn(move || {
                    let mut scratch = MinimizerScratch::default();
                    let mut sample = draw_sample.then(|| DistinctSample::new(DISTINCT_SAMPLE_SIZE));
                    for batch in receiver {
                        for sequence in &batch {
                            scratch.for_each_sampled_minimizer(sequence, |value| {
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
    debug!("Built the depth profile in {:.2?}", started.elapsed());

    Ok(DepthProfile {
        sketch,
        median_count,
    })
}

/// The buffers a pass reuses to walk the minimizers of one read after another.
///
/// Both buffers are scratch, and reusing them from read to read is most of why walking every read
/// of an input is affordable at all. They live together because a walk needs both, so a pass holds
/// one of these rather than two loose buffers, and each thread of a parallel pass takes its own.
#[derive(Default)]
struct MinimizerScratch {
    /// Where `simd_minimizers` puts the positions it finds.
    positions: Vec<u32>,
    /// The piece of read being walked, two bits to a base.
    packed: PackedSeqVec,
}

impl MinimizerScratch {
    /// Feed the sampled canonical minimizers of `sequence` to `observe`.
    ///
    /// Only the sampled fraction of the minimizers is passed on, so the pass that builds the
    /// sketch and the pass that reads it agree on which minimizers exist.
    ///
    /// Each piece is packed before it is walked. Packing reads the piece once and writes a quarter
    /// as many bytes, after which recovering a minimizer's value costs a load, a shift and a mask.
    /// Walking the ASCII bytes instead leaves every value to be re-read from the k-mer and packed
    /// a base at a time, forward and reverse-complement, and a base belongs to several minimizers.
    /// That is where most of a normalizing run used to go. `simd_minimizers` asks for this: its
    /// input-types documentation says ASCII DNA should usually be converted to a `PackedSeqVec`
    /// first.
    fn for_each_sampled_minimizer(&mut self, sequence: &[u8], mut observe: impl FnMut(u64)) {
        // Only the pieces a walk will use are packed, so the non-ACGT runs a read is split by, and
        // the pieces too short to hold a window, cost nothing to pack.
        for segment in walkable_segments(sequence) {
            // `push_ascii` maps whatever byte it is given into two bits without complaint, and
            // `as_slice` indexes past the packed bases into padding the vector only has once
            // something has been pushed. Both are safe here only because the segments are non-empty
            // and hold nothing but ACGT.
            self.packed.clear();
            self.packed.push_ascii(segment);
            self.positions.clear();
            let minimizers = simd_minimizers::canonical_minimizers(KMER_LENGTH, WINDOW_WIDTH)
                .run(self.packed.as_slice(), &mut self.positions);
            for value in minimizers.values_u64() {
                if is_sampled(value) {
                    observe(value);
                }
            }
        }
    }
}

/// The pieces of `sequence` a minimizer walk has anything to say about.
///
/// A run of non-ACGT bases splits the read, because a k-mer that spans one has no two-bit
/// encoding. A piece too short to hold a window has no minimizers and is dropped, which leaves
/// every piece that comes back non-empty and ACGT throughout.
fn walkable_segments(sequence: &[u8]) -> impl Iterator<Item = &[u8]> {
    sequence
        .split(|base| !matches!(base, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
        .filter(|segment| segment.len() >= KMER_LENGTH + WINDOW_WIDTH - 1)
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
}

impl DepthProfile {
    /// A scorer that reads this profile, with scratch of its own.
    pub(crate) fn scorer(&self) -> RetentionScorer<'_> {
        RetentionScorer {
            profile: self,
            scratch: MinimizerScratch::default(),
            counts: Vec::new(),
        }
    }
}

/// Scores reads against a finished [`DepthProfile`].
///
/// The scratch a scoring pass reuses lives here rather than on the profile, so that threads
/// scoring against one profile each take their own.
pub(crate) struct RetentionScorer<'a> {
    profile: &'a DepthProfile,
    scratch: MinimizerScratch,
    counts: Vec<u32>,
}

impl RetentionScorer<'_> {
    /// The probability normalization keeps a read of this sequence.
    pub(crate) fn probability(&mut self, sequence: &[u8]) -> f64 {
        let sketch = &self.profile.sketch;
        let counts = &mut self.counts;
        counts.clear();
        self.scratch.for_each_sampled_minimizer(sequence, |value| {
            counts.push(sketch.estimate(value));
        });

        if counts.is_empty() {
            return 1.0;
        }

        let middle = (counts.len() - 1) / 2;
        let (_, read_depth, _) = counts.select_nth_unstable(middle);
        let read_depth = (*read_depth).max(1);
        let target = self
            .profile
            .median_count
            .saturating_mul(RETENTION_TARGET_MULTIPLIER);
        (target as f64 / read_depth as f64).min(1.0)
    }
}

/// A read, and the probability the profile keeps it.
pub(crate) struct ScoredRead {
    pub(crate) id: Vec<u8>,
    pub(crate) sequence: Vec<u8>,
    pub(crate) retention_probability: f64,
}

/// A read on its way to a scorer.
///
/// The identifier rides along even though only the buffered selection path ever reads it. It is
/// well under a percent of the bytes a read's sequence takes, and the pipeline has to copy the
/// sequence either way, so carrying it is what lets one pipeline serve both paths.
struct Read {
    id: Vec<u8>,
    sequence: Vec<u8>,
}

/// Score every read of `input` for retention, and hand each to `consume` in record order.
///
/// A read's retention probability depends only on that read and the finished profile, so scoring
/// parallelises. The reservoir that consumes the scores does not: it draws from a seeded stream in
/// record order, so a seed only keeps its meaning if it is offered the same probabilities in the
/// same order. Having workers score ahead of a sequential consumer gives both.
///
/// Batches go to the workers in turn and are collected back in the same turn, so record order
/// needs no reconstructing and a worker that dies is noticed rather than waited on forever. It
/// also bounds what is held: a worker takes another batch only once it has somewhere to put the
/// one it has, so the pipeline holds [`SCORING_PIPELINE_BATCHES`] batches a thread, a couple of
/// megabytes a thread whatever the input. The read that closed a batch and the reads' identifiers
/// ride over that, so the figure logged below is approximate. These are reads in transit rather
/// than reads selected, so they sit outside what `--max-read-buffer` caps.
///
/// The cost of collecting in turn is that a worker held up by the machine holds up the pass once
/// its queue fills, where handing every batch to whoever is free would not. Batches are equal in
/// bases and the work is the same for each, so the workers finish in much the same time and the
/// queues carry the rest.
pub(crate) fn score_reads(
    input: &Path,
    profile: &DepthProfile,
    threads: usize,
    mut consume: impl FnMut(usize, ScoredRead),
) -> crate::Result<()> {
    let threads = threads.max(1);
    debug!(
        "Scoring reads for retention on {threads} thread(s), holding up to about {} bytes of reads in flight...",
        SCORING_PIPELINE_BATCHES * threads * SCORING_BATCH_BASES
    );
    let started = Instant::now();

    let read_result = std::thread::scope(|scope| {
        let mut queues = Vec::with_capacity(threads);
        let mut scored = Vec::with_capacity(threads);
        for _ in 0..threads {
            let (queue, unscored) = channel::bounded::<Vec<Read>>(SCORING_QUEUE_BATCHES);
            let (done, finished) = channel::bounded::<Vec<ScoredRead>>(SCORING_QUEUE_BATCHES);
            queues.push(queue);
            scored.push(finished);
            scope.spawn(move || {
                let mut scorer = profile.scorer();
                for reads in unscored {
                    let batch: Vec<ScoredRead> = reads
                        .into_iter()
                        .map(|read| {
                            let retention_probability = scorer.probability(&read.sequence);
                            ScoredRead {
                                id: read.id,
                                sequence: read.sequence,
                                retention_probability,
                            }
                        })
                        .collect();
                    // A send fails only once the consumer has gone, which ends the pass anyway.
                    if done.send(batch).is_err() {
                        break;
                    }
                }
            });
        }

        let reader = scope.spawn(move || {
            let mut reads = Vec::new();
            let mut batch_bases = 0_usize;
            let mut turn = 0_usize;
            let read_result = io::iter_records(input, |id, sequence| {
                batch_bases += sequence.len();
                reads.push(Read {
                    id: id.to_vec(),
                    sequence: sequence.to_vec(),
                });
                if batch_bases >= SCORING_BATCH_BASES {
                    batch_bases = 0;
                    hand_over(&queues, &mut turn, std::mem::take(&mut reads))?;
                }
                Ok(())
            });
            if read_result.is_ok() && !reads.is_empty() {
                hand_over(&queues, &mut turn, reads)?;
            }
            read_result
        });

        let mut index = 0_usize;
        let mut turn = 0_usize;
        // Each worker gives its batches up in the order it was handed them, and they were handed
        // out in turn, so taking them back in turn walks the input in record order. A worker with
        // nothing left to give ends the pass, whether it ran out of reads or died holding some.
        while let Ok(batch) = scored[turn % threads].recv() {
            turn += 1;
            for read in batch {
                consume(index, read);
                index += 1;
            }
        }
        // A worker the loop above stopped short of may still be holding a finished batch it cannot
        // put down. Letting go of every receiver frees it, and freeing it frees the reader waiting
        // to hand it more, so the join below has something to report rather than something to wait
        // on.
        drop(scored);

        reader.join().expect("read scoring thread panicked")
    });
    read_result?;

    debug!("Scored reads for retention in {:.2?}", started.elapsed());
    Ok(())
}

/// Hand a batch to the worker whose turn it is.
///
/// A closed queue means that worker has gone, which on this side of the pass only happens when it
/// died. Stopping here keeps the reader from walking the rest of the input into a pipeline nothing
/// is draining, and the panic the scope carries out says what went wrong.
fn hand_over(
    queues: &[channel::Sender<Vec<Read>>],
    turn: &mut usize,
    reads: Vec<Read>,
) -> std::io::Result<()> {
    let queue = &queues[*turn % queues.len()];
    *turn += 1;
    queue.send(reads).map_err(|_| {
        std::io::Error::new(
            std::io::ErrorKind::BrokenPipe,
            "a read scoring thread stopped early",
        )
    })
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

/// The reads detection falls back on when its draw comes up short of the floor.
///
/// Holds at most `capacity` reads, the ones with the smallest hashes among those the draw passed
/// over, so which reads it holds turns on the seed rather than on their order in the file. It
/// stops taking reads and releases what it holds the moment the draw reaches the floor on its own,
/// so an input large enough to get there carries it only over its first reads.
///
/// The reads are held as sequences rather than as minimizers because most of them are never used:
/// a read enters and is evicted again, and computing its minimizers on the way in would roughly
/// double what detection costs a large input. Evicting reuses the buffer it frees, so what is held
/// is the reads kept rather than one allocation for every read that ever passed through, and
/// [`DETECTION_RESERVE_BYTES`] caps that however long the reads are.
struct ReadReserve {
    capacity: usize,
    byte_cap: usize,
    salt: u64,
    /// Ordered so the largest hash, the next to be evicted, is on top.
    entries: BinaryHeap<(u64, usize, Vec<u8>)>,
    /// What the buffers held below have reserved, which is what the cap is against.
    bytes: usize,
    open: bool,
}

impl ReadReserve {
    fn new(capacity: usize, salt: u64) -> Self {
        Self {
            capacity,
            byte_cap: DETECTION_RESERVE_BYTES,
            salt,
            entries: BinaryHeap::with_capacity(capacity),
            bytes: 0,
            open: capacity > 0,
        }
    }

    /// Offer a read the draw passed over.
    fn offer(&mut self, index: usize, sequence: &[u8]) {
        if !self.open {
            return;
        }

        let key = (mix64(index as u64 ^ self.salt), index);
        if self.entries.len() < self.capacity {
            let held = sequence.to_vec();
            self.bytes += held.capacity();
            self.entries.push((key.0, key.1, held));
        } else if self
            .entries
            .peek()
            .is_some_and(|(hash, index, _)| (*hash, *index) > key)
        {
            // Take the evicted read's buffer rather than allocating another. Far more reads pass
            // through the reserve than sit in it, and allocating for each one leaves the process
            // holding the high-water mark of all of them rather than of the ones it kept.
            let (_, _, mut buffer) = self.entries.pop().expect("the reserve is full");
            self.bytes -= buffer.capacity();
            buffer.clear();
            buffer.extend_from_slice(sequence);
            self.bytes += buffer.capacity();
            self.entries.push((key.0, key.1, buffer));
        } else {
            return;
        }

        // Drop the largest hashes until what is held fits. The smallest hashes that fit are still
        // an unbiased sample of the reads offered, just a smaller one, and lowering the capacity
        // with them keeps later reads on the branch that reuses a buffer instead of allocating.
        while self.bytes > self.byte_cap && self.entries.len() > 1 {
            let (_, _, dropped) = self.entries.pop().expect("more than one read is held");
            self.bytes -= dropped.capacity();
            self.capacity = self.entries.len();
        }
    }

    /// Stop taking reads and release the ones held.
    fn close(&mut self) {
        self.open = false;
        self.entries = BinaryHeap::new();
        self.bytes = 0;
    }

    /// Give up `wanted` of the reads held, smallest hash first.
    fn take(&mut self, wanted: usize) -> Vec<Vec<u8>> {
        self.bytes = 0;
        let mut held: Vec<_> = std::mem::take(&mut self.entries).into_vec();
        held.sort_unstable_by_key(|(hash, index, _)| (*hash, *index));
        held.truncate(wanted);
        held.into_iter().map(|(_, _, sequence)| sequence).collect()
    }
}

/// The salt that decides which reads the reserve holds.
///
/// Drawn separately from the detector's own generator, so that holding reads in reserve does not
/// shift the stream the draw comes from. Were it to shift, every input large enough to ignore the
/// reserve would still sample different reads than it did before the reserve existed.
fn reserve_salt(seed: Option<u64>) -> u64 {
    match seed {
        Some(seed) => mix64(seed ^ 0xa076_1d64_78bd_642f),
        None => rand::rng().random(),
    }
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

    /// Run `detector` over every record of `input`.
    fn run_detection(input: &std::path::Path, mut detector: DepthSkewDetector) -> DepthDetection {
        io::count_records(input, |sequence| {
            detector.observe(sequence);
            Ok(())
        })
        .unwrap();
        detector.finish()
    }

    /// The minimizer sample a detection pass over `input` leaves behind.
    fn detect(input: &std::path::Path) -> DistinctSample {
        run_detection(input, DepthSkewDetector::sampling_all(Some(3))).sample
    }

    /// A file of reads drawn from two sequences at different depths, so a minimizer's count
    /// depends on which reads were sampled rather than being 1 for all of them.
    fn two_depth_fasta(records: usize) -> tempfile::NamedTempFile {
        let chromosome = pseudo_random_dna(20_000, 1);
        let element = pseudo_random_dna(2_000, 2);
        let mut input = tempfile::NamedTempFile::new().unwrap();
        for index in 0..records {
            let source = if index % 2 == 0 {
                &chromosome
            } else {
                &element
            };
            let start = index * 137 % source.len();
            let read: Vec<u8> = (0..800)
                .map(|offset| source[(start + offset) % source.len()])
                .collect();
            writeln!(input, ">read{index}\n{}", String::from_utf8(read).unwrap()).unwrap();
        }
        input.flush().unwrap();
        input
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

    /// The minimizer values a walk over the ASCII bytes of `sequence` would yield.
    ///
    /// This is the shape the pass had before it packed the read, kept here as the reference the
    /// packed walk is measured against. It splits the read with the same [`walkable_segments`] the
    /// pass does, so the two differ in how a value is read back and in nothing else.
    fn ascii_minimizer_values(sequence: &[u8]) -> Vec<u64> {
        use simd_minimizers::packed_seq::AsciiSeq;

        let mut values = Vec::new();
        let mut positions = Vec::new();
        for segment in walkable_segments(sequence) {
            positions.clear();
            let minimizers = simd_minimizers::canonical_minimizers(KMER_LENGTH, WINDOW_WIDTH)
                .run(AsciiSeq(segment), &mut positions);
            values.extend(minimizers.values_u64().filter(|value| is_sampled(*value)));
        }
        values
    }

    /// The minimizer values `scratch` yields for `sequence`.
    fn packed_minimizer_values(scratch: &mut MinimizerScratch, sequence: &[u8]) -> Vec<u64> {
        let mut values = Vec::new();
        scratch.for_each_sampled_minimizer(sequence, |value| values.push(value));
        values
    }

    /// Reading a minimizer's value out of a packed copy of the read has to give the value reading
    /// it out of the ASCII bytes gives. A value that moves lands in a different counter, and every
    /// count the profile is built from moves with it.
    ///
    /// The cases cover what packing has to get right and a walk over ASCII never had to: a length
    /// that is not a whole number of packed bytes, lowercase bases, and reads that non-ACGT runs
    /// split into pieces. Some of those pieces are too short to hold a window and are skipped, and
    /// skipping one must not disturb the pieces after it. One scratch walks every case, in an
    /// order that puts shorter reads after longer ones, because the buffers a walk leaves behind
    /// are the next walk's.
    #[test]
    fn packing_a_read_leaves_its_minimizer_values_alone() {
        let mut cases: Vec<(String, Vec<u8>)> = Vec::new();

        // Four bases share a packed byte, so a read of each length modulo four packs a different
        // number of bases into its last one.
        for length in 4_000..4_004 {
            cases.push((format!("{length} bases"), pseudo_random_dna(length, 5)));
        }

        let read = pseudo_random_dna(3_000, 6);
        cases.push(("lowercase".to_string(), read.to_ascii_lowercase()));
        cases.push(("uppercase".to_string(), read));

        let head = pseudo_random_dna(1_500, 7);
        let tail = pseudo_random_dna(1_500, 8);
        let stub = pseudo_random_dna(KMER_LENGTH + WINDOW_WIDTH - 2, 9);
        let mut split = head.clone();
        split.extend_from_slice(b"NNNNN");
        split.extend_from_slice(&stub);
        split.push(b'N');
        split.extend_from_slice(&tail);
        cases.push(("split by an N run, around a stub".to_string(), split));

        let mut edged = vec![b'N'];
        edged.extend_from_slice(&head);
        edged.push(b'n');
        cases.push(("N at each end".to_string(), edged));

        cases.push(("too short for a window".to_string(), stub));
        cases.push(("empty".to_string(), Vec::new()));

        let mut scratch = MinimizerScratch::default();
        let mut compared = 0;
        for (name, sequence) in cases {
            let expected = ascii_minimizer_values(&sequence);
            compared += expected.len();
            assert_eq!(
                packed_minimizer_values(&mut scratch, &sequence),
                expected,
                "{name} gave different minimizer values packed than as ASCII"
            );
        }

        // Two walks that both find nothing agree on nothing, which would prove nothing.
        assert!(compared > 100, "only {compared} values were compared");
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

    /// One read in a hundred of a small input is too few for the 99.9th percentile to sit still,
    /// so the sample is topped up to the floor from reads the draw passed over.
    #[test]
    fn a_small_input_is_topped_up_to_the_floor() {
        let input = fasta(2_000, 400);

        let report = run_detection(input.path(), DepthSkewDetector::new(Some(11))).report;

        assert_eq!(report.sampled_records, DETECTION_READ_FLOOR);
    }

    /// An input smaller than the floor has every read sampled, and the top-up cannot double count
    /// a read the draw already took.
    #[test]
    fn an_input_below_the_floor_has_every_read_sampled() {
        let records = 300;
        let input = fasta(records, 400);

        let report = run_detection(input.path(), DepthSkewDetector::new(Some(11))).report;

        assert_eq!(report.sampled_records, records);
    }

    /// An input whose own draw clears the floor must sample exactly what it sampled before the
    /// floor existed. Holding reads in reserve must not disturb the draw, or every large input's
    /// verdict would move.
    #[test]
    fn a_draw_that_clears_the_floor_is_left_alone() {
        let input = fasta(6_000, 400);
        let report = |floor| {
            run_detection(
                input.path(),
                DepthSkewDetector::sampling(Some(11), 10, floor),
            )
            .report
        };

        let floored = report(DETECTION_READ_FLOOR);
        let unfloored = report(0);

        assert!(
            floored.sampled_records > DETECTION_READ_FLOOR,
            "the draw took {} reads, which does not clear the floor",
            floored.sampled_records
        );
        assert_eq!(floored.sampled_records, unfloored.sampled_records);
        assert_eq!(floored.score, unfloored.score);
        assert_eq!(floored.skewed, unfloored.skewed);
    }

    /// Reads have no bounded length, so the reserve has to stop at what it holds rather than at
    /// how many it holds. What survives must still be the smallest hashes of everything offered,
    /// or the top-up would favour whichever reads arrived first.
    #[test]
    fn the_reserve_stops_at_its_byte_cap_and_stays_a_bottom_k_sample() {
        let offered = 2_000;
        let salt = 7;
        let mut reserve = ReadReserve::new(DETECTION_READ_FLOOR, salt);
        reserve.byte_cap = 1 << 16;
        let read = vec![b'A'; 4_000];

        for index in 0..offered {
            reserve.offer(index, &read);
        }

        assert!(
            reserve.bytes <= reserve.byte_cap,
            "the reserve holds {} bytes against a cap of {}",
            reserve.bytes,
            reserve.byte_cap
        );
        assert!(!reserve.entries.is_empty());
        assert!(reserve.entries.len() < DETECTION_READ_FLOOR);

        let mut kept: Vec<u64> = reserve.entries.iter().map(|(hash, _, _)| *hash).collect();
        kept.sort_unstable();
        let mut all: Vec<u64> = (0..offered as u64)
            .map(|index| mix64(index ^ salt))
            .collect();
        all.sort_unstable();
        assert_eq!(kept, all[..kept.len()]);
    }

    /// A seed has to mean one sample, top-up included, or a rerun could change the verdict.
    /// Reseeding has to redraw it, or the floor would hand every run the same reads.
    #[test]
    fn the_topped_up_sample_is_the_same_for_a_seed_and_differs_between_seeds() {
        let input = two_depth_fasta(2_000);
        let report = |seed| run_detection(input.path(), DepthSkewDetector::new(Some(seed))).report;

        assert_eq!(report(11).score, report(11).score);
        let baseline = report(11).score;
        assert!(
            (12..20).any(|seed| report(seed).score != baseline),
            "eight reseeds all drew a sample scoring {baseline:?}"
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

        let single = profile(input.path(), 1, None).unwrap();
        let many = profile(input.path(), 4, None).unwrap();

        assert_eq!(single.median_count, many.median_count);
        let mut single_scorer = single.scorer();
        let mut many_scorer = many.scorer();
        for read in &reads {
            assert_eq!(
                single_scorer.probability(read),
                many_scorer.probability(read)
            );
        }
    }

    #[test]
    fn profiling_scores_reads_against_the_detector_sample() {
        let input = fasta(200, 2_000);

        let profile = profile(input.path(), 2, Some(detect(input.path()))).unwrap();

        assert!(profile.median_count >= 1);
    }

    /// Reads of uneven length, enough of them to fill several batches, with one long enough on its
    /// own to close the batch it lands in, so batches of very different sizes travel the pipeline
    /// together and the workers do not get equal shares of the work.
    fn uneven_fasta() -> tempfile::NamedTempFile {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        for index in 0..600_u64 {
            let length = if index == 137 {
                SCORING_BATCH_BASES + 4_096
            } else {
                500 + (index as usize % 7) * 900
            };
            writeln!(
                input,
                ">read{index}\n{}",
                String::from_utf8(pseudo_random_dna(length, index + 1)).unwrap()
            )
            .unwrap();
        }
        input.flush().unwrap();
        input
    }

    /// The scores the reservoir is offered, and the order it is offered them in, decide what a
    /// seed selects. Spreading the scoring across threads has to leave both alone.
    #[test]
    fn scoring_reads_in_parallel_matches_a_sequential_scorer() {
        let input = uneven_fasta();
        let profile = profile(input.path(), 1, None).unwrap();
        let mut scorer = profile.scorer();
        let mut expected = Vec::new();
        io::iter_records(input.path(), |id, sequence| {
            let probability = scorer.probability(sequence);
            expected.push((id.to_vec(), sequence.to_vec(), probability));
            Ok(())
        })
        .unwrap();
        assert_eq!(expected.len(), 600);

        for threads in [1, 3, 8] {
            let mut seen = Vec::new();
            score_reads(input.path(), &profile, threads, |index, read| {
                seen.push((index, read.id, read.sequence, read.retention_probability));
            })
            .unwrap();

            assert_eq!(seen.len(), expected.len(), "at {threads} thread(s)");
            for (position, (index, id, sequence, probability)) in seen.into_iter().enumerate() {
                assert_eq!(index, position, "at {threads} thread(s)");
                assert_eq!(
                    (id, sequence, probability),
                    expected[position],
                    "read {position} at {threads} thread(s)"
                );
            }
        }
    }
}
