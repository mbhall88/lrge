use std::collections::{BinaryHeap, HashSet};
use std::fmt;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use simd_minimizers::packed_seq::AsciiSeq;

const KMER_LENGTH: usize = 15;
const WINDOW_WIDTH: usize = 9;
const READ_SAMPLE_DENOMINATOR: u32 = 100;
const SKETCH_ROWS: usize = 4;
const SKETCH_WIDTH: usize = 1 << 20;
const DISTINCT_SAMPLE_SIZE: usize = 1 << 16;
const MIN_DISTINCT_MINIMIZERS: usize = 128;
const HIGH_QUANTILE_PER_MILLE: usize = 999;
const SKEW_THRESHOLD: f64 = 16.0;

pub(crate) struct DepthSkewReport {
    pub(crate) score: Option<f64>,
    pub(crate) skewed: bool,
    pub(crate) sampled_records: usize,
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
                    "{verdict} (p99.9/median minimizer count: {score:.2}; sampled reads: {})",
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

pub(crate) struct DepthSkewDetector {
    rng: StdRng,
    #[cfg(test)]
    sample_all: bool,
    sketch: CountMinSketch,
    distinct: DistinctSample,
    positions: Vec<u32>,
    sampled_records: usize,
}

impl DepthSkewDetector {
    pub(crate) fn new(seed: Option<u64>) -> Self {
        let rng = match seed {
            Some(seed) => StdRng::seed_from_u64(seed),
            None => StdRng::from_rng(&mut rand::rng()),
        };
        Self {
            rng,
            #[cfg(test)]
            sample_all: false,
            sketch: CountMinSketch::new(),
            distinct: DistinctSample::new(DISTINCT_SAMPLE_SIZE),
            positions: Vec::new(),
            sampled_records: 0,
        }
    }

    pub(crate) fn observe(&mut self, sequence: &[u8]) {
        #[cfg(test)]
        let selected = self.sample_all || self.rng.random_ratio(1, READ_SAMPLE_DENOMINATOR);
        #[cfg(not(test))]
        let selected = self.rng.random_ratio(1, READ_SAMPLE_DENOMINATOR);
        if !selected {
            return;
        }

        self.sampled_records += 1;
        for segment in sequence
            .split(|base| !matches!(base, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
        {
            if segment.len() < KMER_LENGTH + WINDOW_WIDTH - 1 {
                continue;
            }

            self.positions.clear();
            let minimizers = simd_minimizers::canonical_minimizers(KMER_LENGTH, WINDOW_WIDTH)
                .run(AsciiSeq(segment), &mut self.positions);
            for value in minimizers.values_u64() {
                self.sketch.increment(value);
                self.distinct.observe(value);
            }
        }
    }

    pub(crate) fn finish(self) -> DepthSkewReport {
        let sampled_records = self.sampled_records;
        let mut counts: Vec<u32> = self
            .distinct
            .into_keys()
            .map(|value| self.sketch.estimate(value))
            .collect();

        if counts.len() < MIN_DISTINCT_MINIMIZERS {
            return DepthSkewReport {
                score: None,
                skewed: false,
                sampled_records,
            };
        }

        counts.sort_unstable();
        let median = counts[(counts.len() - 1) / 2];
        let high_index = (counts.len() * HIGH_QUANTILE_PER_MILLE)
            .div_ceil(1000)
            .saturating_sub(1);
        let high = counts[high_index];
        let score = high as f64 / median.max(1) as f64;

        DepthSkewReport {
            score: Some(score),
            skewed: score >= SKEW_THRESHOLD,
            sampled_records,
        }
    }

    #[cfg(test)]
    fn sampling_all() -> Self {
        let mut detector = Self::new(Some(0));
        detector.sample_all = true;
        detector
    }
}

struct CountMinSketch {
    counters: Box<[u32]>,
}

impl CountMinSketch {
    fn new() -> Self {
        Self {
            counters: vec![0; SKETCH_ROWS * SKETCH_WIDTH].into_boxed_slice(),
        }
    }

    fn increment(&mut self, value: u64) {
        for row in 0..SKETCH_ROWS {
            let index = row * SKETCH_WIDTH + sketch_index(value, row);
            self.counters[index] = self.counters[index].saturating_add(1);
        }
    }

    fn estimate(&self, value: u64) -> u32 {
        (0..SKETCH_ROWS)
            .map(|row| self.counters[row * SKETCH_WIDTH + sketch_index(value, row)])
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

struct DistinctSample {
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

    fn into_keys(self) -> impl Iterator<Item = u64> {
        self.keys.into_iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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

    #[test]
    fn distinguishes_high_copy_sequence_from_even_coverage() {
        let chromosome = pseudo_random_dna(20_000, 1);
        let plasmid = pseudo_random_dna(2_000, 2);

        let mut even = DepthSkewDetector::sampling_all();
        for _ in 0..40 {
            even.observe(&chromosome);
        }
        let even = even.finish();

        let mut skewed = DepthSkewDetector::sampling_all();
        for _ in 0..40 {
            skewed.observe(&chromosome);
        }
        for _ in 0..800 {
            skewed.observe(&plasmid);
        }
        let skewed = skewed.finish();

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
            "Depth skew detected (p99.9/median minimizer count: 21.25; sampled reads: 87)"
        );
    }
}
