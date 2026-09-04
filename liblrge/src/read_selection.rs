use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use log::{debug, info, warn};
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};

use crate::depth_skew::{self, DepthProfile, DepthSkewDetector, DepthSkewReport, DistinctSample};
use crate::{io, Normalization, Result, DEFAULT_MAX_READ_BUFFER};

/// Bytes charged to a buffered read on top of its identifier and sequence, covering the two
/// heap allocations, the record itself, and a typical read identifier.
const BUFFERED_RECORD_OVERHEAD: u64 = 128;

pub(crate) struct ReadSelector {
    input: PathBuf,
    num_records: usize,
    total_bases: usize,
    seed: Option<u64>,
    normalization: Normalization,
    depth_skew: DepthSkewReport,
    /// The minimizers the detector sampled, held until the profiling pass scores them.
    minimizer_sample: Option<DistinctSample>,
    depth_profile: Option<DepthProfile>,
    threads: usize,
    max_read_buffer: u64,
}

pub(crate) struct SelectionResult {
    pub(crate) lengths: Vec<usize>,
    pub(crate) output_records: Vec<usize>,
    pub(crate) retained_records: usize,
    pub(crate) normalized: bool,
    /// Whether the selection fell back to the low-memory path.
    pub(crate) low_memory: bool,
    /// Peak bytes the normalized paths held to sample with, counting the reads or positions
    /// they buffered and nothing else. Zero when normalization did not run.
    pub(crate) buffer_bytes: u64,
}

impl ReadSelector {
    pub(crate) fn new<P: AsRef<Path>>(
        input: P,
        seed: Option<u64>,
        normalization: Normalization,
    ) -> Result<Self> {
        let input = input.as_ref().to_path_buf();
        let started = Instant::now();
        let mut total_bases = 0_usize;
        // Refusing normalization makes both the verdict and the sample irrelevant, so those runs
        // do no depth work in the counting pass at all. Forcing it still runs detection, for the
        // sample rather than the verdict: see [`ensure_depth_profile`].
        let mut detector = match normalization {
            Normalization::Auto | Normalization::Always => Some(DepthSkewDetector::new(seed)),
            Normalization::Never => None,
        };
        let num_records = io::count_records(&input, |sequence| {
            total_bases += sequence.len();
            if let Some(detector) = &mut detector {
                detector.observe(sequence);
            }
            Ok(())
        })?;
        let (depth_skew, minimizer_sample) = match detector {
            Some(detector) => {
                let detection = detector.finish();
                (detection.report, Some(detection.sample))
            }
            None => (DepthSkewReport::not_assessed(), None),
        };
        debug!("Counted {num_records} reads in {:.2?}", started.elapsed());

        Ok(Self {
            input,
            num_records,
            total_bases,
            seed,
            normalization,
            depth_skew,
            minimizer_sample,
            depth_profile: None,
            threads: 1,
            max_read_buffer: DEFAULT_MAX_READ_BUFFER,
        })
    }

    /// Set the number of threads the depth profiling pass may use.
    ///
    /// The pass only runs when normalization engages, and the profile it builds is the same
    /// whatever thread count builds it.
    pub(crate) fn threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    /// Set the cap on the bytes of selected reads normalization may hold in memory at once.
    ///
    /// A request projected to exceed the cap is served by a low-memory path that buffers read
    /// positions instead of read sequences and pays one extra pass over the input. Both paths
    /// select the same reads for a given seed.
    pub(crate) fn max_read_buffer(mut self, bytes: u64) -> Self {
        self.max_read_buffer = bytes;
        self
    }

    /// Bytes the buffered path would hold if every requested read were selected.
    fn projected_buffer_bytes(&self, capacity: usize) -> u64 {
        let mean_record_bases = match self.num_records {
            0 => 0,
            count => (self.total_bases / count) as u64,
        };
        (capacity as u64).saturating_mul(mean_record_bases + BUFFERED_RECORD_OVERHEAD)
    }

    pub(crate) fn num_records(&self) -> usize {
        self.num_records
    }

    pub(crate) fn depth_skew(&self) -> &DepthSkewReport {
        &self.depth_skew
    }

    pub(crate) fn normalization_message(&self, selection: &SelectionResult) -> Option<String> {
        if !selection.normalized {
            return None;
        }

        match self.normalization {
            Normalization::Auto => Some(format!(
                "{}; depth normalization retained {} of {} reads",
                self.depth_skew, selection.retained_records, self.num_records
            )),
            Normalization::Always => Some(format!(
                "Depth normalization forced; retained {} of {} reads",
                selection.retained_records, self.num_records
            )),
            Normalization::Never => None,
        }
    }

    pub(crate) fn ensure_supported_record_count(&self) -> Result<()> {
        if self.num_records > u32::MAX as usize {
            let msg = format!(
                "Number of reads in input file ({}) exceeds maximum allowed value ({})",
                self.num_records,
                u32::MAX
            );
            return Err(crate::error::LrgeError::TooManyReadsError(msg));
        }

        Ok(())
    }

    /// Select the reads each output asked for and write them.
    ///
    /// The work is a method of its own so that one timer covers every way out of it, the depth
    /// profile the normalized paths build included.
    pub(crate) fn write_selected(
        &mut self,
        outputs: &[(&Path, usize)],
        weights: Option<&[f64]>,
    ) -> Result<SelectionResult> {
        let started = Instant::now();
        let selection = self.select_and_write(outputs, weights)?;
        debug!(
            "Selected and wrote the reads in {:.2?}, the depth profile above included",
            started.elapsed()
        );
        Ok(selection)
    }

    fn select_and_write(
        &mut self,
        outputs: &[(&Path, usize)],
        weights: Option<&[f64]>,
    ) -> Result<SelectionResult> {
        if self.normalization_engaged() {
            assert!(
                weights.is_none(),
                "normalization computes its own read weights"
            );
            let selection = self.write_normalized(outputs)?;
            if selection.low_memory {
                info!(
                    "Selected reads by position to stay under the {} byte read-buffer cap; this read the input one extra time and held {} bytes",
                    self.max_read_buffer, selection.buffer_bytes
                );
                if selection.buffer_bytes > self.max_read_buffer {
                    warn!(
                        "Tracking this many selected reads takes {} bytes even without their sequences, over the {} byte cap; request fewer reads or raise --max-read-buffer",
                        selection.buffer_bytes, self.max_read_buffer
                    );
                }
            } else if selection.buffer_bytes > self.max_read_buffer {
                // The cap is applied to a projection from the mean read length, so a run whose
                // longer reads survive normalization can land here.
                warn!(
                    "Buffering the selected reads took {} bytes, over the {} byte cap; lower --max-read-buffer to force the low-memory path",
                    selection.buffer_bytes, self.max_read_buffer
                );
            } else {
                debug!(
                    "Buffered the selected reads in {} bytes, under the {} byte cap",
                    selection.buffer_bytes, self.max_read_buffer
                );
            }
            return Ok(selection);
        }

        let num_selected = outputs.iter().map(|(_, count)| count).sum();
        let indices =
            crate::sample_unique_indices(num_selected, self.num_records as u32, self.seed, weights);
        let mut selected_indices = match outputs.len() {
            2 => {
                let (first, second) = split_into_hashsets(indices, outputs[0].1);
                vec![first, second]
            }
            _ => partition_into_hashsets(indices, outputs.iter().map(|(_, count)| *count)),
        };

        let mut writers = Vec::with_capacity(outputs.len());
        for (path, _) in outputs {
            writers.push(File::create(path).map(BufWriter::new)?);
        }

        let mut lengths = vec![0; outputs.len()];
        let mut index = 0;
        io::iter_records(&self.input, |id, seq| {
            for (output_index, selected) in selected_indices.iter_mut().enumerate() {
                if selected.remove(&index) {
                    write_record(&mut writers[output_index], id, seq)?;
                    lengths[output_index] += seq.len();
                    break;
                }
            }

            index += 1;
            Ok(())
        })?;

        Ok(SelectionResult {
            lengths,
            output_records: outputs.iter().map(|(_, count)| *count).collect(),
            retained_records: self.num_records,
            normalized: false,
            low_memory: false,
            buffer_bytes: 0,
        })
    }

    /// Build the full depth profile unless a previous selection already did.
    ///
    /// This is the pass that counts every minimizer of every read. Holding it back until
    /// normalization is going to happen is what keeps it off an unskewed run.
    ///
    /// The profile normalizes against the median count of the minimizers detection sampled, which
    /// is why every mode that normalizes has to have run detection. Drawing the sample here
    /// instead, from every minimizer of every read, would draw it from a population sequencing
    /// error dominates: nearly every distinct minimizer of a long-read set is seen in one read
    /// only, so the median count over them is the sketch's noise floor rather than the input's
    /// coverage depth. Detection looks at about one read in a hundred, a population small enough
    /// that the minimizers the genome repeats are the majority of it.
    fn ensure_depth_profile(&mut self) -> Result<()> {
        if self.depth_profile.is_some() {
            return Ok(());
        }

        let sample = self
            .minimizer_sample
            .take()
            .expect("normalization only engages on a run that ran detection");
        self.depth_profile = Some(depth_skew::profile(&self.input, self.threads, sample)?);
        Ok(())
    }

    fn normalization_engaged(&self) -> bool {
        match self.normalization {
            Normalization::Auto => self.depth_skew.skewed,
            Normalization::Always => true,
            Normalization::Never => false,
        }
    }

    fn write_normalized(&mut self, outputs: &[(&Path, usize)]) -> Result<SelectionResult> {
        self.ensure_depth_profile()?;
        let capacity = outputs.iter().map(|(_, count)| count).sum();
        let projected = self.projected_buffer_bytes(capacity);
        debug!(
            "Buffering {capacity} selected reads projects to about {projected} bytes against a {} byte cap",
            self.max_read_buffer
        );

        // Lending the profile to the writer is what lets them take it as a plain argument, so
        // neither has to restate that normalization cannot run without one.
        let profile = self
            .depth_profile
            .as_ref()
            .expect("ensure_depth_profile built one");
        if projected > self.max_read_buffer {
            self.write_normalized_streamed(profile, outputs, capacity)
        } else {
            self.write_normalized_buffered(profile, outputs, capacity)
        }
    }

    /// Select and buffer the reads in one pass, then write them out.
    fn write_normalized_buffered(
        &self,
        profile: &DepthProfile,
        outputs: &[(&Path, usize)],
        capacity: usize,
    ) -> Result<SelectionResult> {
        let mut sampler = NormalizedSampler::new(capacity, self.seed);
        let mut reservoir: Vec<BufferedRecord> = Vec::with_capacity(capacity);
        let mut buffer_bytes = 0_u64;
        let mut peak_buffer_bytes = 0_u64;

        depth_skew::score_reads(&self.input, profile, self.threads, |input_index, read| {
            if let Some(slot) = sampler.offer(input_index, read.retention_probability) {
                // The scoring pass already owns a copy of the read, so buffering it takes the
                // buffers it holds rather than making another pair.
                let record = BufferedRecord {
                    input_index,
                    id: read.id,
                    sequence: read.sequence,
                };
                buffer_bytes += record.buffered_bytes();
                match slot {
                    ReservoirSlot::Append => reservoir.push(record),
                    ReservoirSlot::Replace(slot) => {
                        buffer_bytes -= reservoir[slot].buffered_bytes();
                        reservoir[slot] = record;
                    }
                }
                peak_buffer_bytes = peak_buffer_bytes.max(buffer_bytes);
            }
        })?;

        let order = sampler.finish();
        let mut slots: Vec<Option<BufferedRecord>> = reservoir.into_iter().map(Some).collect();
        let mut selected: Vec<BufferedRecord> = order
            .iter()
            .map(|slot| slots[*slot].take().expect("each slot is taken once"))
            .collect();

        let output_records = allocate_output_counts(
            outputs.iter().map(|(_, count)| *count).collect(),
            selected.len(),
        );
        let mut lengths = vec![0; outputs.len()];
        for (output_index, (records, (path, _))) in split_by_output(&mut selected, &output_records)
            .into_iter()
            .zip(outputs)
            .enumerate()
        {
            records.sort_unstable_by_key(|record| record.input_index);
            let mut writer = File::create(path).map(BufWriter::new)?;
            for record in records {
                write_record(&mut writer, &record.id, &record.sequence)?;
                lengths[output_index] += record.sequence.len();
            }
        }

        Ok(SelectionResult {
            lengths,
            output_records,
            retained_records: sampler.retained_records,
            normalized: true,
            low_memory: false,
            buffer_bytes: peak_buffer_bytes + vec_bytes::<usize>(capacity),
        })
    }

    /// Select read positions in one pass, then write the reads out in a second.
    ///
    /// This holds an input position per selected read rather than the read itself, which is what
    /// keeps a large request inside the read-buffer cap. It drives the same [`NormalizedSampler`]
    /// in the same order as the buffered path, so a seed picks the same reads either way.
    fn write_normalized_streamed(
        &self,
        profile: &DepthProfile,
        outputs: &[(&Path, usize)],
        capacity: usize,
    ) -> Result<SelectionResult> {
        let mut sampler = NormalizedSampler::new(capacity, self.seed);

        depth_skew::score_reads(&self.input, profile, self.threads, |input_index, read| {
            sampler.offer(input_index, read.retention_probability);
        })?;

        let order = sampler.finish();
        let mut selected = sampler.selected_indices(order);
        let output_records = allocate_output_counts(
            outputs.iter().map(|(_, count)| *count).collect(),
            selected.len(),
        );

        // Assignments are walked in input order during the write pass, matching the order the
        // buffered path writes each output in.
        let mut assignments: Vec<(usize, usize)> = Vec::with_capacity(selected.len());
        for (output_index, block) in split_by_output(&mut selected, &output_records)
            .into_iter()
            .enumerate()
        {
            assignments.extend(block.iter().map(|input_index| (*input_index, output_index)));
        }
        assignments.sort_unstable();

        let buffer_bytes = sampler.slot_bytes()
            + vec_bytes::<usize>(selected.capacity())
            + vec_bytes::<(usize, usize)>(assignments.capacity());

        let mut writers = Vec::with_capacity(outputs.len());
        for (path, _) in outputs {
            writers.push(File::create(path).map(BufWriter::new)?);
        }

        let mut lengths = vec![0; outputs.len()];
        let mut cursor = 0;
        let mut input_index = 0_usize;
        io::iter_records(&self.input, |id, sequence| {
            if let Some((_, output_index)) = assignments
                .get(cursor)
                .filter(|(index, _)| *index == input_index)
                .copied()
            {
                write_record(&mut writers[output_index], id, sequence)?;
                lengths[output_index] += sequence.len();
                cursor += 1;
            }
            input_index += 1;
            Ok(())
        })?;

        Ok(SelectionResult {
            lengths,
            output_records,
            retained_records: sampler.retained_records,
            normalized: true,
            low_memory: true,
            buffer_bytes,
        })
    }
}

/// The reservoir decisions shared by both normalized selection paths.
///
/// Both paths offer the same records in the same order with the same retention probabilities, so
/// they draw the same random numbers and end with the same reservoir. That is what lets a seed
/// mean the same thing whichever path a request lands on.
struct NormalizedSampler {
    rng: StdRng,
    capacity: usize,
    retained_records: usize,
    slots: Vec<usize>,
}

impl NormalizedSampler {
    fn new(capacity: usize, seed: Option<u64>) -> Self {
        let rng = match seed {
            Some(seed) => StdRng::seed_from_u64(seed),
            None => StdRng::from_rng(&mut rand::rng()),
        };
        Self {
            rng,
            capacity,
            retained_records: 0,
            slots: Vec::with_capacity(capacity),
        }
    }

    /// Offer the record at `input_index`, which survives normalization with `probability`.
    ///
    /// Returns where the record belongs in the reservoir, or `None` if it was dropped.
    fn offer(&mut self, input_index: usize, probability: f64) -> Option<ReservoirSlot> {
        if probability < 1.0 && !self.rng.random_bool(probability) {
            return None;
        }

        self.retained_records += 1;
        if self.capacity == 0 {
            return None;
        }

        if self.slots.len() < self.capacity {
            self.slots.push(input_index);
            return Some(ReservoirSlot::Append);
        }

        let slot = self.rng.random_range(0..self.retained_records);
        if slot < self.capacity {
            self.slots[slot] = input_index;
            Some(ReservoirSlot::Replace(slot))
        } else {
            None
        }
    }

    /// Shuffle the reservoir, returning the slots in their final order.
    fn finish(&mut self) -> Vec<usize> {
        let mut order: Vec<usize> = (0..self.slots.len()).collect();
        order.shuffle(&mut self.rng);
        order
    }

    /// The input positions of the selected reads, in the order `finish` settled on.
    fn selected_indices(&self, order: Vec<usize>) -> Vec<usize> {
        order.into_iter().map(|slot| self.slots[slot]).collect()
    }

    fn slot_bytes(&self) -> u64 {
        vec_bytes::<usize>(self.slots.capacity())
    }
}

/// Where an offered record belongs in the reservoir.
enum ReservoirSlot {
    /// After the records already held, which are still fewer than the reservoir takes.
    Append,
    /// In place of the record already at this position.
    Replace(usize),
}

fn vec_bytes<T>(capacity: usize) -> u64 {
    (capacity * std::mem::size_of::<T>()) as u64
}

/// Split the selected reads into one block per output, in the order the outputs were requested.
fn split_by_output<'a, T>(selected: &'a mut [T], counts: &[usize]) -> Vec<&'a mut [T]> {
    let mut remaining = selected;
    let mut blocks = Vec::with_capacity(counts.len());
    for count in counts {
        let (block, rest) = remaining.split_at_mut(*count);
        blocks.push(block);
        remaining = rest;
    }
    blocks
}

fn write_record<W: Write>(writer: &mut W, id: &[u8], sequence: &[u8]) -> std::io::Result<()> {
    writer.write_all(b">")?;
    writer.write_all(id)?;
    writer.write_all(b"\n")?;
    writer.write_all(sequence)?;
    writer.write_all(b"\n")
}

struct BufferedRecord {
    input_index: usize,
    id: Vec<u8>,
    sequence: Vec<u8>,
}

impl BufferedRecord {
    fn buffered_bytes(&self) -> u64 {
        (self.id.len() + self.sequence.len()) as u64 + BUFFERED_RECORD_OVERHEAD
    }
}

fn allocate_output_counts(mut requested: Vec<usize>, selected: usize) -> Vec<usize> {
    let total_requested = requested.iter().sum::<usize>();
    let mut shortage = total_requested.saturating_sub(selected);
    if requested.len() == 2
        && requested.iter().all(|count| *count > 0)
        && shortage > 0
        && selected >= 2
        && selected <= requested[1]
    {
        let first = ((selected as u128 * requested[0] as u128) / total_requested as u128)
            .clamp(1, (selected - 1) as u128) as usize;
        return vec![first, selected - first];
    }

    let nonempty_outputs = requested.iter().filter(|count| **count > 0).count();
    let keep_one_per_output = selected >= nonempty_outputs;
    for count in &mut requested {
        let minimum = usize::from(keep_one_per_output && *count > 0);
        let reduction = shortage.min(*count - minimum);
        *count -= reduction;
        shortage -= reduction;
    }
    requested
}

fn partition_into_hashsets<T: std::hash::Hash + Eq>(
    mut original: Vec<T>,
    sizes: impl IntoIterator<Item = usize>,
) -> Vec<HashSet<T>> {
    sizes
        .into_iter()
        .map(|size| {
            let mut set = HashSet::with_capacity(size.min(original.len()));
            for _ in 0..size.min(original.len()) {
                if let Some(element) = original.pop() {
                    set.insert(element);
                }
            }
            set
        })
        .collect()
}

pub(crate) fn split_into_hashsets<T: std::hash::Hash + Eq>(
    original: Vec<T>,
    size_first: usize,
) -> (HashSet<T>, HashSet<T>) {
    let mut sets = partition_into_hashsets(original, [size_first, usize::MAX]);
    let second = sets.pop().expect("partition has two sets");
    let first = sets.pop().expect("partition has two sets");
    (first, second)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Normalization;
    use std::io::Write;

    /// One step of the generator the test inputs are built from.
    fn next_draw(state: &mut u64) -> u64 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        *state >> 32
    }

    fn pseudo_random_dna(len: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        (0..len)
            .map(|_| b"ACGT"[(next_draw(&mut state) & 3) as usize])
            .collect()
    }

    /// What running an input through one normalization mode came to.
    struct Normalized {
        selector: ReadSelector,
        selection: SelectionResult,
        /// The bytes of the reads it wrote.
        reads: Vec<u8>,
    }

    /// Select twenty reads from `input` under `mode`, on `threads` threads.
    ///
    /// The tests that compare one normalizing mode against the other compare the whole answer:
    /// what the run measured, how many reads it kept, and which reads it wrote.
    fn normalized_selection(input: &Path, mode: Normalization, threads: usize) -> Normalized {
        let output = tempfile::NamedTempFile::new().unwrap();
        let mut selector = ReadSelector::new(input, Some(42), mode)
            .unwrap()
            .threads(threads);
        let selection = selector
            .write_selected(&[(output.path(), 20)], None)
            .unwrap();
        let reads = std::fs::read(output.path()).unwrap();
        Normalized {
            selector,
            selection,
            reads,
        }
    }

    #[test]
    fn read_selector_writes_seeded_groups() {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            input,
            ">read0\nAAA\n>read1\nCCC\n>read2\nGGG\n>read3\nTTT\n>read4\nACG\n>read5\nTAC"
        )
        .unwrap();
        let tempdir = tempfile::tempdir().unwrap();
        let target = tempdir.path().join("target.fa");
        let query = tempdir.path().join("query.fa");

        let mut selector = ReadSelector::new(input.path(), Some(42), Normalization::Auto).unwrap();
        assert_eq!(selector.num_records(), 6);

        let selection = selector
            .write_selected(&[(target.as_path(), 2), (query.as_path(), 1)], None)
            .unwrap();

        assert_eq!(selection.lengths, vec![6, 3]);
        assert_eq!(
            std::fs::read_to_string(target).unwrap(),
            ">read1\nCCC\n>read2\nGGG\n"
        );
        assert_eq!(std::fs::read_to_string(query).unwrap(), ">read0\nAAA\n");

        let reads = tempdir.path().join("reads.fa");
        let selection = selector
            .write_selected(&[(reads.as_path(), 3)], None)
            .unwrap();

        assert_eq!(selection.lengths, vec![9]);
        assert_eq!(
            std::fs::read_to_string(reads).unwrap(),
            ">read0\nAAA\n>read1\nCCC\n>read2\nGGG\n"
        );
    }

    #[test]
    fn auto_preserves_legacy_output_when_depth_is_not_skewed() {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        for index in 0..100 {
            writeln!(
                input,
                ">read{index}\n{}",
                String::from_utf8(pseudo_random_dna(250, index + 1)).unwrap()
            )
            .unwrap();
        }
        let tempdir = tempfile::tempdir().unwrap();
        let auto_path = tempdir.path().join("auto.fa");
        let never_path = tempdir.path().join("never.fa");

        let mut auto = ReadSelector::new(input.path(), Some(42), Normalization::Auto).unwrap();
        auto.write_selected(&[(auto_path.as_path(), 20)], None)
            .unwrap();
        let mut never = ReadSelector::new(input.path(), Some(42), Normalization::Never).unwrap();
        never
            .write_selected(&[(never_path.as_path(), 20)], None)
            .unwrap();

        assert!(!auto.depth_skew().skewed);
        assert_eq!(
            std::fs::read(auto_path).unwrap(),
            std::fs::read(never_path).unwrap()
        );
    }

    #[test]
    fn always_normalizes_when_auto_does_not_engage() {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        for index in 0..100 {
            writeln!(
                input,
                ">read{index}\n{}",
                String::from_utf8(pseudo_random_dna(250, index + 1)).unwrap()
            )
            .unwrap();
        }
        let tempdir = tempfile::tempdir().unwrap();
        let auto_path = tempdir.path().join("auto.fa");
        let always_path = tempdir.path().join("always.fa");

        let mut auto = ReadSelector::new(input.path(), Some(42), Normalization::Auto).unwrap();
        let auto_selection = auto
            .write_selected(&[(auto_path.as_path(), 20)], None)
            .unwrap();
        let mut always = ReadSelector::new(input.path(), Some(42), Normalization::Always).unwrap();
        let always_selection = always
            .write_selected(&[(always_path.as_path(), 20)], None)
            .unwrap();

        assert!(!auto_selection.normalized);
        assert!(auto.normalization_message(&auto_selection).is_none());
        assert!(always_selection.normalized);
        assert!(always
            .normalization_message(&always_selection)
            .unwrap()
            .contains("forced"));
        assert!(
            !always.depth_skew().skewed,
            "the input was called skewed, so this does not test forcing"
        );
    }

    /// An input whose distinct minimizers are mostly ones that occur in a single read.
    ///
    /// This is the shape a real long-read set has and the small inputs above do not. Reads carry
    /// substitutions, so most of what the whole input has ever seen occurs in one read only,
    /// while the minimizers of the sequence the reads came from recur once per read that covers
    /// them. A high-copy element gives the input a skew for detection to find.
    ///
    /// It holds more reads than the detection floor times the read denominator, so detection's
    /// own draw clears the floor. Below that the top-up hands detection every read, the sample is
    /// the whole input's either way, and where it came from cannot matter.
    fn error_prone_skewed_input() -> tempfile::NamedTempFile {
        const READ_LENGTH: usize = 40;
        const RECORDS: usize = 60_000;
        // Fourteen reads in sixty come from the element. It is five hundred times shorter than
        // the chromosome, so those reads cover it two orders of magnitude more deeply, which is
        // the skew detection has to find.
        const ELEMENT_SHARE: usize = 14;
        // One base in this many is redrawn. Enough that most of the input's distinct minimizers
        // are singletons, few enough that a read still carries the sequence it came from.
        const ERROR_DENOMINATOR: u64 = 200;

        let chromosome = pseudo_random_dna(50_000, 21);
        let element = pseudo_random_dna(100, 22);
        let mut state = 99_u64;
        let mut input = tempfile::NamedTempFile::new().unwrap();
        let mut read = Vec::with_capacity(READ_LENGTH);
        for index in 0..RECORDS {
            let source = if index % 60 < ELEMENT_SHARE {
                &element
            } else {
                &chromosome
            };
            let start = next_draw(&mut state) as usize % source.len();
            read.clear();
            for offset in 0..READ_LENGTH {
                let mut base = source[(start + offset) % source.len()];
                if next_draw(&mut state).is_multiple_of(ERROR_DENOMINATOR) {
                    base = b"ACGT"[(next_draw(&mut state) & 3) as usize];
                }
                read.push(base);
            }
            writeln!(input, ">read{index}").unwrap();
            input.write_all(&read).unwrap();
            input.write_all(b"\n").unwrap();
        }
        input.flush().unwrap();
        input
    }

    /// A forced normalization has to reach the answer the detector's own run reaches.
    ///
    /// Both modes normalize against the median count of a bottom-k sample of minimizers, and what
    /// that median measures turns on which reads the sample was drawn from. Drawn from every read
    /// of an error-prone input it lands among the singletons, which is the sketch's noise floor
    /// rather than the genome's depth, and normalization then throws away nearly every read.
    /// Drawn from the hundredth of the reads detection looks at, it lands on the coverage depth.
    /// Forcing normalization must not change which of the two the run measures.
    #[test]
    fn forcing_normalization_measures_the_same_depth_as_detecting_it() {
        let input = error_prone_skewed_input();

        // Both passes over the input parallelise, and neither depends on how the work is divided,
        // so the threads only buy the suite back some of the time this input costs.
        let auto = normalized_selection(input.path(), Normalization::Auto, 4);
        let always = normalized_selection(input.path(), Normalization::Always, 4);
        let detection = auto.selector.depth_skew();

        assert!(
            detection.skewed,
            "the input was not called skewed: {detection}"
        );
        assert!(
            detection.sampled_records < auto.selector.num_records() / 10,
            "detection sampled {} of {} reads, so the sample is the whole input either way",
            detection.sampled_records,
            auto.selector.num_records()
        );
        assert!(auto.selection.normalized);
        assert_eq!(
            auto.selection.retained_records, always.selection.retained_records,
            "forcing normalization retained a different number of reads than detecting skew did"
        );
        assert_eq!(auto.reads, always.reads);
    }

    /// An input holding fewer reads than the detection floor has every one of them sampled.
    ///
    /// The floor is there because a hundredth of a small input is too few reads to place a
    /// percentile on, and it means the sample a small input's profile is built on is drawn from
    /// the whole input rather than a hundredth of it. That is the one case where a forced
    /// normalization measured the same depth before this pass ran detection as it does now, so
    /// the two modes agree here by a different route than they do above.
    #[test]
    fn an_input_below_the_detection_floor_is_sampled_whole() {
        let input = skewed_input();

        let auto = normalized_selection(input.path(), Normalization::Auto, 1);
        let always = normalized_selection(input.path(), Normalization::Always, 1);
        let records = always.selector.num_records();
        let detection = always.selector.depth_skew();

        assert!(
            records < crate::depth_skew::DETECTION_READ_FLOOR,
            "the input holds {records} reads, which is not below the floor"
        );
        assert_eq!(
            detection.sampled_records, records,
            "detection left some of a below-floor input unsampled"
        );
        assert!(
            detection.skewed,
            "the input was not called skewed: {detection}"
        );
        assert!(auto.selection.normalized);
        assert_eq!(
            auto.selection.retained_records,
            always.selection.retained_records
        );
        assert_eq!(auto.reads, always.reads);
    }

    #[test]
    fn uniform_weights_write_the_same_reads_as_unweighted_selection() {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            input,
            ">read0\nAAA\n>read1\nCCC\n>read2\nGGG\n>read3\nTTT\n>read4\nACG\n>read5\nTAC"
        )
        .unwrap();
        let tempdir = tempfile::tempdir().unwrap();
        let unweighted = tempdir.path().join("unweighted.fa");
        let weighted = tempdir.path().join("weighted.fa");
        let mut selector = ReadSelector::new(input.path(), Some(42), Normalization::Auto).unwrap();

        selector
            .write_selected(&[(unweighted.as_path(), 3)], None)
            .unwrap();
        selector
            .write_selected(&[(weighted.as_path(), 3)], Some(&[1.0; 6]))
            .unwrap();

        assert_eq!(
            std::fs::read(unweighted).unwrap(),
            std::fs::read(weighted).unwrap()
        );
    }

    #[test]
    fn zero_weight_reads_are_not_written() {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        writeln!(input, ">read0\nAAA\n>read1\nCCC\n>read2\nGGG").unwrap();
        let tempdir = tempfile::tempdir().unwrap();
        let selected = tempdir.path().join("selected.fa");
        let mut selector = ReadSelector::new(input.path(), Some(7), Normalization::Auto).unwrap();
        let weights = [1.0, 0.0, 1.0];

        selector
            .write_selected(&[(selected.as_path(), 2)], Some(&weights))
            .unwrap();

        let expected = b">read0\nAAA\n>read2\nGGG\n";
        assert_eq!(std::fs::read(selected).unwrap(), expected);
    }

    /// The full depth profile is the expensive pass. An unskewed run must not build one.
    #[test]
    fn auto_does_not_profile_an_unskewed_input() {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        for index in 0..100 {
            writeln!(
                input,
                ">read{index}\n{}",
                String::from_utf8(pseudo_random_dna(250, index + 1)).unwrap()
            )
            .unwrap();
        }
        let tempdir = tempfile::tempdir().unwrap();
        let selected = tempdir.path().join("selected.fa");

        let mut auto = ReadSelector::new(input.path(), Some(42), Normalization::Auto).unwrap();
        auto.write_selected(&[(selected.as_path(), 20)], None)
            .unwrap();

        assert!(!auto.depth_skew().skewed);
        assert!(
            auto.depth_profile.is_none(),
            "an unskewed run built the full depth profile"
        );
    }

    /// Reads are scored in parallel and sampled sequentially, so a thread count has to be free to
    /// change without moving what a seed selects. Both normalized paths score every read, so both
    /// are checked.
    #[test]
    fn a_seed_selects_the_same_reads_at_any_thread_count() {
        let input = skewed_input();
        let tempdir = tempfile::tempdir().unwrap();

        for (path_prefix, cap) in [("buffered", u64::MAX), ("streamed", 0)] {
            let one_thread = tempdir.path().join(format!("{path_prefix}_one.fa"));
            let many_threads = tempdir.path().join(format!("{path_prefix}_many.fa"));

            let mut single = ReadSelector::new(input.path(), Some(42), Normalization::Always)
                .unwrap()
                .threads(1)
                .max_read_buffer(cap);
            let single_selection = single
                .write_selected(&[(one_thread.as_path(), 20)], None)
                .unwrap();
            let mut many = ReadSelector::new(input.path(), Some(42), Normalization::Always)
                .unwrap()
                .threads(4)
                .max_read_buffer(cap);
            let many_selection = many
                .write_selected(&[(many_threads.as_path(), 20)], None)
                .unwrap();

            assert!(single_selection.normalized);
            assert_eq!(single_selection.low_memory, cap == 0);
            assert_eq!(many_selection.low_memory, cap == 0);
            assert_eq!(
                single_selection.retained_records, many_selection.retained_records,
                "{path_prefix} path"
            );
            assert_eq!(
                std::fs::read(one_thread).unwrap(),
                std::fs::read(many_threads).unwrap(),
                "{path_prefix} path"
            );
        }
    }

    /// One read in a hundred of a thousand is about ten, too few to place a percentile, so the
    /// selector's detection pass tops the sample up to the floor.
    #[test]
    fn read_selector_samples_records_for_depth_skew_detection() {
        let records = 1_000;
        let mut input = tempfile::NamedTempFile::new().unwrap();
        for index in 0..records {
            writeln!(input, ">read{index}\nACGTACGTACGTACGTACGTACGTACGT").unwrap();
        }

        let selector = ReadSelector::new(input.path(), Some(42), Normalization::Auto).unwrap();
        let sampled = selector.depth_skew().sampled_records;

        assert_eq!(sampled, crate::depth_skew::DETECTION_READ_FLOOR);
    }

    #[test]
    fn normalization_reduces_reads_from_a_high_copy_element() {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        for index in 0..20 {
            writeln!(
                input,
                ">chromosome{index}\n{}",
                String::from_utf8(pseudo_random_dna(250, index + 1)).unwrap()
            )
            .unwrap();
        }
        let element = String::from_utf8(pseudo_random_dna(250, 10_000)).unwrap();
        for index in 0..400 {
            writeln!(input, ">element{index}\n{element}").unwrap();
        }

        let tempdir = tempfile::tempdir().unwrap();
        let legacy_path = tempdir.path().join("legacy.fa");
        let normalized_path = tempdir.path().join("normalized.fa");

        let mut legacy = ReadSelector::new(input.path(), Some(42), Normalization::Never).unwrap();
        legacy
            .write_selected(&[(legacy_path.as_path(), 20)], None)
            .unwrap();
        let mut normalized =
            ReadSelector::new(input.path(), Some(42), Normalization::Always).unwrap();
        normalized
            .write_selected(&[(normalized_path.as_path(), 20)], None)
            .unwrap();

        let legacy = std::fs::read_to_string(legacy_path).unwrap();
        let normalized = std::fs::read_to_string(normalized_path).unwrap();
        let legacy_chromosome_reads = legacy.matches(">chromosome").count();
        let normalized_chromosome_reads = normalized.matches(">chromosome").count();

        assert!(
            legacy_chromosome_reads <= 5,
            "legacy selected {legacy_chromosome_reads} chromosome reads"
        );
        assert!(
            normalized_chromosome_reads >= 15,
            "normalized selected {normalized_chromosome_reads} chromosome reads"
        );
    }

    #[test]
    fn normalized_pool_can_be_smaller_than_the_request() {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        for index in 0..20 {
            writeln!(
                input,
                ">chromosome{index}\n{}",
                String::from_utf8(pseudo_random_dna(250, index + 1)).unwrap()
            )
            .unwrap();
        }
        let element = String::from_utf8(pseudo_random_dna(250, 10_000)).unwrap();
        for index in 0..400 {
            writeln!(input, ">element{index}\n{element}").unwrap();
        }

        let tempdir = tempfile::tempdir().unwrap();
        let selected_path = tempdir.path().join("selected.fa");
        let mut selector =
            ReadSelector::new(input.path(), Some(42), Normalization::Always).unwrap();

        let result = selector
            .write_selected(&[(selected_path.as_path(), 100)], None)
            .unwrap();

        assert!(result.output_records[0] < 100);
        assert_eq!(result.output_records[0], result.retained_records);
        assert_eq!(
            std::fs::read_to_string(selected_path)
                .unwrap()
                .matches('>')
                .count(),
            result.output_records[0]
        );
    }

    fn skewed_input() -> tempfile::NamedTempFile {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        for index in 0..40 {
            writeln!(
                input,
                ">chromosome{index}\n{}",
                String::from_utf8(pseudo_random_dna(600, index + 1)).unwrap()
            )
            .unwrap();
        }
        let element = String::from_utf8(pseudo_random_dna(600, 10_000)).unwrap();
        for index in 0..400 {
            writeln!(input, ">element{index}\n{element}").unwrap();
        }
        input
    }

    #[test]
    fn both_paths_select_the_same_reads_for_a_given_seed() {
        let input = skewed_input();
        let tempdir = tempfile::tempdir().unwrap();
        let buffered_target = tempdir.path().join("buffered_target.fa");
        let buffered_query = tempdir.path().join("buffered_query.fa");
        let streamed_target = tempdir.path().join("streamed_target.fa");
        let streamed_query = tempdir.path().join("streamed_query.fa");

        let mut buffered = ReadSelector::new(input.path(), Some(42), Normalization::Always)
            .unwrap()
            .max_read_buffer(u64::MAX);
        let buffered_selection = buffered
            .write_selected(
                &[
                    (buffered_target.as_path(), 20),
                    (buffered_query.as_path(), 10),
                ],
                None,
            )
            .unwrap();

        let mut streamed = ReadSelector::new(input.path(), Some(42), Normalization::Always)
            .unwrap()
            .max_read_buffer(0);
        let streamed_selection = streamed
            .write_selected(
                &[
                    (streamed_target.as_path(), 20),
                    (streamed_query.as_path(), 10),
                ],
                None,
            )
            .unwrap();

        assert!(!buffered_selection.low_memory);
        assert!(streamed_selection.low_memory);
        assert_eq!(
            buffered_selection.retained_records,
            streamed_selection.retained_records
        );
        assert_eq!(
            buffered_selection.output_records,
            streamed_selection.output_records
        );
        assert_eq!(buffered_selection.lengths, streamed_selection.lengths);
        assert_eq!(
            std::fs::read(buffered_target).unwrap(),
            std::fs::read(streamed_target).unwrap()
        );
        assert_eq!(
            std::fs::read(buffered_query).unwrap(),
            std::fs::read(streamed_query).unwrap()
        );
    }

    #[test]
    fn both_paths_agree_when_the_pool_is_smaller_than_the_request() {
        let input = skewed_input();
        let tempdir = tempfile::tempdir().unwrap();
        let buffered_path = tempdir.path().join("buffered.fa");
        let streamed_path = tempdir.path().join("streamed.fa");

        let mut buffered = ReadSelector::new(input.path(), Some(7), Normalization::Always)
            .unwrap()
            .max_read_buffer(u64::MAX);
        let buffered_selection = buffered
            .write_selected(&[(buffered_path.as_path(), 400)], None)
            .unwrap();
        let mut streamed = ReadSelector::new(input.path(), Some(7), Normalization::Always)
            .unwrap()
            .max_read_buffer(0);
        let streamed_selection = streamed
            .write_selected(&[(streamed_path.as_path(), 400)], None)
            .unwrap();

        assert!(buffered_selection.output_records[0] < 400);
        assert_eq!(
            buffered_selection.output_records,
            streamed_selection.output_records
        );
        assert_eq!(
            std::fs::read(buffered_path).unwrap(),
            std::fs::read(streamed_path).unwrap()
        );
    }

    #[test]
    fn the_low_memory_path_holds_the_read_buffer_within_the_cap() {
        let input = skewed_input();
        let tempdir = tempfile::tempdir().unwrap();
        let buffered_path = tempdir.path().join("buffered.fa");
        let streamed_path = tempdir.path().join("streamed.fa");
        // Enough for 30 read indices, nowhere near enough for 30 six-hundred-base reads.
        let cap = 8 * 1024;

        let mut buffered = ReadSelector::new(input.path(), Some(42), Normalization::Always)
            .unwrap()
            .max_read_buffer(u64::MAX);
        let buffered_selection = buffered
            .write_selected(&[(buffered_path.as_path(), 30)], None)
            .unwrap();
        let mut streamed = ReadSelector::new(input.path(), Some(42), Normalization::Always)
            .unwrap()
            .max_read_buffer(cap);
        let streamed_selection = streamed
            .write_selected(&[(streamed_path.as_path(), 30)], None)
            .unwrap();

        assert!(buffered_selection.buffer_bytes > cap);
        assert!(streamed_selection.low_memory);
        assert!(
            streamed_selection.buffer_bytes <= cap,
            "low-memory path held {} bytes against a {cap} byte cap",
            streamed_selection.buffer_bytes
        );
    }

    #[test]
    fn a_request_that_fits_the_cap_stays_on_the_buffered_path() {
        let input = skewed_input();
        let tempdir = tempfile::tempdir().unwrap();
        let selected = tempdir.path().join("selected.fa");

        let mut selector = ReadSelector::new(input.path(), Some(42), Normalization::Always)
            .unwrap()
            .max_read_buffer(crate::DEFAULT_MAX_READ_BUFFER);
        let selection = selector
            .write_selected(&[(selected.as_path(), 30)], None)
            .unwrap();

        assert!(!selection.low_memory);
    }

    #[test]
    fn the_cap_does_not_change_unnormalized_selection() {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        for index in 0..100 {
            writeln!(
                input,
                ">read{index}\n{}",
                String::from_utf8(pseudo_random_dna(250, index + 1)).unwrap()
            )
            .unwrap();
        }
        let tempdir = tempfile::tempdir().unwrap();
        let generous = tempdir.path().join("generous.fa");
        let stingy = tempdir.path().join("stingy.fa");

        let mut selector = ReadSelector::new(input.path(), Some(42), Normalization::Never)
            .unwrap()
            .max_read_buffer(u64::MAX);
        selector
            .write_selected(&[(generous.as_path(), 20)], None)
            .unwrap();
        let mut selector = ReadSelector::new(input.path(), Some(42), Normalization::Never)
            .unwrap()
            .max_read_buffer(0);
        let selection = selector
            .write_selected(&[(stingy.as_path(), 20)], None)
            .unwrap();

        assert!(!selection.low_memory);
        assert_eq!(
            std::fs::read(generous).unwrap(),
            std::fs::read(stingy).unwrap()
        );
    }

    #[test]
    fn undersized_two_set_pool_preserves_the_query_when_possible() {
        assert_eq!(allocate_output_counts(vec![800, 200], 500), vec![300, 200]);
    }

    #[test]
    fn pool_smaller_than_the_query_scales_both_sets() {
        assert_eq!(allocate_output_counts(vec![800, 200], 100), vec![80, 20]);
    }
}
