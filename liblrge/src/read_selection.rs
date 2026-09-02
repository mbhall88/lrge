use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};

use crate::depth_skew::{DepthProfile, DepthSkewDetector, DepthSkewReport};
use crate::{io, Normalization, Result};

pub(crate) struct ReadSelector {
    input: PathBuf,
    num_records: usize,
    seed: Option<u64>,
    normalization: Normalization,
    depth_skew: DepthSkewReport,
    depth_profile: Option<DepthProfile>,
}

pub(crate) struct SelectionResult {
    pub(crate) lengths: Vec<usize>,
    pub(crate) output_records: Vec<usize>,
    pub(crate) retained_records: usize,
    pub(crate) normalized: bool,
}

impl ReadSelector {
    pub(crate) fn new<P: AsRef<Path>>(
        input: P,
        seed: Option<u64>,
        normalization: Normalization,
    ) -> Result<Self> {
        let input = input.as_ref().to_path_buf();
        let (num_records, depth_skew, depth_profile) = if normalization == Normalization::Never {
            let num_records = io::count_records(&input, |_| Ok(()))?;
            (
                num_records,
                DepthSkewReport {
                    score: None,
                    skewed: false,
                    sampled_records: 0,
                },
                None,
            )
        } else {
            let mut detector = match normalization {
                Normalization::Auto => DepthSkewDetector::detecting(seed),
                Normalization::Always => DepthSkewDetector::profiling(seed),
                Normalization::Never => unreachable!("handled above"),
            };
            let num_records = io::count_records(&input, |sequence| {
                detector.observe(sequence);
                Ok(())
            })?;
            let profile = detector.finish();
            let report = profile.report.clone();
            (num_records, report, Some(profile))
        };

        Ok(Self {
            input,
            num_records,
            seed,
            normalization,
            depth_skew,
            depth_profile,
        })
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

    pub(crate) fn write_selected(
        &mut self,
        outputs: &[(&Path, usize)],
        weights: Option<&[f64]>,
    ) -> Result<SelectionResult> {
        if self.normalization_engaged() {
            assert!(
                weights.is_none(),
                "normalization computes its own read weights"
            );
            return self.write_normalized(outputs);
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
                    let writer = &mut writers[output_index];
                    writer.write_all(b">")?;
                    writer.write_all(id)?;
                    writer.write_all(b"\n")?;
                    writer.write_all(seq)?;
                    writer.write_all(b"\n")?;
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
        })
    }

    fn normalization_engaged(&self) -> bool {
        match self.normalization {
            Normalization::Auto => self.depth_skew.skewed,
            Normalization::Always => true,
            Normalization::Never => false,
        }
    }

    fn write_normalized(&mut self, outputs: &[(&Path, usize)]) -> Result<SelectionResult> {
        let capacity = outputs.iter().map(|(_, count)| count).sum();
        let profile = self
            .depth_profile
            .as_mut()
            .expect("normalization requires a depth profile");
        let mut rng = match self.seed {
            Some(seed) => StdRng::seed_from_u64(seed),
            None => StdRng::from_rng(&mut rand::rng()),
        };
        let mut reservoir = Vec::with_capacity(capacity);
        let mut retained_records = 0_usize;
        let mut input_index = 0_usize;

        io::iter_records(&self.input, |id, sequence| {
            let probability = profile.retention_probability(sequence);
            if probability < 1.0 && !rng.random_bool(probability) {
                input_index += 1;
                return Ok(());
            }

            retained_records += 1;
            if capacity > 0 {
                let record = BufferedRecord {
                    input_index,
                    id: id.to_vec(),
                    sequence: sequence.to_vec(),
                };
                if reservoir.len() < capacity {
                    reservoir.push(record);
                } else {
                    let slot = rng.random_range(0..retained_records);
                    if slot < capacity {
                        reservoir[slot] = record;
                    }
                }
            }
            input_index += 1;
            Ok(())
        })?;

        reservoir.shuffle(&mut rng);
        let output_records = allocate_output_counts(
            outputs.iter().map(|(_, count)| *count).collect(),
            reservoir.len(),
        );
        let mut lengths = vec![0; outputs.len()];
        let mut offset = 0;
        for (output_index, ((path, _), count)) in outputs
            .iter()
            .zip(output_records.iter().copied())
            .enumerate()
        {
            let end = offset + count;
            let records = &mut reservoir[offset..end];
            records.sort_unstable_by_key(|record| record.input_index);
            let mut writer = File::create(path).map(BufWriter::new)?;
            for record in records {
                writer.write_all(b">")?;
                writer.write_all(&record.id)?;
                writer.write_all(b"\n")?;
                writer.write_all(&record.sequence)?;
                writer.write_all(b"\n")?;
                lengths[output_index] += record.sequence.len();
            }
            offset = end;
        }

        Ok(SelectionResult {
            lengths,
            output_records,
            retained_records,
            normalized: true,
        })
    }
}

struct BufferedRecord {
    input_index: usize,
    id: Vec<u8>,
    sequence: Vec<u8>,
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
        assert_eq!(always.depth_skew().sampled_records, 0);
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

    #[test]
    fn read_selector_samples_records_for_depth_skew_detection() {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        for index in 0..1_000 {
            writeln!(input, ">read{index}\nACGTACGTACGTACGTACGTACGTACGT").unwrap();
        }

        let selector = ReadSelector::new(input.path(), Some(42), Normalization::Auto).unwrap();
        let sampled = selector.depth_skew().sampled_records;

        assert!((3..=20).contains(&sampled), "sampled {sampled} reads");
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

    #[test]
    fn undersized_two_set_pool_preserves_the_query_when_possible() {
        assert_eq!(allocate_output_counts(vec![800, 200], 500), vec![300, 200]);
    }

    #[test]
    fn pool_smaller_than_the_query_scales_both_sets() {
        assert_eq!(allocate_output_counts(vec![800, 200], 100), vec![80, 20]);
    }
}
