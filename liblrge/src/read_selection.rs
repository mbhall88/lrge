use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::{io, Result};

pub(crate) struct ReadSelector {
    input: PathBuf,
    num_records: usize,
    seed: Option<u64>,
}

impl ReadSelector {
    pub(crate) fn new<P: AsRef<Path>>(input: P, seed: Option<u64>) -> Result<Self> {
        let input = input.as_ref().to_path_buf();
        let num_records = io::count_records(&input)?;

        Ok(Self {
            input,
            num_records,
            seed,
        })
    }

    pub(crate) fn num_records(&self) -> usize {
        self.num_records
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

    pub(crate) fn write_selected(&self, outputs: &[(&Path, usize)]) -> Result<Vec<usize>> {
        let num_selected = outputs.iter().map(|(_, count)| count).sum();
        let indices = crate::unique_random_set(num_selected, self.num_records as u32, self.seed);
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

        Ok(lengths)
    }
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
    use std::io::Write;

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

        let selector = ReadSelector::new(input.path(), Some(42)).unwrap();
        assert_eq!(selector.num_records(), 6);

        let lengths = selector
            .write_selected(&[(target.as_path(), 2), (query.as_path(), 1)])
            .unwrap();

        assert_eq!(lengths, vec![6, 3]);
        assert_eq!(
            std::fs::read_to_string(target).unwrap(),
            ">read1\nCCC\n>read2\nGGG\n"
        );
        assert_eq!(std::fs::read_to_string(query).unwrap(), ">read0\nAAA\n");

        let reads = tempdir.path().join("reads.fa");
        let lengths = selector.write_selected(&[(reads.as_path(), 3)]).unwrap();

        assert_eq!(lengths, vec![9]);
        assert_eq!(
            std::fs::read_to_string(reads).unwrap(),
            ">read0\nAAA\n>read1\nCCC\n>read2\nGGG\n"
        );
    }
}
