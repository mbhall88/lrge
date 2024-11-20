//! Data structure for PAF records along with serialization and deserialization methods.
use std::str::FromStr;

use serde::{Deserialize, Deserializer, Serialize, Serializer};

/// Mapping result - i.e., PafRecord
/// See https://lh3.github.io/minimap2/minimap2.html for full details of the PAF format provided by minimap2
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub(crate) struct PafRecord {
    #[serde(
        serialize_with = "serialize_bytes",
        deserialize_with = "deserialize_bytes"
    )]
    pub query_name: Vec<u8>,
    pub query_len: i32,
    /// Query start coordinate (0-based)
    pub query_start: i32,
    /// Query end coordinate (0-based)
    pub query_end: i32,
    /// ‘+’ if query/target on the same strand; ‘-’ if opposite
    pub strand: char,
    #[serde(
        serialize_with = "serialize_bytes",
        deserialize_with = "deserialize_bytes"
    )]
    pub target_name: Vec<u8>,
    pub target_len: i32,
    /// Target start coordinate on the original strand
    pub target_start: i32,
    /// Target end coordinate on the original strand
    pub target_end: i32,
    /// Number of matching bases in the mapping
    pub match_len: i32,
    /// Number bases, including gaps, in the mapping
    pub block_len: i32,
    /// Mapping quality (0-255 with 255 for missing)
    pub mapq: u32,
    /// Type of aln: P/primary, S/secondary and I,i/inversion
    #[serde(serialize_with = "serialize_tp", deserialize_with = "deserialize_tag")]
    pub tp: char,
    /// Number of minimizers on the chain
    #[serde(serialize_with = "serialize_cm", deserialize_with = "deserialize_tag")]
    pub cm: i32,
    /// Number of residues in the matching chain (chaining score)
    #[serde(serialize_with = "serialize_s1", deserialize_with = "deserialize_tag")]
    pub s1: i32,
    /// Approximate per-base sequence divergence
    #[serde(serialize_with = "serialize_dv", deserialize_with = "deserialize_tag")]
    pub dv: f32,
    /// Length of query regions harboring repetitive seeds
    #[serde(serialize_with = "serialize_rl", deserialize_with = "deserialize_tag")]
    pub rl: i32,
}

/// Serialize `Vec<u8>` as a UTF-8 string
fn serialize_bytes<S>(bytes: &[u8], serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let s = String::from_utf8_lossy(bytes);
    serializer.serialize_str(&s)
}

/// Deserialize a UTF-8 string into `Vec<u8>`
fn deserialize_bytes<'de, D>(deserializer: D) -> Result<Vec<u8>, D::Error>
where
    D: Deserializer<'de>,
{
    let s: &str = Deserialize::deserialize(deserializer)?;
    Ok(s.as_bytes().to_vec())
}

/// Serialize the tp tag
fn serialize_tp<S, T>(value: &T, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
    T: std::fmt::Display,
{
    serialize_tag_with_name("tp", value, serializer)
}

/// Serialize the cm tag
fn serialize_cm<S, T>(value: &T, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
    T: std::fmt::Display,
{
    serialize_tag_with_name("cm", value, serializer)
}

/// Serialize the s1 tag
fn serialize_s1<S, T>(value: &T, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
    T: std::fmt::Display,
{
    serialize_tag_with_name("s1", value, serializer)
}

/// Serialize the dv tag - format the float with 4 decimal places
fn serialize_dv<S>(value: &f32, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    // format the float with 4 decimal places, or if the value is zero, just serialize it as an integer
    let value = if *value < f32::EPSILON {
        "0".to_string()
    } else {
        format!("{:.4}", value)
    };
    serialize_tag_with_name("dv", &value, serializer)
}

/// Serialize the rl tag
fn serialize_rl<S, T>(value: &T, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
    T: std::fmt::Display,
{
    serialize_tag_with_name("rl", value, serializer)
}

/// Generic serialization for fields like `cm:i:123`
fn serialize_tag_with_name<S, T>(name: &str, value: &T, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
    T: std::fmt::Display,
{
    let mut prefix = match std::any::type_name::<T>() {
        "char" => "A",
        "i32" => "i",
        "f32" => "f",
        s => s,
    };

    if name == "dv" {
        prefix = "f";
    }

    let formatted = format!("{}:{}:{}", name, prefix, value);
    serializer.serialize_str(&formatted)
}

/// Generic deserialization for fields like `cm:i:123`
fn deserialize_tag<'de, T, D>(deserializer: D) -> Result<T, D::Error>
where
    T: FromStr,
    T::Err: std::fmt::Display,
    D: Deserializer<'de>,
{
    let s: &str = Deserialize::deserialize(deserializer)?;
    s.split(':')
        .last()
        .ok_or_else(|| serde::de::Error::custom("Invalid field format"))
        .and_then(|val| val.parse::<T>().map_err(serde::de::Error::custom))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_deserialize_mapping() {
        let buf = b"SRR28370649.1\t4402\t40\t237\t-\tSRR28370649.7311\t5094\t41\t238\t190\t197\t0\ttp:A:S\tcm:i:59\ts1:i:190\tdv:f:0.0022\trl:i:56";
        let expected = PafRecord {
            query_name: b"SRR28370649.1".to_vec(),
            query_len: 4402,
            query_start: 40,
            query_end: 237,
            strand: '-',
            target_name: b"SRR28370649.7311".to_vec(),
            target_len: 5094,
            target_start: 41,
            target_end: 238,
            match_len: 190,
            block_len: 197,
            mapq: 0,
            tp: 'S',
            cm: 59,
            s1: 190,
            dv: 0.0022,
            rl: 56,
        };
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(&buf[..]);
        for result in rdr.deserialize() {
            // Notice that we need to provide a type hint for automatic
            // deserialization.
            let mapping: PafRecord = result.unwrap();
            assert_eq!(mapping, expected);
        }
    }

    #[test]
    fn test_serialize_mapping() {
        let mapping = PafRecord {
            query_name: b"SRR28370649.1".to_vec(),
            query_len: 4402,
            query_start: 40,
            query_end: 237,
            strand: '-',
            target_name: b"SRR28370649.7311".to_vec(),
            target_len: 5094,
            target_start: 41,
            target_end: 238,
            match_len: 190,
            block_len: 197,
            mapq: 0,
            tp: 'S',
            cm: 59,
            s1: 190,
            dv: 0.0022,
            rl: 56,
        };
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_writer(vec![]);
        wtr.serialize(mapping).unwrap();
        let result = wtr.into_inner().unwrap();
        let result = String::from_utf8(result).unwrap();
        let expected = "SRR28370649.1\t4402\t40\t237\t-\tSRR28370649.7311\t5094\t41\t238\t190\t197\t0\ttp:A:S\tcm:i:59\ts1:i:190\tdv:f:0.0022\trl:i:56\n";
        assert_eq!(result, expected);
    }

    #[test]
    fn test_serialize_mapping_dv_round_down() {
        let mapping = PafRecord {
            query_name: b"SRR28370649.1".to_vec(),
            query_len: 4402,
            query_start: 40,
            query_end: 237,
            strand: '-',
            target_name: b"SRR28370649.7311".to_vec(),
            target_len: 5094,
            target_start: 41,
            target_end: 238,
            match_len: 190,
            block_len: 197,
            mapq: 0,
            tp: 'S',
            cm: 59,
            s1: 190,
            dv: 0.0022111,
            rl: 56,
        };
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_writer(vec![]);
        wtr.serialize(mapping).unwrap();
        let result = wtr.into_inner().unwrap();
        let result = String::from_utf8(result).unwrap();
        let expected = "SRR28370649.1\t4402\t40\t237\t-\tSRR28370649.7311\t5094\t41\t238\t190\t197\t0\ttp:A:S\tcm:i:59\ts1:i:190\tdv:f:0.0022\trl:i:56\n";
        assert_eq!(result, expected);
    }

    #[test]
    fn test_serialize_mapping_dv_round_up() {
        let mapping = PafRecord {
            query_name: b"SRR28370649.1".to_vec(),
            query_len: 4402,
            query_start: 40,
            query_end: 237,
            strand: '-',
            target_name: b"SRR28370649.7311".to_vec(),
            target_len: 5094,
            target_start: 41,
            target_end: 238,
            match_len: 190,
            block_len: 197,
            mapq: 0,
            tp: 'S',
            cm: 59,
            s1: 190,
            dv: 0.0021999,
            rl: 56,
        };
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_writer(vec![]);
        wtr.serialize(mapping).unwrap();
        let result = wtr.into_inner().unwrap();
        let result = String::from_utf8(result).unwrap();
        let expected = "SRR28370649.1\t4402\t40\t237\t-\tSRR28370649.7311\t5094\t41\t238\t190\t197\t0\ttp:A:S\tcm:i:59\ts1:i:190\tdv:f:0.0022\trl:i:56\n";
        assert_eq!(result, expected);
    }

    #[test]
    fn test_serialize_mapping_dv_fill_to_decimal_places() {
        let mapping = PafRecord {
            query_name: b"SRR28370649.1".to_vec(),
            query_len: 4402,
            query_start: 40,
            query_end: 237,
            strand: '-',
            target_name: b"SRR28370649.7311".to_vec(),
            target_len: 5094,
            target_start: 41,
            target_end: 238,
            match_len: 190,
            block_len: 197,
            mapq: 0,
            tp: 'S',
            cm: 59,
            s1: 190,
            dv: 0.004,
            rl: 56,
        };
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_writer(vec![]);
        wtr.serialize(mapping).unwrap();
        let result = wtr.into_inner().unwrap();
        let result = String::from_utf8(result).unwrap();
        let expected = "SRR28370649.1\t4402\t40\t237\t-\tSRR28370649.7311\t5094\t41\t238\t190\t197\t0\ttp:A:S\tcm:i:59\ts1:i:190\tdv:f:0.0040\trl:i:56\n";
        assert_eq!(result, expected);
    }

    #[test]
    fn test_serialize_mapping_dv_zero() {
        let mapping = PafRecord {
            query_name: b"SRR28370649.1".to_vec(),
            query_len: 4402,
            query_start: 40,
            query_end: 237,
            strand: '-',
            target_name: b"SRR28370649.7311".to_vec(),
            target_len: 5094,
            target_start: 41,
            target_end: 238,
            match_len: 190,
            block_len: 197,
            mapq: 0,
            tp: 'S',
            cm: 59,
            s1: 190,
            dv: 0.0000,
            rl: 56,
        };
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_writer(vec![]);
        wtr.serialize(mapping).unwrap();
        let result = wtr.into_inner().unwrap();
        let result = String::from_utf8(result).unwrap();
        let expected = "SRR28370649.1\t4402\t40\t237\t-\tSRR28370649.7311\t5094\t41\t238\t190\t197\t0\ttp:A:S\tcm:i:59\ts1:i:190\tdv:f:0\trl:i:56\n";
        assert_eq!(result, expected);
    }
}
