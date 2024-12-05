use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::Path;

#[cfg(feature = "bzip2")]
use bzip2::bufread::BzDecoder;
#[cfg(feature = "gzip")]
use flate2::bufread::MultiGzDecoder;
#[cfg(feature = "xz")]
use liblzma::read::XzDecoder;
use needletail::parse_fastx_reader;
#[cfg(feature = "zstd")]
use zstd::stream::read::Decoder as ZstdDecoder;

/// The compression format of a file.
#[derive(Debug, PartialEq, Copy, Clone, Default)]
enum CompressionFormat {
    #[cfg(feature = "bzip2")]
    Bzip2,
    #[cfg(feature = "gzip")]
    Gzip,
    #[default]
    None,
    #[cfg(feature = "xz")]
    Xz,
    #[cfg(feature = "zstd")]
    Zstd,
}

/// Detects the compression format of a file by reading the magic bytes at the start of the file.
fn detect_compression_format<R: Read + Seek>(reader: &mut R) -> io::Result<CompressionFormat> {
    let original_position = reader.stream_position()?;

    // move the reader to the start of the file
    reader.seek(SeekFrom::Start(0))?;

    let mut magic = [0; 5];
    reader
        .read_exact(&mut magic)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let format = match magic {
        #[cfg(feature = "gzip")]
        [0x1f, 0x8b, ..] => CompressionFormat::Gzip,
        #[cfg(feature = "bzip2")]
        [0x42, 0x5a, ..] => CompressionFormat::Bzip2,
        #[cfg(feature = "zstd")]
        [0x28, 0xb5, 0x2f, 0xfd, ..] => CompressionFormat::Zstd,
        #[cfg(feature = "xz")]
        [0xfd, 0x37, 0x7a, 0x58, 0x5a] => CompressionFormat::Xz,
        _ => CompressionFormat::None,
    };

    // Seek back to the original position
    reader.seek(SeekFrom::Start(original_position))?;

    Ok(format)
}

/// Opens a file and returns a reader. Supports gzip and zstd compression if the corresponding
/// feature is enabled. If the file is not compressed, a regular file reader is returned. If the
/// file is compressed with an unsupported format, an error is returned.
pub(crate) fn open_file<P: AsRef<Path>>(path: P) -> io::Result<Box<dyn Read + Send>> {
    let mut buf = File::open(&path).map(BufReader::new)?;
    let compression_format = detect_compression_format(&mut buf)?;

    let reader: Box<dyn Read + Send> = match compression_format {
        #[cfg(feature = "gzip")]
        CompressionFormat::Gzip => Box::new(MultiGzDecoder::new(buf)),

        #[cfg(feature = "zstd")]
        CompressionFormat::Zstd => Box::new(ZstdDecoder::new(buf)?),

        #[cfg(feature = "bzip2")]
        CompressionFormat::Bzip2 => Box::new(BzDecoder::new(buf)),

        #[cfg(feature = "xz")]
        CompressionFormat::Xz => Box::new(XzDecoder::new(buf)),

        CompressionFormat::None => Box::new(buf),
    };

    Ok(reader)
}

pub(crate) fn count_fastq_records<R: Read + Send>(reader: R) -> io::Result<usize> {
    let mut count = 0;

    let mut fastx_reader =
        parse_fastx_reader(reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    while let Some(r) = fastx_reader.next() {
        let _ = r.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        count += 1;
    }

    Ok(count)
}

/// A message that can be sent in a channel.
pub(crate) enum Message {
    /// The intention is to send a read ID and a read sequence.
    Data((Vec<u8>, Vec<u8>)),
}

pub(crate) trait FastqRecordExt {
    fn read_id(&self) -> &[u8];
}

impl FastqRecordExt for needletail::parser::SequenceRecord<'_> {
    /// The needletail FastxRecord `id` method returns the whole header line, including the comment
    /// and the read ID. This method returns only the read ID.
    fn read_id(&self) -> &[u8] {
        let id = self.id();
        id.split(|&x| x.is_ascii_whitespace())
            .next()
            .unwrap_or_default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_detect_gzip_format() {
        let data = vec![
            0x1f, 0x8b, 0x08, 0x08, 0x1c, 0x6b, 0xe2, 0x66, 0x00, 0x03, 0x74, 0x65, 0x78, 0x74,
            0x2e, 0x74, 0x78, 0x74, 0x00, 0x4b, 0xcb, 0xcf, 0x57, 0x48, 0x4a, 0x2c, 0xe2, 0x02,
            0x00, 0x27, 0xb4, 0xdd, 0x13, 0x08, 0x00, 0x00, 0x00,
        ];
        let mut reader = Cursor::new(data);
        // position the reader at the original position
        let original_position = reader.position();
        let format = detect_compression_format(&mut reader).unwrap();
        assert_eq!(format, CompressionFormat::Gzip);
        assert_eq!(reader.position(), original_position);
    }

    #[test]
    fn test_detect_bzip2_format() {
        let data = vec![
            0x42, 0x5a, 0x68, 0x39, 0x31, 0x41, 0x59, 0x26, 0x53, 0x59, 0x7b, 0x6e, 0xa8, 0x38,
            0x00, 0x00, 0x02, 0x51, 0x80, 0x00, 0x10, 0x40, 0x00, 0x31, 0x00, 0x90, 0x00, 0x20,
            0x00, 0x22, 0x1a, 0x63, 0x50, 0x86, 0x00, 0x2c, 0x8c, 0x3c, 0x5d, 0xc9, 0x14, 0xe1,
            0x42, 0x41, 0xed, 0xba, 0xa0, 0xe0,
        ];
        let mut reader = Cursor::new(data);
        // position the reader at the original position
        let original_position = reader.position();
        let format = detect_compression_format(&mut reader).unwrap();
        assert_eq!(format, CompressionFormat::Bzip2);
        assert_eq!(reader.position(), original_position);
    }

    #[test]
    fn test_detect_zstd_format() {
        let data = vec![
            0x28, 0xb5, 0x2f, 0xfd, 0x24, 0x08, 0x41, 0x00, 0x00, 0x66, 0x6f, 0x6f, 0x20, 0x62,
            0x61, 0x72, 0x0a, 0x37, 0x17, 0xa5, 0xec,
        ];
        let mut reader = Cursor::new(data);
        // position the reader at the original position
        let original_position = reader.position();
        let format = detect_compression_format(&mut reader).unwrap();
        assert_eq!(format, CompressionFormat::Zstd);
        assert_eq!(reader.position(), original_position);
    }

    #[test]
    fn test_detect_xz_format() {
        let data = vec![
            0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00, 0x00, 0x04, 0xe6, 0xd6, 0xb4, 0x46, 0x02, 0x00,
            0x21, 0x01, 0x16, 0x00, 0x00, 0x00, 0x74, 0x2f, 0xe5, 0xa3, 0x01, 0x00, 0x07, 0x66,
            0x6f, 0x6f, 0x20, 0x62, 0x61, 0x72, 0x0a, 0x00, 0xfd, 0xbb, 0xfb, 0x3b, 0x8e, 0xcc,
            0x32, 0x13, 0x00, 0x01, 0x20, 0x08, 0xbb, 0x19, 0xd9, 0xbb, 0x1f, 0xb6, 0xf3, 0x7d,
            0x01, 0x00, 0x00, 0x00, 0x00, 0x04, 0x59, 0x5a,
        ];
        let mut reader = Cursor::new(data);
        // position the reader at the original position
        let original_position = reader.position();
        let format = detect_compression_format(&mut reader).unwrap();
        assert_eq!(format, CompressionFormat::Xz);

        // confirm that the reader is still at the original position
        assert_eq!(reader.position(), original_position);
    }

    #[test]
    fn test_detect_none_format() {
        let data = b"I'm not compressed";
        let mut reader = Cursor::new(data);
        let format = detect_compression_format(&mut reader).unwrap();
        assert_eq!(format, CompressionFormat::None);
    }

    #[test]
    fn test_detect_format_when_reader_is_part_way_through() {
        let data = vec![
            0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00, 0x00, 0x04, 0xe6, 0xd6, 0xb4, 0x46, 0x02, 0x00,
            0x21, 0x01, 0x16, 0x00, 0x00, 0x00, 0x74, 0x2f, 0xe5, 0xa3, 0x01, 0x00, 0x07, 0x66,
            0x6f, 0x6f, 0x20, 0x62, 0x61, 0x72, 0x0a, 0x00, 0xfd, 0xbb, 0xfb, 0x3b, 0x8e, 0xcc,
            0x32, 0x13, 0x00, 0x01, 0x20, 0x08, 0xbb, 0x19, 0xd9, 0xbb, 0x1f, 0xb6, 0xf3, 0x7d,
            0x01, 0x00, 0x00, 0x00, 0x00, 0x04, 0x59, 0x5a,
        ];
        let mut reader = Cursor::new(data);
        reader.seek(SeekFrom::Start(10)).unwrap();
        // position the reader at the original position
        let original_position = reader.position();
        let format = detect_compression_format(&mut reader).unwrap();
        assert_eq!(format, CompressionFormat::Xz);

        // confirm that the reader is still at the original position
        assert_eq!(reader.position(), original_position);
    }

    #[test]
    fn test_count_fastq_records() {
        let data = b"@SEQ_ID\nGATTA\n+\n!!!!!\n@SEQ_ID2\nGATTA\n+\n!!!!!\n";
        let reader = Cursor::new(data);
        let count = count_fastq_records(reader).unwrap();
        assert_eq!(count, 2);
    }

    #[test]
    fn test_count_fastq_records_with_empty_file() {
        let data = b"";
        let reader = Cursor::new(data);
        let err = count_fastq_records(reader).unwrap_err();
        assert!(err.to_string().contains("Is the file empty?"));
    }

    #[test]
    fn test_count_fastq_records_with_invalid_data() {
        let data = b"@SEQ_ID\nGATTA\n+\n!!!!\n@SEQ_ID2\nGATTA\n+\n!!!!!";
        let reader = Cursor::new(data);
        let err = count_fastq_records(reader).unwrap_err();
        assert!(err.to_string().contains("but quality length is"));
    }

    #[test]
    fn test_read_id_no_comment() {
        let data = b"@SEQ_ID\nGATTA\n+\n!!!!!\n";
        let reader = Cursor::new(data);
        let mut fastx_reader = parse_fastx_reader(reader).unwrap();
        let record = fastx_reader.next().unwrap().unwrap();
        assert_eq!(record.read_id(), b"SEQ_ID");
    }

    #[test]
    fn test_read_id_with_comment() {
        let data = b"@SEQ_ID comment\nGATTA\n+\n!!!!!\n";
        let reader = Cursor::new(data);
        let mut fastx_reader = parse_fastx_reader(reader).unwrap();
        let record = fastx_reader.next().unwrap().unwrap();
        assert_eq!(record.read_id(), b"SEQ_ID");
    }

    #[test]
    fn test_read_id_with_empty_comment() {
        let data = b"@SEQ_ID \nGATTA\n+\n!!!!!\n";
        let reader = Cursor::new(data);
        let mut fastx_reader = parse_fastx_reader(reader).unwrap();
        let record = fastx_reader.next().unwrap().unwrap();
        assert_eq!(record.read_id(), b"SEQ_ID");
    }

    #[test]
    fn test_read_id_with_multiple_spaces() {
        let data = b"@SEQ_ID   comment\nGATTA\n+\n!!!!!\n";
        let reader = Cursor::new(data);
        let mut fastx_reader = parse_fastx_reader(reader).unwrap();
        let record = fastx_reader.next().unwrap().unwrap();
        assert_eq!(record.read_id(), b"SEQ_ID");
    }

    #[test]
    fn test_read_id_with_tabs() {
        let data = b"@SEQ_ID\tst:Z:2024-06-05T11:34:21.517+00:00\tRG:Z:0e9626940687df5718807f8d3dcf3c2d2b2e49c6_dna_r10.4.1_e8.2_400bps_sup@v5.0.0_SQK-RBK114-96_barcode58\nGATTA\n+\n!!!!!\n";
        let reader = Cursor::new(data);
        let mut fastx_reader = parse_fastx_reader(reader).unwrap();
        let record = fastx_reader.next().unwrap().unwrap();
        assert_eq!(record.read_id(), b"SEQ_ID");
    }
}
