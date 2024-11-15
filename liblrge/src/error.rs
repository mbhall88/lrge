use std::fmt;

/// A custom error type to represent various errors in liblrge.
#[derive(Debug)]
pub enum LrgeError {
    /// An IO error occurred.
    IoError(std::io::Error),

    /// A FASTQ parsing error occurred.
    FastqParseError(String),

    /// Too many reads were requested.
    TooManyReadsError(String),

    /// Too few reads were requested.
    TooFewReadsError(String),
}

impl fmt::Display for LrgeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LrgeError::IoError(err) => write!(f, "IO error: {}", err),
            LrgeError::FastqParseError(msg) => write!(f, "FASTQ parse error: {}", msg),
            LrgeError::TooManyReadsError(msg) => write!(f, "Too many reads requested: {}", msg),
            LrgeError::TooFewReadsError(msg) => write!(f, "Too few reads requested: {}", msg),
        }
    }
}

impl std::error::Error for LrgeError {}

/// Converts a `std::io::Error` into an [`LrgeError`].
impl From<std::io::Error> for LrgeError {
    fn from(error: std::io::Error) -> Self {
        LrgeError::IoError(error)
    }
}
