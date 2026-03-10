use std::error::Error;
use std::fmt;

/// S2Error represents errors that can occur in S2 operations
#[derive(Debug)]
pub enum S2Error {
    /// Invalid loop (e.g., self-intersection)
    InvalidLoop(String),
    /// Invalid argument provided to a function
    InvalidArgument(String),
    /// Error during encoding/decoding
    EncodingError(String),
    /// Generic error with message
    Other(String),
}

impl fmt::Display for S2Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            S2Error::InvalidLoop(msg) => write!(f, "Invalid loop: {}", msg),
            S2Error::InvalidArgument(msg) => write!(f, "Invalid argument: {}", msg),
            S2Error::EncodingError(msg) => write!(f, "Encoding error: {}", msg),
            S2Error::Other(msg) => write!(f, "{}", msg),
        }
    }
}

impl From<S2Error> for S2Result<()> {
    fn from(e: S2Error) -> Self {
        Err(e)
    }
}

impl Error for S2Error {}

/// Result type for S2 operations
pub type S2Result<T> = Result<T, S2Error>;
