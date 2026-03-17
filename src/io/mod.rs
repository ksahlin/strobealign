pub mod fasta;
pub mod fastq;
pub mod paf;
pub mod reads;
pub mod record;
pub mod sam;
pub mod xopen;

use std::io;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum SequenceIOError {
    #[error("IO")]
    IO(#[from] io::Error),

    #[error("FASTA file cannot be parsed: {0}")]
    Fasta(String),

    #[error("FASTQ file cannot be parsed: {0}")]
    Fastq(String),

    #[error("Invalid character in record name")]
    Name,

    #[error("Duplicate record name {0}")]
    DuplicateName(String),
}

/// Split header into name and comment
pub fn split_header(header: &str) -> (String, Option<String>) {
    let header = header.trim_ascii_end();
    match header.split_once(&[' ', '\t']) {
        Some((name, comment)) => (name.to_string(), Some(comment.to_string())),
        None => (header.to_string(), None),
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_split_header() {
        assert_eq!(
            split_header("id comment"),
            ("id".to_string(), Some("comment".to_string()))
        );
        assert_eq!(
            split_header("id comment "),
            ("id".to_string(), Some("comment".to_string()))
        );
        assert_eq!(split_header("id"), ("id".to_string(), None));
        assert_eq!(split_header("id "), ("id".to_string(), None));
        assert_eq!(split_header("id  "), ("id".to_string(), None));
    }
}
