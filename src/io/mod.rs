pub mod fasta;
pub mod fastq;
pub mod paf;
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
