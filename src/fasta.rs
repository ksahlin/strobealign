use std::io::BufRead;
use std::io;
use std::string::FromUtf8Error;
use thiserror::Error;

#[derive(Debug)]
pub struct RefSequence {
    pub name: String,
    pub sequence: Vec<u8>,
}

#[derive(Error, Debug)]
pub enum FastaError {
    #[error("IO")]
    IO(#[from] io::Error),

    #[error("FASTA file cannot be parsed: {0}")]
    Parse(String),

    #[error("Not UTF-8")]
    FromUtf8Error(#[from] FromUtf8Error),
}

pub fn read_fasta<R: BufRead>(reader: &mut R) -> Result<Vec<RefSequence>, FastaError> {
    let mut records = Vec::<RefSequence>::new();
    let mut name = String::new();
    let mut sequence = Vec::new();
    let mut has_record = false;
    for line in reader.lines() {
        let line = line?;
        let line = line.as_bytes();
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            if has_record {
                records.push(RefSequence {name, sequence});
            }
            name = String::from_utf8(line[1..].to_vec())?;
            if let Some(i) = name.find(|c: char| c.is_ascii_whitespace()) {
                name = name[..i].to_string();
            }
            sequence = Vec::new();
            has_record = true;
        } else {
            if !has_record {
                return Err(FastaError::Parse("FASTA file must start with '>'".to_string()));
            }
            sequence.extend(line.iter().map(|&c| c.to_ascii_uppercase()));
        }
    }
    if has_record {
        records.push(RefSequence {name, sequence});
    }

    Ok(records)
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::BufReader;
    use super::*;

    #[test]
    fn test_read_fasta() {
        let f = File::open("tests/phix.fasta").unwrap();
        let mut reader = BufReader::new(f);
        let records = read_fasta(&mut reader).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, "NC_001422.1");
        assert_eq!(records[0].sequence.len(), 5386);
        assert_eq!(&records[0].sequence[..5], b"GAGTT");
    }

    #[test]
    fn test_invalid_fasta() {
        let f = File::open("tests/phix.1.fastq").unwrap();
        let mut reader = BufReader::new(f);
        let result = read_fasta(&mut reader);
        assert!(result.is_err());
    }
}