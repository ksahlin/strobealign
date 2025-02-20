use std::collections::HashSet;
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

    #[error("Invalid UTF-8")]
    Utf8(#[from] FromUtf8Error),
    
    #[error("Invalid character in record name")]
    Name,
    
    #[error("Duplicate record name {0}")]
    DuplicateName(String),
}

/// Check whether a name is fine to use in SAM output.
/// The SAM specification is much stricter than this and forbids these
/// characters: "\'()*,<=>[\\]`{}
/// However, because even samtools itself does not complain when it encounters
/// one of them, contig names with these characters *are* used in practice, so
/// we only do some basic checks.
fn is_valid_name(name: &[u8]) -> bool{
    if name.is_empty() {
        return false;
    }
    for c in name {
        if *c < 33 || *c > 126 {
            return false;
        }
    }

    true
}

pub fn read_fasta<R: BufRead>(reader: &mut R) -> Result<Vec<RefSequence>, FastaError> {
    let mut records = Vec::<RefSequence>::new();
    let mut name = String::new();
    let mut sequence = Vec::new();
    let mut has_record = false;
    let mut names = HashSet::new();
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
            let mut name_bytes = &line[1..];
            if let Some(i) = name_bytes.iter().position(|c| c.is_ascii_whitespace()) {
                name_bytes = &name_bytes[..i];
            }
            if !is_valid_name(name_bytes) {
                return Err(FastaError::Name);
            }
            name = String::from_utf8(name_bytes.to_vec())?;
            if names.contains(&name) {
                return Err(FastaError::DuplicateName(name));
            }
            names.insert(name.clone());
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
    fn test_parse() {
        let tmp = temp_file::with_contents(b">ref1 a comment\nacgt\n\n>ref2\naacc\ngg\n\ntt\n>empty\n>empty_at_end_of_file");
        let mut reader = BufReader::new(File::open(tmp.path()).unwrap());
        let records = read_fasta(&mut reader).unwrap();
        assert_eq!(records.len(), 4);
        assert_eq!(records[0].name, "ref1");
        assert_eq!(records[0].sequence, b"ACGT");
        assert_eq!(records[1].name, "ref2");
        assert_eq!(records[1].sequence, b"AACCGGTT");
        assert_eq!(records[2].name, "empty");
        assert_eq!(records[2].sequence, b"");
        assert_eq!(records[3].name, "empty_at_end_of_file");
        assert_eq!(records[3].sequence, b"");
    }

    #[test]
    fn test_some_special_characters() {
        let mut reader = BufReader::new(b"><>;abc\nAAAA\n>abc\nCCCC\n".as_slice());
        let records = read_fasta(&mut reader).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "<>;abc");
        assert_eq!(records[1].name, "abc");
    }

    #[test]
    fn test_whitespace_in_header() {
        let mut reader = BufReader::new(b">ab c\nACGT\n>de\tf\nACGT\n".as_slice());
        let records = read_fasta(&mut reader).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "ab");
        assert_eq!(records[1].name, "de");
    }

    #[test]
    fn test_fasta_error_on_fastq() {
        let f = File::open("tests/phix.1.fastq").unwrap();
        let mut reader = BufReader::new(f);
        let result = read_fasta(&mut reader);
        assert!(result.is_err());
    }
    
    #[test]
    fn test_invalid_contig_name() {
        let mut reader = BufReader::new(b">name\x08\n".as_slice());
        let result = read_fasta(&mut reader);
        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_contig_name2() {
        let mut reader = BufReader::new(b">\0x80\nACGT\n".as_slice());
        let result = read_fasta(&mut reader);
        assert!(result.is_err());
    }
    
    #[test]
    fn test_duplicate_contig_names() {
        let mut reader = BufReader::new(b">abc\nAAAA\n>abc\nCCCC\n".as_slice());
        let result = read_fasta(&mut reader);
        assert!(result.is_err());
    }
}