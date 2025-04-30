use std::io::{BufRead, BufReader, Read};
use std::io;
use std::string::FromUtf8Error;
use thiserror::Error;
use crate::fastq::SequenceRecord;

#[derive(Debug, Clone)]
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
/// The SAM specification is quite strict and forbids these
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

fn check_duplicate_names(records: &[RefSequence]) -> Result<(), FastaError> {
    let mut names: Vec<_> = records.iter().map(|r| &r.name).collect();
    names.sort();
    for window in names.windows(2) {
        if window[0] == window[1] {
            return Err(FastaError::DuplicateName(window[0].to_string()));
        }
    }

    Ok(())
}

#[derive(Debug)]
struct FastaReader<R> {
    reader: BufReader<R>,
    err: bool,
    header: Option<String>,
    sequence: Vec<u8>,
}

impl<R: Read> FastaReader<R> {
    pub fn new(reader: R) -> FastaReader<R> {
        FastaReader {
            reader: BufReader::new(reader),
            err: false,
            header: None,
            sequence: vec![],
        }
    }
}

impl<R: Read> Iterator for FastaReader<R> {
    type Item = Result<SequenceRecord, FastaError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.err {
            return None;
        }
        loop {
            let mut line = String::new();
            match self.reader.read_line(&mut line) {
                Ok(n) => {
                    if n == 0 || line.starts_with('>') {
                        if let Some(header) = &self.header {
                            let (name, comment) = split_header(header[1..].trim_ascii_end());
                            self.header = if n > 0 { Some(line) } else { None };
                            let record = SequenceRecord {
                                name: name.to_string(),
                                comment,
                                sequence: std::mem::take(&mut self.sequence),
                                qualities: None,
                            };

                            return Some(Ok(record))
                        }
                        if n == 0 {
                            return None;
                        }
                        self.header = Some(line)
                    } else {
                        if self.header.is_none() {
                            return Some(Err(FastaError::Parse("FASTA file must start with '>'".to_string())));
                        }
                        self.sequence.extend(line.trim_ascii_end().bytes().map(|c| c.to_ascii_uppercase()));
                    }
                }
                Err(e) => {
                    self.err = true;
                    return Some(Err(FastaError::IO(e)));
                }
            }
        }
    }
}

/// Split header into name and comment
pub fn split_header(header: &str) -> (String, Option<String>) {
    match header.split_once(' ') {
        Some((name, comment)) => (name.to_string(), Some(comment.to_string())),
        None => (header.to_string(), None),
    }
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
            let mut name_bytes = &line[1..];
            if let Some(i) = name_bytes.iter().position(|c| c.is_ascii_whitespace()) {
                name_bytes = &name_bytes[..i];
            }
            if !is_valid_name(name_bytes) {
                return Err(FastaError::Name);
            }
            name = String::from_utf8(name_bytes.to_vec())?;
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

    check_duplicate_names(&records)?;

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
    fn test_fasta_reader_empty_file() {
        let buf = BufReader::new(b"".as_slice());
        let reader = FastaReader::new(buf);
        let records = reader.collect::<Result<Vec<_>, _>>().unwrap();

        assert_eq!(records.len(), 0);
    }

    #[test]
    fn test_fasta_reader() {
        let f = File::open("tests/phix.fasta").unwrap();
        let mut bufreader = BufReader::new(f);
        let reader = FastaReader::new(&mut bufreader);

        let records = reader.collect::<Result<Vec<_>, _>>().unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, "NC_001422.1");
        assert_eq!(records[0].sequence.len(), 5386);
        assert_eq!(&records[0].sequence[..5], b"GAGTT");
    }

    #[test]
    fn test_fasta_reader_only_name() {
        let mut buf = BufReader::new(b">name\n".as_slice());
        let records = FastaReader::new(&mut buf).collect::<Result<Vec<_>, _>>().unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, "name");
        assert_eq!(records[0].sequence.len(), 0);
    }

    #[test]
    fn test_fasta_reader_uppercase() {
        let mut buf = BufReader::new(b">name\r\nacgt\n>name2\nTgCa".as_slice());
        let records = FastaReader::new(&mut buf).collect::<Result<Vec<_>, _>>().unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "name");
        assert_eq!(records[0].sequence, b"ACGT");
        
        assert_eq!(records[1].name, "name2");
        assert_eq!(records[1].sequence, b"TGCA");
    }

    #[test]
    fn test_fasta_must_start_with_greater_than() {
        let mut buf = BufReader::new(b"bla".as_slice());
        let result = FastaReader::new(&mut buf).next().unwrap();
        assert!(result.is_err());
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
