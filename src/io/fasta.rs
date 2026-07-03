use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::io::record::{End, SequenceRecord};
use crate::io::{SequenceIOError, split_header};
use crate::packed_seq::PackedSeq;
use crate::refseq::RefSequence;

/// Check whether a name is fine to use in SAM output.
/// The SAM specification is quite strict and forbids these
/// characters: "\'()*,<=>[\\]`{}
/// However, because even samtools itself does not complain when it encounters
/// one of them, contig names with these characters *are* used in practice, so
/// we only do some basic checks.
fn is_valid_name(name: &[u8]) -> bool {
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

fn check_duplicate_names(refseq: &RefSequence) -> Result<(), SequenceIOError> {
    let mut names = refseq.names.clone();
    names.sort();
    for window in names.windows(2) {
        if window[0] == window[1] {
            return Err(SequenceIOError::DuplicateName(window[0].to_string()));
        }
    }

    Ok(())
}

#[derive(Debug)]
pub struct FastaReader<B: BufRead> {
    reader: B,
    err: bool,
    header: Option<String>,
    sequence: Vec<u8>,
}

impl<B: BufRead> FastaReader<B> {
    pub fn new(reader: B) -> FastaReader<B> {
        FastaReader {
            reader,
            err: false,
            header: None,
            sequence: vec![],
        }
    }
}

impl<B: BufRead> Iterator for FastaReader<B> {
    type Item = Result<SequenceRecord, SequenceIOError>;

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
                            let (name, comment) = split_header(&header[1..]);
                            self.header = if n > 0 { Some(line) } else { None };
                            let record = SequenceRecord {
                                name: name.to_string(),
                                end: End::None,
                                comment,
                                sequence: std::mem::take(&mut self.sequence),
                                qualities: None,
                            };

                            return Some(Ok(record));
                        }
                        if n == 0 {
                            return None;
                        }
                        self.header = Some(line)
                    } else {
                        if self.header.is_none() {
                            return Some(Err(SequenceIOError::Fasta(
                                "FASTA file must start with '>'".to_string(),
                            )));
                        }
                        self.sequence.extend(
                            line.trim_ascii_end()
                                .bytes()
                                .map(|c| c.to_ascii_uppercase()),
                        );
                    }
                }
                Err(e) => {
                    self.err = true;
                    return Some(Err(SequenceIOError::IO(e)));
                }
            }
        }
    }
}

pub fn read_fasta<R: BufRead>(reader: &mut R) -> Result<RefSequence, SequenceIOError> {
    let mut refseq = RefSequence::new();
    let mut current_name: Option<String> = None;
    let mut current_seq = PackedSeq::new();
    let mut line = String::new();

    loop {
        line.clear();
        match reader.read_line(&mut line) {
            Ok(0) => break,
            Ok(_) => {
                if let Some(header) = line.strip_prefix('>') {
                    if let Some(name) = current_name.take() {
                        refseq.push(name, current_seq);
                        current_seq = PackedSeq::new();
                    }
                    let (name, _comment) = split_header(header);
                    current_name = Some(name.to_string());
                } else if current_name.is_some() {
                    for c in line.trim_ascii_end().bytes() {
                        current_seq.push(c);
                    }
                } else if !line.trim_ascii().is_empty() {
                    return Err(SequenceIOError::Fasta(
                        "FASTA file must start with '>'".to_string(),
                    ));
                }
            }
            Err(e) => return Err(SequenceIOError::IO(e)),
        }
    }
    if let Some(name) = current_name {
        refseq.push(name, current_seq);
    }

    for name in &refseq.names {
        if !is_valid_name(name.as_bytes()) {
            return Err(SequenceIOError::Name);
        }
    }
    check_duplicate_names(&refseq)?;

    Ok(refseq)
}

pub fn read_ref<P: AsRef<Path>>(path: P) -> Result<RefSequence, SequenceIOError> {
    let f = File::open(path).unwrap();
    let mut reader = BufReader::new(f);

    read_fasta(&mut reader)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::BufReader;

    #[test]
    fn read_ref_works() {
        let refseq = read_ref("tests/phix.fasta").unwrap();
        assert_eq!(refseq.names.len(), 1);
        assert_eq!(refseq.names[0], "NC_001422.1");
        assert_eq!(refseq.contig_len(0), 5386);
        assert_eq!(refseq.contig(0).decode(0, 5), b"GAGTT");
    }

    #[test]
    fn fasta_reader_empty_file() {
        let buf = BufReader::new(b"".as_slice());
        let reader = FastaReader::new(buf);
        let records = reader.collect::<Result<Vec<_>, _>>().unwrap();

        assert_eq!(records.len(), 0);
    }

    #[test]
    fn fasta_reader() {
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
    fn fasta_reader_only_name() {
        let mut buf = BufReader::new(b">name\n".as_slice());
        let records = FastaReader::new(&mut buf)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, "name");
        assert_eq!(records[0].sequence.len(), 0);
    }

    #[test]
    fn fasta_reader_uppercase() {
        let mut buf = BufReader::new(b">name\r\nacgt\n>name2\nTgCa".as_slice());
        let records = FastaReader::new(&mut buf)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "name");
        assert_eq!(records[0].sequence, b"ACGT");

        assert_eq!(records[1].name, "name2");
        assert_eq!(records[1].sequence, b"TGCA");
    }

    #[test]
    fn fasta_must_start_with_greater_than() {
        let mut buf = BufReader::new(b"bla".as_slice());
        let result = FastaReader::new(&mut buf).next().unwrap();
        assert!(result.is_err());
    }

    #[test]
    fn parse() {
        let tmp = temp_file::with_contents(
            b">ref1 a comment\nacgt\n\n>ref2\naacc\ngg\n\ntt\n>empty\n>empty_at_end_of_file",
        );
        let refseq = read_ref(tmp.path()).unwrap();
        assert_eq!(refseq.names.len(), 4);
        assert_eq!(refseq.names[0], "ref1");
        assert_eq!(refseq.contig(0).decode_all(), b"ACGT");
        assert_eq!(refseq.names[1], "ref2");
        assert_eq!(refseq.contig(1).decode_all(), b"AACCGGTT");
        assert_eq!(refseq.names[2], "empty");
        assert_eq!(refseq.contig(2).decode_all(), b"");
        assert_eq!(refseq.names[3], "empty_at_end_of_file");
        assert_eq!(refseq.contig(3).decode_all(), b"");
    }

    #[test]
    fn some_special_characters() {
        let mut reader = BufReader::new(b"><>;abc\nAAAA\n>abc\nCCCC\n".as_slice());
        let refseq = read_fasta(&mut reader).unwrap();
        assert_eq!(refseq.names.len(), 2);
        assert_eq!(refseq.names[0], "<>;abc");
        assert_eq!(refseq.names[1], "abc");
    }

    #[test]
    fn whitespace_in_header() {
        let mut reader = BufReader::new(b">ab c\nACGT\n>de\tf\nACGT\n".as_slice());
        let refseq = read_fasta(&mut reader).unwrap();
        assert_eq!(refseq.names.len(), 2);
        assert_eq!(refseq.names[0], "ab");
        assert_eq!(refseq.names[1], "de");
    }

    #[test]
    fn fasta_error_on_fastq() {
        assert!(read_ref("tests/phix.1.fastq").is_err());
    }

    #[test]
    fn invalid_contig_name() {
        let mut reader = BufReader::new(b">name\x08\n".as_slice());
        let result = read_fasta(&mut reader);
        assert!(result.is_err());
    }

    #[test]
    fn invalid_contig_name2() {
        let mut reader = BufReader::new(b">\0x80\nACGT\n".as_slice());
        let result = read_fasta(&mut reader);
        assert!(result.is_err());
    }

    #[test]
    fn duplicate_contig_names() {
        let mut reader = BufReader::new(b">abc\nAAAA\n>abc\nCCCC\n".as_slice());
        let result = read_fasta(&mut reader);
        assert!(result.is_err());
    }
}
