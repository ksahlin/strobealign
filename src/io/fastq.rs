use std::io::BufRead;

use crate::io::SequenceIOError;
use crate::io::record::{End, SequenceRecord};
use crate::io::split_header;

#[derive(Debug)]
pub struct FastqReader<B: BufRead> {
    reader: B,
    err: bool,
}

impl<B: BufRead> FastqReader<B> {
    pub fn new(reader: B) -> FastqReader<B> {
        FastqReader { reader, err: false }
    }
}

impl<B: BufRead> Iterator for FastqReader<B> {
    type Item = Result<SequenceRecord, SequenceIOError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.err {
            return None;
        }
        let mut name = String::new();
        match self.reader.read_line(&mut name) {
            Ok(0) => {
                return None;
            }
            Ok(_) => {}
            Err(e) => {
                self.err = true;
                return Some(Err(SequenceIOError::IO(e)));
            }
        }
        if !name.starts_with('@') {
            let start = name.bytes().nth(0).unwrap() as char;
            let msg = format!("Record must start with '@', but found '{}'.", start);
            self.err = true;
            return Some(Err(SequenceIOError::Fastq(msg)));
        }
        let name = name[1..].trim_end();
        let (mut name, comment) = split_header(name);

        if name.is_empty() {
            self.err = true;
            return Some(Err(SequenceIOError::Fastq(
                "Record identifier is empty".to_string(),
            )));
        }

        let mut end = End::None;
        if name.ends_with("/1") {
            name.truncate(name.len() - 2);
            end = End::One
        } else if name.ends_with("/2") {
            name.truncate(name.len() - 2);
            end = End::Two
        }

        let mut sequence = Vec::new();
        self.reader.read_until(b'\n', &mut sequence).unwrap();
        sequence.pop();
        if !sequence.is_empty() && sequence[sequence.len() - 1] == b'\r' {
            sequence.pop();
        }

        let mut name2 = String::new();
        self.reader.read_line(&mut name2).unwrap();

        let mut qualities = Vec::new();
        self.reader.read_until(b'\n', &mut qualities).unwrap();
        qualities.pop();
        if !qualities.is_empty() && qualities[qualities.len() - 1] == b'\r' {
            qualities.pop();
        }
        if sequence.len() != qualities.len() {
            let msg = format!(
                "Found {} nucleotides but {} quality values in record '{}'",
                sequence.len(),
                qualities.len(),
                name
            );
            self.err = true;
            return Some(Err(SequenceIOError::Fastq(msg)));
        }
        Some(Ok(SequenceRecord {
            name,
            end,
            comment,
            sequence: sequence.iter().map(|&c| c.to_ascii_uppercase()).collect(),
            qualities: Some(qualities),
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::{FastqReader, SequenceIOError, SequenceRecord};

    use std::io::Cursor;

    #[test]
    fn test_invalid_record_start() {
        let f = Cursor::new(b"@a\nA\n+\n#\n>b");
        let reader = FastqReader::new(f);
        let result = reader.collect::<Result<Vec<SequenceRecord>, SequenceIOError>>();

        assert!(matches!(result, Err(SequenceIOError::Fastq(_))));
    }

    #[test]
    fn test_too_few_quality_values() {
        let f = Cursor::new(b"@a\nA\n+\n#\n>b");
        let reader = FastqReader::new(f);
        let result = reader.collect::<Result<Vec<SequenceRecord>, SequenceIOError>>();

        assert!(matches!(result, Err(SequenceIOError::Fastq(_))));
    }
}
