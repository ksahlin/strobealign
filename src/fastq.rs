use std::fmt;
use std::io::{BufRead, BufReader, BufWriter, Read, Result, Write};
use std::str;

#[derive(Debug, Clone)]
pub struct SequenceRecord {
    pub name: String,
    pub comment: Option<String>,
    pub sequence: Vec<u8>,
    pub qualities: Vec<u8>,
}

impl SequenceRecord {
    pub fn new(name: String, comment: Option<String>, sequence: Vec<u8>, qualities: Vec<u8>) -> Self {
        SequenceRecord { name, comment, sequence, qualities }
    }

    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

impl fmt::Display for SequenceRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s = str::from_utf8(&self.sequence).unwrap();
        let q = str::from_utf8(&self.qualities).unwrap();
        write!(
            f,
            "name={}, length={}, sequence={}, qualities={}",
            self.name,
            self.len(),
            s,
            q
        )
    }
}

#[derive(Debug)]
pub struct FastqReader<R: Read> {
    reader: BufReader<R>,
    err: bool,
}

impl<R: Read> FastqReader<R> {
    pub fn new(reader: R) -> FastqReader<R> {
        FastqReader {
            reader: BufReader::new(reader),
            err: false,
        }
    }
}

impl<R: Read> Iterator for FastqReader<R> {
    type Item = Result<SequenceRecord>;

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
                return Some(Err(e));
            }
        }
        let name = name[1..].trim_end();
        let (mut name, comment) = match name.split_once(' ') {
            Some((name, comment)) => (name, Some(comment)),
            None => (name, None),
        };
        if name.ends_with("/1") || name.ends_with("/2") {
            name = &name[..name.len() - 2];
        }
        let name = name;

        if name.is_empty() {
            //self.err = true;
            //return Some(Err("error"));
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
        assert_eq!(sequence.len(), qualities.len());
        Some(Ok(SequenceRecord {
            name: name.to_string(),
            comment: comment.map(|s| s.to_string()),
            sequence: sequence.iter().map(|&c| c.to_ascii_uppercase()).collect(),
            qualities,
        }))
    }
}

pub struct FastqWriter<W: Write> {
    writer: BufWriter<W>,
}

impl<W: Write> FastqWriter<W> {
    pub fn new(writer: W) -> Self {
        FastqWriter {
            writer: BufWriter::new(writer),
        }
    }

    pub fn write_record(&mut self, record: &SequenceRecord) {
        self.writer.write_all(b"@").unwrap();
        self.writer.write_all(record.name.as_bytes()).unwrap();
        self.writer.write_all(b"\n").unwrap();
        self.writer.write_all(&record.sequence).unwrap();
        self.writer.write_all(b"\n+\n").unwrap();
        self.writer.write_all(&record.qualities).unwrap();
        self.writer.write_all(b"\n").unwrap();
        //write!(self.writer, "@{}\n{:?}\n+\n{:?}\n", record.name, record.sequence, record.qualities);
    }
}
