use std::collections::VecDeque;
use std::fmt;
use std::io::{BufRead, BufReader, BufWriter, Read, Result, Write};
use std::str;
use crate::io::xopen;

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
pub struct PeekableFastqReader<R: Read + Send> {
    fastq_reader: FastqReader<R>,
    buffer: VecDeque<SequenceRecord>,
}

impl<R: Read + Send> PeekableFastqReader<R> {
    pub fn new(reader: R) -> PeekableFastqReader<R> {
        let fastq_reader = FastqReader::new(reader);
        PeekableFastqReader {
            fastq_reader,
            buffer: VecDeque::new(),
        }
    }

    // Retrieve up to n records
    pub fn peek(&mut self, n: usize) -> Result<Vec<SequenceRecord>> {
        while self.buffer.len() < n {
            if let Some(item) = self.fastq_reader.next() {
                self.buffer.push_back(item?)
            } else {
                break;
            }
        }

        Ok(self.buffer.make_contiguous().iter().cloned().collect())
    }
}

impl<R: Read + Send> Iterator for PeekableFastqReader<R> {
    type Item = Result<SequenceRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(item) = self.buffer.pop_front() {
            Some(Ok(item))
        } else {
            self.fastq_reader.next()
        }
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
        let (name, comment) = match name.split_once(' ') {
            Some((name, comment)) => (name, Some(comment)),
            None => (name, None),
        };

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

/// Iterate over paired-end *or* single-end reads
pub fn record_iterator<'a, R: Read + Send + 'a>(
    fastq_reader1: PeekableFastqReader<R>, path_r2: Option<&str>
) -> Result<Box<dyn Iterator<Item=Result<(SequenceRecord, Option<SequenceRecord>)>> + Send + 'a>> 
{
    if let Some(r2_path) = path_r2 {
        let fastq_reader2 = FastqReader::new(xopen(r2_path)?);
        Ok(
            Box::new(
                fastq_reader1.zip(fastq_reader2).map(
                    |p|
                        match p {
                            (Ok(r1), Ok(r2)) => { Ok((r1, Some(r2))) }
                            (Err(e), _) => Err(e),
                            (_, Err(e)) => Err(e),
                        }
                )
            )
        )
    } else {
        Ok(
            Box::new(
                fastq_reader1.map(
                    |r| r.map(|sr| (sr, None))
                )
            )
        )
    }
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use crate::fastq::PeekableFastqReader;

    #[test]
    fn test_peekable_fastq_reader() {
        let f = File::open("tests/phix.1.fastq").unwrap();
        let mut reader = PeekableFastqReader::new(f);

        assert_eq!(reader.next().unwrap().unwrap().name, "SRR1377138.1");
        assert_eq!(
            reader.peek(2).unwrap().iter().map(|record| record.name.clone()).collect::<Vec<String>>(),
            vec![String::from("SRR1377138.2"), String::from("SRR1377138.3/1")]
        );
        assert_eq!(reader.next().unwrap().unwrap().name, "SRR1377138.2");
    }
}