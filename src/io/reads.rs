use std::collections::VecDeque;
use std::io::{self, BufReader};

use super::SequenceIOError;
use super::fastq::FastqReader;
use super::record::{RecordPair, SequenceRecord};
use super::xopen::xopen;

type SequenceReader = Box<dyn Iterator<Item = Result<SequenceRecord, SequenceIOError>> + Send>;

pub struct PeekableSequenceReader {
    reader: SequenceReader,
    buffer: VecDeque<SequenceRecord>,
}

impl PeekableSequenceReader {
    pub fn new(reader: SequenceReader) -> Self {
        PeekableSequenceReader {
            reader,
            buffer: VecDeque::new(),
        }
    }

    // Retrieve up to n records
    pub fn peek(&mut self, n: usize) -> Result<Vec<SequenceRecord>, SequenceIOError> {
        while self.buffer.len() < n {
            if let Some(item) = self.reader.next() {
                self.buffer.push_back(item?)
            } else {
                break;
            }
        }

        Ok(self.buffer.make_contiguous().to_vec())
    }
}

impl Iterator for PeekableSequenceReader {
    type Item = Result<SequenceRecord, SequenceIOError>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(item) = self.buffer.pop_front() {
            Some(Ok(item))
        } else {
            self.reader.next()
        }
    }
}

/// Iterate over paired-end *or* single-end reads
pub fn record_iterator(
    reader1: PeekableSequenceReader,
    path_r2: Option<&str>,
) -> io::Result<Box<dyn Iterator<Item = Result<RecordPair, SequenceIOError>> + Send>> {
    if let Some(r2_path) = path_r2 {
        let fastq_reader2 = FastqReader::new(BufReader::new(xopen(r2_path)?));
        Ok(Box::new(reader1.zip(fastq_reader2).map(|p| match p {
            (Ok(r1), Ok(r2)) => Ok((r1, Some(r2))),
            (Err(e), _) => Err(e),
            (_, Err(e)) => Err(e),
        })))
    } else {
        Ok(Box::new(reader1.map(|r| r.map(|sr| (sr, None)))))
    }
}

struct InterleavedIterator {
    reader: PeekableSequenceReader,
    next_record: Option<SequenceRecord>,
}

impl InterleavedIterator {
    pub fn new(reader: PeekableSequenceReader) -> Self {
        InterleavedIterator {
            reader,
            next_record: None,
        }
    }
}

impl Iterator for InterleavedIterator {
    type Item = Result<RecordPair, SequenceIOError>;

    fn next(&mut self) -> Option<<Self as Iterator>::Item> {
        let record1 = if let Some(record) = self.next_record.take() {
            record
        } else if let Some(record) = self.reader.next() {
            match record {
                Ok(record) => record,
                Err(e) => {
                    return Some(Err(e));
                }
            }
        } else {
            return None;
        };

        match self.reader.next() {
            Some(Err(e)) => Some(Err(e)),
            Some(Ok(record2)) => {
                if record1.name != record2.name {
                    self.next_record = Some(record2);

                    Some(Ok((record1, None)))
                } else {
                    Some(Ok((record1, Some(record2))))
                }
            }
            None => Some(Ok((record1, None))),
        }
    }
}

pub fn interleaved_record_iterator(
    reader: PeekableSequenceReader,
) -> Box<dyn Iterator<Item = Result<RecordPair, SequenceIOError>> + Send> {
    Box::new(InterleavedIterator::new(reader))
}

#[cfg(test)]
mod test {
    use std::{
        fs::File,
        io::{BufReader, Cursor},
    };

    use super::*;

    #[test]
    fn test_peekable_sequence_reader() {
        let f = File::open("tests/phix.1.fastq").unwrap();
        let mut reader = PeekableSequenceReader::new(Box::new(FastqReader::new(BufReader::new(f))));

        assert_eq!(reader.next().unwrap().unwrap().name, "SRR1377138.1");
        assert_eq!(
            reader
                .peek(2)
                .unwrap()
                .iter()
                .map(|record| record.name.clone())
                .collect::<Vec<String>>(),
            vec![String::from("SRR1377138.2"), String::from("SRR1377138.3")]
        );
        assert_eq!(reader.next().unwrap().unwrap().name, "SRR1377138.2");
    }

    #[test]
    fn test_interleaved_record_iterator() {
        let f = File::open("tests/interleaved.fq").unwrap();
        let reader = PeekableSequenceReader::new(Box::new(FastqReader::new(BufReader::new(f))));
        let it = interleaved_record_iterator(reader);
        let record_pairs: Vec<RecordPair> = it
            .collect::<Result<Vec<RecordPair>, SequenceIOError>>()
            .unwrap();

        assert_eq!(record_pairs.len(), 6);
        assert_eq!(record_pairs[0].0.name, "SRR4052021.2");
        assert_eq!(record_pairs[0].1.as_ref().unwrap().name, "SRR4052021.2");
        assert!(record_pairs[0].0.sequence.starts_with(b"GTCGCCCA"));
        assert!(
            record_pairs[0]
                .1
                .as_ref()
                .unwrap()
                .sequence
                .starts_with(b"ATGTATTA")
        );

        assert_eq!(record_pairs[1].0.name, "SRR4052021.3");
        assert_eq!(record_pairs[1].1.as_ref().unwrap().name, "SRR4052021.3");

        assert_eq!(record_pairs[2].0.name, "SRR4052021.13852607");
        assert!(record_pairs[2].1.is_none());
    }

    #[test]
    fn test_interleaved_records_with_slash12() {
        // Record name a/1 followed by a/2 should be recognized as a
        // paired-end read
        let f = Cursor::new(b"@a/1\nAACC\n+\n####\n@a/2\nGGTT\n+n\n####\n");
        let reader = PeekableSequenceReader::new(Box::new(FastqReader::new(BufReader::new(f))));
        let it = interleaved_record_iterator(reader);

        let record_pairs: Vec<RecordPair> = it
            .collect::<Result<Vec<RecordPair>, SequenceIOError>>()
            .unwrap();

        assert_eq!(record_pairs.len(), 1);
        assert_eq!(record_pairs[0].0.name, "a");
        assert_eq!(record_pairs[0].1.as_ref().unwrap().name, "a");
    }
}
