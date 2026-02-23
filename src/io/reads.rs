use std::collections::VecDeque;
use std::io::{self, BufRead, BufReader};

use super::SequenceIOError;
use super::fastq::FastqReader;
use super::record::{RecordPair, SequenceRecord};
use super::xopen::xopen;

#[derive(Debug)]
pub struct PeekableSequenceReader<B: BufRead + Send> {
    fastq_reader: FastqReader<B>,
    buffer: VecDeque<SequenceRecord>,
}

impl<B: BufRead + Send> PeekableSequenceReader<B> {
    pub fn new(reader: B) -> Self {
        let fastq_reader = FastqReader::new(reader);
        PeekableSequenceReader {
            fastq_reader,
            buffer: VecDeque::new(),
        }
    }

    // Retrieve up to n records
    pub fn peek(&mut self, n: usize) -> Result<Vec<SequenceRecord>, SequenceIOError> {
        while self.buffer.len() < n {
            if let Some(item) = self.fastq_reader.next() {
                self.buffer.push_back(item?)
            } else {
                break;
            }
        }

        Ok(self.buffer.make_contiguous().to_vec())
    }
}

impl<B: BufRead + Send> Iterator for PeekableSequenceReader<B> {
    type Item = Result<SequenceRecord, SequenceIOError>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(item) = self.buffer.pop_front() {
            Some(Ok(item))
        } else {
            self.fastq_reader.next()
        }
    }
}

/// Iterate over paired-end *or* single-end reads
pub fn record_iterator<'a, B: BufRead + Send + 'a>(
    fastq_reader1: PeekableSequenceReader<B>,
    path_r2: Option<&str>,
) -> io::Result<Box<dyn Iterator<Item = Result<RecordPair, SequenceIOError>> + Send + 'a>> {
    if let Some(r2_path) = path_r2 {
        let fastq_reader2 = FastqReader::new(BufReader::new(xopen(r2_path)?));
        Ok(Box::new(fastq_reader1.zip(fastq_reader2).map(
            |p| match p {
                (Ok(r1), Ok(r2)) => Ok((r1, Some(r2))),
                (Err(e), _) => Err(e),
                (_, Err(e)) => Err(e),
            },
        )))
    } else {
        Ok(Box::new(fastq_reader1.map(|r| r.map(|sr| (sr, None)))))
    }
}

struct InterleavedIterator<B: BufRead + Send> {
    fastq_reader: PeekableSequenceReader<B>,
    next_record: Option<SequenceRecord>,
}

impl<B: BufRead + Send> InterleavedIterator<B> {
    pub fn new(fastq_reader: PeekableSequenceReader<B>) -> Self {
        InterleavedIterator {
            fastq_reader,
            next_record: None,
        }
    }
}

impl<B: BufRead + Send> Iterator for InterleavedIterator<B> {
    type Item = Result<RecordPair, SequenceIOError>;

    fn next(&mut self) -> Option<<Self as Iterator>::Item> {
        let record1 = if let Some(record) = self.next_record.take() {
            record
        } else if let Some(record) = self.fastq_reader.next() {
            match record {
                Ok(record) => record,
                Err(e) => {
                    return Some(Err(e));
                }
            }
        } else {
            return None;
        };

        match self.fastq_reader.next() {
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

pub fn interleaved_record_iterator<'a, B: BufRead + Send + 'a>(
    fastq_reader: PeekableSequenceReader<B>,
) -> Box<dyn Iterator<Item = Result<RecordPair, SequenceIOError>> + Send + 'a> {
    Box::new(InterleavedIterator::new(fastq_reader))
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
        let mut reader = PeekableSequenceReader::new(BufReader::new(f));

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
        let reader = PeekableSequenceReader::new(BufReader::new(f));
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
        let reader = PeekableSequenceReader::new(f);
        let it = interleaved_record_iterator(reader);

        let record_pairs: Vec<RecordPair> = it
            .collect::<Result<Vec<RecordPair>, SequenceIOError>>()
            .unwrap();

        assert_eq!(record_pairs.len(), 1);
        assert_eq!(record_pairs[0].0.name, "a");
        assert_eq!(record_pairs[0].1.as_ref().unwrap().name, "a");
    }
}
