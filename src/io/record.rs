use std::fmt;

/// Which end a read is from according to the /1 or /2 suffix of its name.
///
/// Reads should normally not have this suffix, but some tools add them.
/// We don’t consider the actual /1 or /2 suffix to be part of the read name
/// and strip it. To retain the information, it is stored as part of a
/// SequenceRecord. This is only used when pairing up 'mixed' reads
/// (readname/1 followed by readname/2 is considered to be a paired-end read).
#[derive(Debug, Clone)]
pub enum End {
    /// First read in a pair
    One,

    /// Second read in a pair
    Two,

    /// Read did not have a /1 or /2 suffix in the input file
    None,
}

impl Default for End {
    fn default() -> Self {
        End::None
    }
}

#[derive(Debug, Clone)]
pub struct SequenceRecord {
    pub name: String,
    pub end: End,
    pub comment: Option<String>,
    pub sequence: Vec<u8>,
    pub qualities: Option<Vec<u8>>,
}

impl SequenceRecord {
    pub fn new(
        name: String,
        comment: Option<String>,
        sequence: Vec<u8>,
        qualities: Option<Vec<u8>>,
    ) -> Self {
        SequenceRecord {
            name,
            end: End::None,
            comment,
            sequence,
            qualities,
        }
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
        let e = vec![];
        let q = str::from_utf8(self.qualities.as_ref().unwrap_or(&e)).unwrap();
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

pub type RecordPair = (SequenceRecord, Option<SequenceRecord>);
