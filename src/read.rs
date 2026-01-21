use crate::revcomp::reverse_complement;

/// A sequence and its reverse complement
#[derive(Debug)]
pub struct Read<'a> {
    seq: &'a [u8],
    rc: Vec<u8>,
}

impl<'a> Read<'a> {
    pub fn new(seq: &'a [u8]) -> Self {
        Read {
            seq,
            rc: reverse_complement(seq),
        }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    pub fn seq(&self) -> &[u8] {
        self.seq
    }

    pub fn rc(&self) -> &[u8] {
        &self.rc
    }
}
