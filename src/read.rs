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
            rc: reverse_complement(&seq),
        }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn seq(&self) -> &[u8] {
        &self.seq
    }

    pub fn rc(&self) -> &[u8] {
        &self.rc
    }
}
