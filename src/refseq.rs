use crate::packed_seq::{PackedSeq, PackedSeqSlice};

/// A position on a contig. This separate type is here to prevent confusion
/// with "flat" reference coordinates, which are used everywhere else.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ContigPosition(pub u32);

#[derive(Debug, Clone)]
pub struct ContigStarts(
    /// Start positions for all contigs. `starts[i]` is the start position of
    /// contig with index `i` within `sequence`. This contains one more item
    /// than there are contigs, which is set to the length of `sequence`.
    pub Vec<usize>,
);

impl Default for ContigStarts {
    fn default() -> Self {
        Self(vec![0])
    }
}

/// Contig start positions within the concatenated reference.
///
/// Use this only to lookup *start* coordinates.
/// Otherwise, when trying to lookup an end coordinate that happens to be
/// the end of a contig i, this is found as position 0 of contig i + 1.
impl ContigStarts {
    pub fn new(mut starts: Vec<usize>, total_length: usize) -> Self {
        starts.push(total_length);

        Self(starts)
    }

    /// Returns the index of the contig that contains the given position.
    pub fn index(&self, start: usize) -> usize {
        self.0[1..].partition_point(|&x| x <= start)
    }

    /// Returns the start position of the contig that contains the given position.
    pub fn ref_contig_start(&self, start: usize) -> usize {
        self.0[self.index(start)]
    }

    /// Returns (contig_index, contig_start)
    pub fn unflatten(&self, start: usize) -> (usize, ContigPosition) {
        let ref_index = self.index(start);
        let start_within_contig = start - self.0[ref_index];

        (ref_index, ContigPosition(start_within_contig as u32))
    }

    pub fn total_length(&self) -> usize {
        *self.0.last().unwrap()
    }
}

#[derive(Default, Debug, Clone)]
pub struct RefSequence {
    /// Contig names
    pub names: Vec<String>,

    /// Concatenated sequence of all contigs
    sequence: PackedSeq,

    pub starts: ContigStarts,
}

impl RefSequence {
    pub fn new(sequence: PackedSeq, starts: Vec<usize>, names: Vec<String>) -> Self {
        assert_eq!(starts.len(), names.len());
        let total_length = sequence.len();

        RefSequence {
            sequence,
            starts: ContigStarts::new(starts, total_length),
            names,
        }
    }

    pub fn sequence(&self) -> &PackedSeq {
        &self.sequence
    }

    /// Returns (contig_index, contig_start)
    pub fn unflatten(&self, start: usize) -> (usize, ContigPosition) {
        self.starts.unflatten(start)
    }

    pub fn decode(&self, start: usize, end: usize) -> Vec<u8> {
        self.sequence.decode(start, end)
    }

    pub fn contig<'a>(&'a self, index: usize) -> PackedSeqSlice<'a> {
        let start = self.starts.0[index];
        let len = self.starts.0[index + 1] - start;

        PackedSeqSlice::new(&self.sequence, start, len)
    }

    pub fn contig_start(&self, index: usize) -> usize {
        self.starts.0[index]
    }

    pub fn max_contig_len(&self) -> Option<usize> {
        (0..self.names.len()).map(|i| self.contig(i).len()).max()
    }

    pub fn total_length(&self) -> usize {
        self.starts.total_length()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn testref() -> RefSequence {
        let mut sequence = PackedSeq::new();
        sequence.extend(b"AAAACCCGGT".to_vec());
        let names = vec![
            "n1".to_string(),
            "n2".to_string(),
            "n3".to_string(),
            "n4".to_string(),
        ];
        let starts = vec![0, 4, 7, 9];
        RefSequence::new(sequence, starts, names)
    }

    #[test]
    fn new() {
        let rs = testref();

        assert_eq!(rs.total_length(), 10);
        assert_eq!(rs.max_contig_len(), Some(4));

        assert_eq!(rs.contig(0).len(), 4);
        assert_eq!(rs.contig(1).len(), 3);
        assert_eq!(rs.contig(2).len(), 2);
        assert_eq!(rs.contig(3).len(), 1);
        assert_eq!(rs.contig(0).len(), 4);
    }

    #[test]
    fn decode() {
        let rs = testref();

        assert_eq!(rs.decode(0, 0), b"");
        assert_eq!(rs.decode(0, 1), b"A");
        assert_eq!(rs.decode(0, 2), b"AA");
        assert_eq!(rs.decode(0, 3), b"AAA");
        assert_eq!(rs.decode(4, 5), b"C");
        assert_eq!(rs.decode(4, 6), b"CC");
        assert_eq!(rs.decode(4, 7), b"CCC");
        assert_eq!(rs.decode(9, 10), b"T");
    }

    #[test]
    fn unflatten() {
        let rs = testref();

        assert_eq!(rs.unflatten(0), (0, ContigPosition(0)));
        assert_eq!(rs.unflatten(1), (0, ContigPosition(1)));
        assert_eq!(rs.unflatten(2), (0, ContigPosition(2)));
        assert_eq!(rs.unflatten(3), (0, ContigPosition(3)));
        assert_eq!(rs.unflatten(4), (1, ContigPosition(0)));
        assert_eq!(rs.unflatten(5), (1, ContigPosition(1)));
        assert_eq!(rs.unflatten(6), (1, ContigPosition(2)));
        assert_eq!(rs.unflatten(7), (2, ContigPosition(0)));
        assert_eq!(rs.unflatten(8), (2, ContigPosition(1)));
        assert_eq!(rs.unflatten(9), (3, ContigPosition(0)));
    }
}
