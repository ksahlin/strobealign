use crate::packed_seq::{PackedSeq, PackedSeqSlice};

/// A position on a contig. This separate type is here to prevent confusion
/// with "flat" reference coordinates, which are used everywhere else.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ContigPosition(pub u32);

#[derive(Default, Debug, Clone)]
pub struct RefSequence {
    /// Contig names
    pub names: Vec<String>,

    /// Concatenated sequence of all contigs
    sequence: PackedSeq,

    /// Start positions for all contigs. `starts[i]` is the start position of
    /// contig with index `i` within `sequence`. This contains one more item
    /// than there are contigs, which is set to the length of `sequence`.
    starts: Vec<usize>,
}

impl RefSequence {
    pub fn new(sequence: PackedSeq, starts: Vec<usize>, names: Vec<String>) -> Self {
        assert_eq!(starts.len(), names.len());
        let total_length = sequence.len();
        let mut refseq = RefSequence {
            sequence,
            starts,
            names,
        };
        refseq.starts.push(total_length);

        refseq
    }

    pub fn sequence(&self) -> &PackedSeq {
        &self.sequence
    }

    /// Returns (reference_index, reference_start)
    pub fn unflatten(&self, start: usize) -> (usize, ContigPosition) {
        let ref_index = self.starts[1..].partition_point(|&x| x <= start);
        let start_within_contig = start - self.starts[ref_index];

        (ref_index, ContigPosition(start_within_contig as u32))
    }

    pub fn decode(&self, start: usize, end: usize) -> Vec<u8> {
        self.sequence.decode(start, end)
    }

    pub fn contig<'a>(&'a self, index: usize) -> PackedSeqSlice<'a> {
        let start = self.starts[index];
        let len = self.starts[index + 1] - start;

        PackedSeqSlice::new(&self.sequence, start, len)
    }

    pub fn contig_start(&self, index: usize) -> usize {
        self.starts[index]
    }

    pub fn max_contig_len(&self) -> Option<usize> {
        (0..self.names.len()).map(|i| self.contig(i).len()).max()
    }

    pub fn total_length(&self) -> usize {
        *self.starts.last().unwrap()
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
