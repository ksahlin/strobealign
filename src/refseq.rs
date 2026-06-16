use crate::packed_seq::PackedSeq;

#[derive(Default, Debug, Clone)]
pub struct RefSequence {
    pub names: Vec<String>,
    sequences: Vec<PackedSeq>,
}

impl RefSequence {
    pub fn new() -> Self {
        RefSequence {
            names: vec![],
            sequences: vec![],
        }
    }

    pub fn push(&mut self, name: String, seq: PackedSeq) {
        self.names.push(name);
        self.sequences.push(seq);
    }

    pub fn contig(&self, index: usize) -> &PackedSeq {
        &self.sequences[index]
    }

    pub fn contig_len(&self, index: usize) -> usize {
        self.sequences[index].len()
    }

    pub fn max_contig_len(&self) -> Option<usize> {
        self.sequences.iter().map(|seq| seq.len()).max()
    }

    pub fn total_length(&self) -> usize {
        self.sequences.iter().map(|seq| seq.len()).sum()
    }
}
