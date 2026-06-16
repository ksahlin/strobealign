use crate::packed_seq::PackedSeq;

#[derive(Default, Debug, Clone)]
pub struct RefSequence {
    pub names: Vec<String>,
    pub sequences: Vec<PackedSeq>,
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
}
