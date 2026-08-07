//! Statistics about mapping a single or multiple reads

use std::ops;

use crate::hit::HitsDetails;

#[derive(Default, Debug, Clone)]
pub struct ChainDetails {
    // TODO should be moved out of here and into Details
    pub hits: HitsDetails,

    pub n_reads: usize,

    pub n_randstrobes: usize,

    pub n_anchors: usize,

    /// Number of chains found
    pub n_chains: usize,

    pub time_randstrobes: f64,
    pub time_find_hits: f64,
    pub time_chaining: f64,
    pub time_rescue: f64,
    pub time_recovery: f64,
    pub time_sort_chains: f64,

    /// Whether both orientations were tested
    pub both_orientations: bool,
}

impl ops::AddAssign<ChainDetails> for ChainDetails {
    fn add_assign(&mut self, rhs: ChainDetails) {
        self.hits += rhs.hits;
        self.n_reads += rhs.n_reads;
        self.n_randstrobes += rhs.n_randstrobes;
        self.n_anchors += rhs.n_anchors;
        self.n_chains += rhs.n_chains;
        self.time_randstrobes += rhs.time_randstrobes;
        self.time_find_hits += rhs.time_find_hits;
        self.time_chaining += rhs.time_chaining;
        self.time_rescue += rhs.time_rescue;
        self.time_recovery += rhs.time_recovery;
        self.time_sort_chains += rhs.time_sort_chains;
    }
}

/// Details about aligning a single read
#[derive(Default, Debug, Clone)]
pub struct Details {
    pub chain: ChainDetails,

    pub inconsistent_chains: usize,

    /// No. of times rescue by local alignment was attempted
    pub mate_rescue: usize,

    /// No. of computed alignments (get_alignment or rescue_mate)
    pub tried_alignment: usize,

    /// No. of gapped alignments computed (in get_alignment)
    pub gapped: usize,

    /// No. of best alignments with same score
    pub best_alignments: usize,

    pub time_extend: f64,
}

impl ops::AddAssign<Details> for Details {
    fn add_assign(&mut self, rhs: Details) {
        self.chain += rhs.chain;
        self.inconsistent_chains += rhs.inconsistent_chains;
        self.mate_rescue += rhs.mate_rescue;
        self.tried_alignment += rhs.tried_alignment;
        self.gapped += rhs.gapped;
        self.best_alignments += rhs.best_alignments;
        self.time_extend += rhs.time_extend;
    }
}

impl Details {
    pub fn new() -> Self {
        Details::default()
    }
}

impl From<ChainDetails> for Details {
    fn from(chain_details: ChainDetails) -> Self {
        Details {
            chain: chain_details,
            ..Details::default()
        }
    }
}
