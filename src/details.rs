//! Statistics about mapping a single or multiple reads

use std::ops;

use crate::hit::HitsDetails;

#[derive(Default, Debug, Clone)]
pub struct NamDetails {
    // TODO should be moved out of here and into Details
    pub hits: HitsDetails,

    pub n_reads: usize,

    pub n_randstrobes: usize,

    pub n_anchors: usize,

    /// Number of NAMs found
    pub n_nams: usize,

    /// Number of rescue NAMs found
    pub n_rescue_nams: usize,

    /// Number of times find_nams_rescue() was needed
    pub nam_rescue: usize,

    /// Number of partial hits found during rescue
    pub n_rescue_partial_hits: usize,

    /// Number of rescue hits
    pub n_rescue_hits: usize,

    pub time_randstrobes: f64,
    pub time_find_hits: f64,
    pub time_chaining: f64,
    pub time_rescue: f64,
    pub time_sort_nams: f64,
}

impl ops::AddAssign<NamDetails> for NamDetails {
    fn add_assign(&mut self, rhs: NamDetails) {
        self.hits += rhs.hits;
        self.n_reads += rhs.n_reads;
        self.n_randstrobes += rhs.n_randstrobes;
        self.n_anchors += rhs.n_anchors;
        self.n_nams += rhs.n_nams;
        self.n_rescue_nams += rhs.n_rescue_nams;
        self.nam_rescue += rhs.nam_rescue;
        self.n_rescue_partial_hits += rhs.n_rescue_partial_hits;
        self.n_rescue_hits += rhs.n_rescue_hits;
        self.time_randstrobes += rhs.time_randstrobes;
        self.time_find_hits += rhs.time_find_hits;
        self.time_chaining += rhs.time_chaining;
        self.time_rescue += rhs.time_rescue;
        self.time_sort_nams += rhs.time_sort_nams;
    }
}

/// Details about aligning a single read
#[derive(Default, Debug, Clone)]
pub struct Details {
    pub nam: NamDetails,

    pub inconsistent_nams: usize,

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
        self.nam += rhs.nam;
        self.inconsistent_nams += rhs.inconsistent_nams;
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

impl From<NamDetails> for Details {
    fn from(nam_details: NamDetails) -> Self {
        Details {
            nam: nam_details,
            .. Details::default()
        }
    }
}
