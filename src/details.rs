use std::ops;

#[derive(Default, Debug, Clone)]
pub struct NamDetails {
    pub n_randstrobes: usize,
    pub n_nams: usize,
    pub n_rescue_nams: usize,
    pub nam_rescue: bool,
    /// Number of non-rescue hits
    pub n_hits: usize,
    pub n_rescue_hits: usize,

    pub time_randstrobes: f64,
    pub time_find_nams: f64,
    pub time_rescue: f64,
    pub time_sort_nams: f64,
}

/// Details about aligning a single read
#[derive(Default, Debug, Clone)]
pub struct Details {
    /// find_nams_rescue() was needed
    pub nam_rescue: bool,

    /// Number of non-rescue hits
    pub n_hits: usize,

    /// Number of rescue hits
    pub n_rescue_hits: usize,

    // Number of NAMs found
    pub nams: usize,

    pub nam_inconsistent: usize,

    /// No. of times rescue by local alignment was attempted
    pub mate_rescue: usize,

    /// No. of computed alignments (get_alignment or rescue_mate)
    pub tried_alignment: usize,

    /// No. of gapped alignments computed (in get_alignment)
    pub gapped: usize,

    /// No. of best alignments with same score
    pub best_alignments: usize,

}

impl ops::AddAssign<Details> for Details {
    fn add_assign(&mut self, rhs: Details) {
        self.nam_rescue = self.nam_rescue || rhs.nam_rescue;
        self.nams += rhs.nams;
        self.nam_inconsistent += rhs.nam_inconsistent;
        self.mate_rescue += rhs.mate_rescue;
        self.tried_alignment += rhs.tried_alignment;
        self.gapped += rhs.gapped;
        self.best_alignments += rhs.best_alignments;
    }    
} 

impl Details {
    pub fn new() -> Self {
        Details::default()
    }
}
