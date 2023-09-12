/// Details about aligning a single or paired-end read
#[derive(Default,Debug,Clone)]
pub struct Details {
    /// find_nams_rescue() was needed
    pub nam_rescue: bool,

    // Number of NAMs found
    pub nams: usize,

    pub nam_inconsistent: usize,

    /// No. of times rescue by local alignment was attempted
    pub mate_rescue: usize,

    /// No. of computed alignments (get_alignment or rescue_mate)
    pub tried_alignment: usize,

    /// No. of gapped alignments computed (in get_alignment)
    pub gapped: usize,
}

impl Details {
    pub fn new() -> Self {
        Details::default()
    }
}
