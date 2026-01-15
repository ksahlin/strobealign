use std::fmt::{Display, Formatter};
use std::time::Instant;

use fastrand::Rng;
use log::Level::Trace;
use log::trace;

use crate::chainer::Chainer;
use crate::details::NamDetails;
use crate::fasta::RefSequence;
use crate::index::StrobemerIndex;
use crate::mapper;
use crate::mcsstrategy::McsStrategy;
use crate::read::Read;

/// Non-overlapping approximate match
#[derive(Clone, Debug)]
pub struct Nam {
    pub nam_id: usize,
    pub ref_start: usize,
    pub ref_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub query_prev_match_startpos: usize,
    pub ref_prev_match_startpos: usize,
    pub n_matches: usize,
    pub ref_id: usize,
    pub score: f32,
    pub is_revcomp: bool,
}

impl Nam {
    pub fn ref_span(&self) -> usize {
        self.ref_end - self.ref_start
    }

    pub fn query_span(&self) -> usize {
        self.query_end - self.query_start
    }

    pub fn projected_ref_start(&self) -> usize {
        self.ref_start.saturating_sub(self.query_start)
    }
}

impl Display for Nam {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Nam(ref_id={}, query: {}..{}, ref: {}..{}, rc={}, score={})",
            self.ref_id,
            self.query_start,
            self.query_end,
            self.ref_start,
            self.ref_end,
            self.is_revcomp as u8,
            self.score
        )?;
        Ok(())
    }
}

/// Determine whether the NAM represents a match to the forward or
/// reverse-complemented sequence by checking in which orientation the
/// first and last strobe in the NAM match
///
/// - If first and last strobe match in forward orientation, return true.
/// - If first and last strobe match in reverse orientation, update the NAM
///   in place and return true.
/// - If first and last strobe do not match consistently, return false.
pub fn reverse_nam_if_needed(
    nam: &mut Nam,
    read: &Read,
    references: &[RefSequence],
    k: usize,
) -> bool {
    let ref_start_kmer = &references[nam.ref_id].sequence[nam.ref_start..nam.ref_start + k];
    let ref_end_kmer = &references[nam.ref_id].sequence[nam.ref_end - k..nam.ref_end];

    let (seq, seq_rc) = if nam.is_revcomp {
        (read.rc(), read.seq())
    } else {
        (read.seq(), read.rc())
    };
    let read_start_kmer = &seq[nam.query_start..nam.query_start + k];
    let read_end_kmer = &seq[nam.query_end - k..nam.query_end];
    if ref_start_kmer == read_start_kmer && ref_end_kmer == read_end_kmer {
        return true;
    }

    // False forward or false reverse (possible due to symmetrical hash values)
    // we need two extra checks for this - hopefully this will remove all the false matches we see
    // (true hash collisions should be very few)
    let read_len = read.len();
    let q_start_tmp = read_len - nam.query_end;
    let q_end_tmp = read_len - nam.query_start;
    // false reverse match, change coordinates in nam to forward
    let read_start_kmer = &seq_rc[q_start_tmp..q_start_tmp + k];
    let read_end_kmer = &seq_rc[q_end_tmp - k..q_end_tmp];
    if ref_start_kmer == read_start_kmer && ref_end_kmer == read_end_kmer {
        nam.is_revcomp = !nam.is_revcomp;
        nam.query_start = q_start_tmp;
        nam.query_end = q_end_tmp;
        true
    } else {
        false
    }
}

/// Obtain NAMs for a sequence record, doing rescue if needed.
///
/// NAMs are returned sorted by decreasing score
pub fn get_nams_by_chaining(
    sequence: &[u8],
    index: &StrobemerIndex,
    chainer: &Chainer,
    rescue_distance: usize,
    mcs_strategy: McsStrategy,
    rng: &mut Rng,
) -> (NamDetails, Vec<Nam>) {
    let timer = Instant::now();
    let query_randstrobes = mapper::randstrobes_query(sequence, &index.parameters);
    let time_randstrobes = timer.elapsed().as_secs_f64();

    trace!(
        "we have {} + {} randstrobes",
        query_randstrobes[0].len(),
        query_randstrobes[1].len()
    );
    let (mut nam_details, mut nams) =
        chainer.get_chains(&query_randstrobes, index, rescue_distance, mcs_strategy);

    let timer = Instant::now();

    nams.sort_by(|a, b| b.score.total_cmp(&a.score));
    shuffle_top_nams(&mut nams, rng);
    nam_details.time_sort_nams = timer.elapsed().as_secs_f64();
    nam_details.time_randstrobes = time_randstrobes;

    if log::log_enabled!(Trace) {
        trace!(
            "Found {} NAMs (rescue done: {})",
            nams.len(),
            nam_details.nam_rescue
        );
        let mut printed = 0;
        for nam in &nams {
            if nam.n_matches > 1 || printed < 10 {
                trace!("- {}", nam);
                printed += 1;
            }
        }
        if printed < nams.len() {
            trace!("+ {} single-anchor chains", nams.len() - printed);
        }
    }

    (nam_details, nams)
}

/// Shuffle the top-scoring NAMs. Input must be sorted by score.
/// This helps to ensure we pick a random location in case there are multiple
/// equally good ones.
fn shuffle_top_nams(nams: &mut [Nam], rng: &mut Rng) {
    if let Some(best) = nams.first() {
        let best_score = best.score;

        let pos = nams.iter().position(|nam| nam.score != best_score);
        let end = pos.unwrap_or(nams.len());
        if end > 1 {
            rng.shuffle(&mut nams[0..end]);
        }
    }
}
