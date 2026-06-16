use std::fmt::{Display, Formatter};
use std::time::Instant;

use fastrand::Rng;
use log::Level::Trace;
use log::trace;

use crate::chainer::Anchor;
use crate::chainer::Chainer;
use crate::details::ChainDetails;
use crate::index::StrobemerIndex;
use crate::mcsstrategy::McsStrategy;
use crate::read::Read;
use crate::refseq::RefSequence;
use crate::seeding::randstrobes_query;
use crate::shuffle::shuffle_best;

/// A list of anchors
#[derive(Clone, Debug, Default)]
pub struct Chain {
    pub id: usize,
    /// Start coordinate of the contig containing `ref_start`.
    pub ref_contig_start: usize,
    pub ref_start: usize,
    pub ref_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub matching_bases: usize,
    pub score: f32,
    pub is_revcomp: bool,
    pub anchors: Vec<Anchor>,
}

impl Chain {
    pub fn ref_span(&self) -> usize {
        self.ref_end - self.ref_start
    }

    pub fn query_span(&self) -> usize {
        self.query_end - self.query_start
    }

    pub fn projected_ref_start(&self) -> usize {
        self.ref_start
            .saturating_sub(self.query_start)
            .max(self.ref_contig_start)
    }

    /// Returns whether a chain represents a consistent match between read and
    /// reference by comparing the nucleotide sequences of the first and last
    /// strobe (taking orientation into account).
    pub fn is_consistent(&self, read: &Read, refseq: &RefSequence, k: usize) -> bool {
        let ref_start_kmer = refseq.decode(self.ref_start, self.ref_start + k);
        let ref_end_kmer = refseq.decode(self.ref_end - k, self.ref_end);

        let seq = if self.is_revcomp {
            read.rc()
        } else {
            read.seq()
        };
        let read_start_kmer = &seq[self.query_start..self.query_start + k];
        let read_end_kmer = &seq[self.query_end - k..self.query_end];

        ref_start_kmer == read_start_kmer && ref_end_kmer == read_end_kmer
    }
}

impl Display for Chain {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Chain(query: {}..{}, ref: {}..{}, rc={}, score={})",
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

/// Determine whether the chain represents a match to the forward or
/// reverse-complemented sequence by checking in which orientation the
/// first and last strobe in the chain match
///
/// - If first and last strobe match in forward orientation, return true.
/// - If first and last strobe match in reverse orientation, update the chain
///   in place and return true.
/// - If first and last strobe do not match consistently, return false.
pub fn reverse_chain_if_needed(
    chain: &mut Chain,
    read: &Read,
    refseq: &RefSequence,
    k: usize,
) -> bool {
    let ref_start_kmer = refseq.decode(chain.ref_start, chain.ref_start + k);
    let ref_end_kmer = refseq.decode(chain.ref_end - k, chain.ref_end);

    let (seq, seq_rc) = if chain.is_revcomp {
        (read.rc(), read.seq())
    } else {
        (read.seq(), read.rc())
    };
    let read_start_kmer = &seq[chain.query_start..chain.query_start + k];
    let read_end_kmer = &seq[chain.query_end - k..chain.query_end];
    if ref_start_kmer == read_start_kmer && ref_end_kmer == read_end_kmer {
        return true;
    }

    // False forward or false reverse (possible due to symmetrical hash values)
    // we need two extra checks for this - hopefully this will remove all the false matches we see
    // (true hash collisions should be very few)
    let read_len = read.len();
    let q_start_tmp = read_len - chain.query_end;
    let q_end_tmp = read_len - chain.query_start;
    // false reverse match, change coordinates in chain to forward
    let read_start_kmer = &seq_rc[q_start_tmp..q_start_tmp + k];
    let read_end_kmer = &seq_rc[q_end_tmp - k..q_end_tmp];
    if ref_start_kmer == read_start_kmer && ref_end_kmer == read_end_kmer {
        chain.is_revcomp = !chain.is_revcomp;
        chain.query_start = q_start_tmp;
        chain.query_end = q_end_tmp;
        true
    } else {
        false
    }
}

/// Obtain chains for a sequence record, doing rescue if needed.
pub fn get_chains(
    sequence: &[u8],
    index: &StrobemerIndex,
    chainer: &Chainer,
    rescue_distance: usize,
    mcs_strategy: McsStrategy,
) -> (ChainDetails, Vec<Chain>) {
    let timer = Instant::now();
    let query_randstrobes = randstrobes_query(sequence, &index.parameters);
    let time_randstrobes = timer.elapsed().as_secs_f64();

    trace!(
        "we have {} + {} randstrobes",
        query_randstrobes[0].len(),
        query_randstrobes[1].len()
    );

    let (mut chain_details, chains) = chainer.get_chains(
        &query_randstrobes,
        index,
        rescue_distance,
        mcs_strategy,
        sequence.len(),
    );

    chain_details.time_randstrobes = time_randstrobes;

    (chain_details, chains)
}

pub fn sort_chains(chains: &mut [Chain], rng: &mut Rng) -> f64 {
    let timer = Instant::now();
    chains.sort_by(|a, b| b.score.total_cmp(&a.score));
    shuffle_best(chains, |chain| chain.score, rng);

    if log::log_enabled!(Trace) {
        trace!("Found {} chains", chains.len());
        let mut printed = 0;
        for chain in chains.iter() {
            if chain.anchors.len() > 1 || printed < 10 {
                trace!("- {}", chain);
                printed += 1;
            }
        }
        if printed < chains.len() {
            trace!("+ {} single-anchor chains", chains.len() - printed);
        }
    }

    timer.elapsed().as_secs_f64()
}
