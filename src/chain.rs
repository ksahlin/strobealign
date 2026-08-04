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

/// Obtain chains for a sequence record, doing rescue if needed.
pub fn get_chains(
    read: &Read,
    index: &StrobemerIndex,
    refseq: &RefSequence,
    chainer: &Chainer,
    rescue_distance: usize,
    mcs_strategy: McsStrategy,
) -> (ChainDetails, Vec<Chain>) {
    let timer = Instant::now();
    let query_randstrobes = randstrobes_query(read.seq(), &index.parameters);
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
        read,
        refseq,
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
