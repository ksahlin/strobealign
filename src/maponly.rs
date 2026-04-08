use fastrand::Rng;

use crate::chain::{Chain, get_chains};
use crate::chainer::Chainer;
use crate::details::Details;
use crate::index::StrobemerIndex;
use crate::insertsize::InsertSizeDistribution;
use crate::io::fasta::RefSequence;
use crate::io::paf::PafRecord;
use crate::io::record::{End, SequenceRecord};
use crate::mapper::{ChainPair, get_best_scoring_chain_pairs, mapping_quality};
use crate::mcsstrategy::McsStrategy;

/// Map a single-end read to the reference and return PAF records
///
/// This implements "mapping-only" mode in which no base-level alignments are computed
pub fn map_single_end_read(
    record: &SequenceRecord,
    index: &StrobemerIndex,
    references: &[RefSequence],
    rescue_distance: usize,
    mcs_strategy: McsStrategy,
    chainer: &Chainer,
    rng: &mut Rng,
) -> (Vec<PafRecord>, Details) {
    let (chain_details, chains) = get_chains(
        &record.sequence,
        index,
        chainer,
        rescue_distance,
        mcs_strategy,
        rng,
    );

    if chains.is_empty() {
        (vec![], chain_details.into())
    } else {
        let mapq = mapping_quality(&chains);
        (
            vec![paf_record_from_chain(
                &chains[0],
                &record.name,
                references,
                record.sequence.len(),
                Some(mapq),
                End::None,
            )],
            chain_details.into(),
        )
    }
}

/// Map a single-end read to the reference and estimate abundances
///
/// This implements abundance estimation mode (`--aemb`)
pub fn abundances_single_end_read(
    record: &SequenceRecord,
    index: &StrobemerIndex,
    abundances: &mut [f64],
    rescue_distance: usize,
    mcs_strategy: McsStrategy,
    chainer: &Chainer,
    rng: &mut Rng,
) {
    let (_, chains) = get_chains(
        &record.sequence,
        index,
        chainer,
        rescue_distance,
        mcs_strategy,
        rng,
    );
    let n_best = chains
        .iter()
        .take_while(|chain| chain.score == chains[0].score)
        .count();
    let weight = record.sequence.len() as f64 / n_best as f64;
    for chain in &chains[0..n_best] {
        abundances[chain.ref_id] += weight;
    }
}

/// Convert chain into PAF record
fn paf_record_from_chain(
    chain: &Chain,
    name: &str,
    references: &[RefSequence],
    query_length: usize,
    mapq: Option<u8>,
    end: End,
) -> PafRecord {
    PafRecord {
        query_name: name.into(),
        end,
        query_length: query_length as u64,
        query_start: chain.query_start as u64,
        query_end: chain.query_end as u64,
        is_revcomp: chain.is_revcomp,
        target_name: references[chain.ref_id].name.clone(),
        target_length: references[chain.ref_id].sequence.len() as u64,
        target_start: chain.ref_start as u64,
        target_end: chain.ref_end as u64,
        matching_bases: chain.matching_bases as u64,
        alignment_length: (chain.ref_end - chain.ref_start) as u64,
        mapping_quality: mapq,
    }
}

/// Map a paired-end read pair to the reference and return PAF records
///
/// This implements "mapping-only" mode in which no base-level alignments are computed
pub fn map_paired_end_read(
    r1: &SequenceRecord,
    r2: &SequenceRecord,
    index: &StrobemerIndex,
    references: &[RefSequence],
    rescue_distance: usize,
    insert_size_distribution: &mut InsertSizeDistribution,
    mcs_strategy: McsStrategy,
    chainer: &Chainer,
    rng: &mut Rng,
) -> (Vec<PafRecord>, Details) {
    let (mut chain_details1, chains1) = get_chains(
        &r1.sequence,
        index,
        chainer,
        rescue_distance,
        mcs_strategy,
        rng,
    );
    let (chain_details2, chains2) = get_chains(
        &r2.sequence,
        index,
        chainer,
        rescue_distance,
        mcs_strategy,
        rng,
    );

    let chain_pairs = get_best_scoring_chain_pairs(
        &chains1,
        &chains2,
        insert_size_distribution.mu,
        insert_size_distribution.sigma,
    );
    let mapped_chain =
        get_best_paired_map_location(&chain_pairs, &chains1, &chains2, insert_size_distribution);
    let mut records = vec![];

    match mapped_chain {
        MappedChains::Individual(chain1, chain2) => {
            if let Some(chain) = chain1 {
                records.push(paf_record_from_chain(
                    &chain,
                    &r1.name,
                    references,
                    r1.sequence.len(),
                    None,
                    End::One,
                ))
            }
            if let Some(chain) = chain2 {
                records.push(paf_record_from_chain(
                    &chain,
                    &r2.name,
                    references,
                    r2.sequence.len(),
                    None,
                    End::Two,
                ))
            }
        }
        MappedChains::Pair(chain1, chain2) => {
            records.push(paf_record_from_chain(
                &chain1,
                &r1.name,
                references,
                r1.sequence.len(),
                None,
                End::One,
            ));
            records.push(paf_record_from_chain(
                &chain2,
                &r2.name,
                references,
                r2.sequence.len(),
                None,
                End::Two,
            ));
        }
        MappedChains::Unmapped => {}
    }
    chain_details1 += chain_details2;
    (records, chain_details1.into())
}

/// Map a paired-end read pair to the reference and estimate abundances
///
/// This implements abundance estimation mode (`--aemb`)
pub fn abundances_paired_end_read(
    r1: &SequenceRecord,
    r2: &SequenceRecord,
    index: &StrobemerIndex,
    abundances: &mut [f64],
    rescue_distance: usize,
    insert_size_distribution: &mut InsertSizeDistribution,
    mcs_strategy: McsStrategy,
    chainer: &Chainer,
    rng: &mut Rng,
) {
    let chains1 = get_chains(
        &r1.sequence,
        index,
        chainer,
        rescue_distance,
        mcs_strategy,
        rng,
    )
    .1;
    let chains2 = get_chains(
        &r2.sequence,
        index,
        chainer,
        rescue_distance,
        mcs_strategy,
        rng,
    )
    .1;

    let chain_pairs = get_best_scoring_chain_pairs(
        &chains1,
        &chains2,
        insert_size_distribution.mu,
        insert_size_distribution.sigma,
    );
    let mapped_chain =
        get_best_paired_map_location(&chain_pairs, &chains1, &chains2, insert_size_distribution);

    match mapped_chain {
        MappedChains::Pair(chain1, chain2) => {
            let joint_score = chain1.score + chain2.score;
            let n_best = chain_pairs
                .iter()
                .take_while(|chain_pair| {
                    chain_pair.chain1.as_ref().map_or(0.0, |chain| chain.score)
                        + chain_pair.chain2.as_ref().map_or(0.0, |chain| chain.score)
                        == joint_score
                })
                .count();
            let weight_r1 = r1.sequence.len() as f64 / n_best as f64;
            let weight_r2 = r2.sequence.len() as f64 / n_best as f64;
            for chain_pair in &chain_pairs[..n_best] {
                if let Some(chain) = &chain_pair.chain1 {
                    abundances[chain.ref_id] += weight_r1;
                }
                if let Some(chain) = &chain_pair.chain2 {
                    abundances[chain.ref_id] += weight_r2;
                }
            }
        }
        MappedChains::Individual(_, _) => {
            for (chains, read_len) in [(&chains1, r1.sequence.len()), (&chains2, r2.sequence.len())]
            {
                let n_best = chains
                    .iter()
                    .take_while(|chain| chain.score == chains[0].score)
                    .count();
                let weight = read_len as f64 / n_best as f64;
                for chain in &chains[0..n_best] {
                    abundances[chain.ref_id] += weight;
                }
            }
        }
        MappedChains::Unmapped => {}
    }
}

enum MappedChains {
    Individual(Option<Chain>, Option<Chain>),
    Pair(Chain, Chain),
    Unmapped,
}

/// Given two lists of chains from R1 and R2, find the best location (preferably a proper pair).
/// This is used for mapping-only (PAF) mode and abundances output
fn get_best_paired_map_location(
    chain_pairs: &[ChainPair],
    chains1: &[Chain],
    chains2: &[Chain],
    insert_size_distribution: &mut InsertSizeDistribution,
) -> MappedChains {
    if chain_pairs.is_empty() && chains1.is_empty() && chains2.is_empty() {
        return MappedChains::Unmapped;
    }

    // Find first chain pair that is a proper pair.
    // The first one is also the one with the highest score
    // since chain_pairs is sorted descending by score
    let best_joint_pair = chain_pairs
        .iter()
        .find(|&chain_pair| chain_pair.chain1.is_some() && chain_pair.chain2.is_some());

    let joint_score = if let Some(chain_pair) = best_joint_pair {
        chain_pair.chain1.as_ref().map_or(0.0, |chain| chain.score)
            + chain_pair.chain2.as_ref().map_or(0.0, |chain| chain.score)
    } else {
        0.0
    };

    // Get individual best scores.
    // chains1 and chains2 are also sorted descending by score.
    let best_individual_chain1 = chains1.first();
    let best_individual_chain2 = chains2.first();

    let individual_score = best_individual_chain1.map_or(0.0, |chain| chain.score)
        + best_individual_chain2.map_or(0.0, |chain| chain.score);

    // Divisor 2 is penalty for being mapped individually
    if joint_score > individual_score / 2.0 {
        let best_joint_pair = best_joint_pair.unwrap();
        let best = (
            best_joint_pair.chain1.clone(),
            best_joint_pair.chain2.clone(),
        );
        if insert_size_distribution.sample_size < 400 {
            insert_size_distribution.update(
                best.0
                    .as_ref()
                    .unwrap()
                    .ref_start
                    .abs_diff(best.1.as_ref().unwrap().ref_start),
            );
        }

        MappedChains::Pair(best.0.unwrap(), best.1.unwrap())
    } else {
        MappedChains::Individual(
            best_individual_chain1.cloned(),
            best_individual_chain2.cloned(),
        )
    }
}
