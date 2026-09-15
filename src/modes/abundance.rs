//! Abundance estimation mode (`--aemb`)

use crate::{
    chain::Chain,
    insertsize::InsertSizeDistribution,
    io::record::SequenceRecord,
    modes::map::{ChainPair, MappedChains, get_best_paired_mapping_location, get_chain_pairs},
    refseq::RefSequence,
};

/// Map a single-end read to the reference and estimate abundances
///
/// This implements abundance estimation mode (`--aemb`)
pub fn abundances_single_end_read(
    record: &SequenceRecord,
    refseq: &RefSequence,
    abundances: &mut [f64],
    chains: &[Chain],
) {
    let n_best = chains
        .iter()
        .take_while(|chain| chain.score == chains[0].score)
        .count();
    let weight = record.sequence.len() as f64 / n_best as f64;
    for chain in &chains[0..n_best] {
        let contig_id = refseq.unflatten(chain.ref_start).0;
        abundances[contig_id] += weight;
    }
}

/// Map a paired-end read pair to the reference and estimate abundances
///
/// This implements abundance estimation mode (`--aemb`)
pub fn abundances_paired_end_read(
    r1: &SequenceRecord,
    r2: &SequenceRecord,
    refseq: &RefSequence,
    abundances: &mut [f64],
    insert_size_distribution: &mut InsertSizeDistribution,
    chains_pair: &mut [Vec<Chain>; 2],
) {
    let [chains1, chains2] = chains_pair;
    if chains1.is_empty() && chains2.is_empty() {
        return;
    }

    let chain_pairs = get_chain_pairs(
        chains1,
        chains2,
        &refseq.starts,
        insert_size_distribution.mu,
        insert_size_distribution.sigma,
    );

    match get_best_paired_mapping_location(&chain_pairs, chains1, chains2, insert_size_distribution)
    {
        MappedChains::Individual(_, _) => {
            for (chains, read_len) in [(&chains1, r1.sequence.len()), (&chains2, r2.sequence.len())]
            {
                let n_best = chains
                    .iter()
                    .take_while(|chain| chain.score == chains[0].score)
                    .count();
                let weight = read_len as f64 / n_best as f64;
                for chain in &chains[0..n_best] {
                    let contig_id = refseq.unflatten(chain.ref_start).0;
                    abundances[contig_id] += weight;
                }
            }
        }
        MappedChains::Pair(_, _, joint_score) => {
            let n_best = chain_pairs
                .iter()
                .take_while(|chain_pair| chain_pair.score == joint_score)
                .count();
            let weight_r1 = r1.sequence.len() as f64 / n_best as f64;
            let weight_r2 = r2.sequence.len() as f64 / n_best as f64;
            for ChainPair {
                chain1,
                chain2,
                score: _,
            } in &chain_pairs[..n_best]
            {
                abundances[refseq.unflatten(chain1.ref_start).0] += weight_r1;
                abundances[refseq.unflatten(chain2.ref_start).0] += weight_r2;
            }
        }
    }
}
