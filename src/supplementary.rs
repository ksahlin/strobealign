use crate::{nam::Nam, shuffle::shuffle_best};
use fastrand::Rng;

pub struct ChainGroup {
    pub chains: Vec<Nam>,
    pub mapq: u8,
}

pub fn forward_query_coords(nam: &Nam, read_length: usize) -> (usize, usize) {
    if nam.is_revcomp {
        (read_length - nam.query_end, read_length - nam.query_start)
    } else {
        (nam.query_start, nam.query_end)
    }
}

pub fn group_chains(chains: Vec<Nam>, read_length: usize, rng: &mut Rng) -> Vec<ChainGroup> {
    if chains.is_empty() {
        return vec![];
    }

    let mut chains = chains;
    chains.sort_unstable_by(|a, b| {
        let (a_start, a_end) = forward_query_coords(a, read_length);
        let (b_start, b_end) = forward_query_coords(b, read_length);
        a_start.cmp(&b_start).then(b_end.cmp(&a_end))
    });

    let mut groups: Vec<ChainGroup> = vec![];
    let mut current: Vec<Nam> = vec![];
    let mut group_start = 0usize;
    let mut group_end = 0usize;

    for nam in chains {
        let (s, e) = forward_query_coords(&nam, read_length);
        if current.is_empty() {
            group_start = s;
            group_end = e;
            current.push(nam);
        } else if s >= group_start && e <= group_end {
            // fully contained in the anchor's interval → same group
            current.push(nam);
        } else {
            // starts a new group
            groups.push(build_chain_group(&mut current, rng));
            group_start = s;
            group_end = e;
            current.push(nam);
        }
    }
    if !current.is_empty() {
        groups.push(build_chain_group(&mut current, rng));
    }

    groups
}

fn build_chain_group(group: &mut Vec<Nam>, rng: &mut Rng) -> ChainGroup {
    group.sort_unstable_by(|a, b| b.score.partial_cmp(&a.score).unwrap());
    let n_best = shuffle_best(group, |n| n.score, rng);

    let mapq = if group.len() == 1 {
        60
    } else if n_best > 1 {
        0
    } else {
        // exactly one best; use minimap2 formula with the runner-up
        let s1 = group[0].score;
        let s2 = group[1].score;
        let min_matches = group[0].anchors.len().min(10) as f32 / 10.0;
        let uncapped = 40.0 * (1.0 - s2 / s1) * min_matches * s1.ln();
        uncapped.min(60.0) as u8
    };

    ChainGroup {
        chains: std::mem::take(group),
        mapq,
    }
}

/// Computes the normalized query overlap between two NAMs.
pub fn query_overlap(a: &Nam, b: &Nam, read_length: usize) -> f32 {
    let (a_start, a_end) = forward_query_coords(a, read_length);
    let (b_start, b_end) = forward_query_coords(b, read_length);
    let start = a_start.max(b_start);
    let end = a_end.min(b_end);
    let overlap = end.saturating_sub(start);
    if overlap == 0 {
        return 0.0;
    }
    let len_a = a.query_span();
    let len_b = b.query_span();
    let min_len = len_a.min(len_b);
    overlap as f32 / min_len as f32
}

pub fn get_primary_chains(
    candidates: Vec<Nam>,
    read_length: usize,
    rng: &mut Rng,
) -> (Vec<(Nam, u8)>, Vec<Nam>) {
    if candidates.is_empty() {
        return (vec![], vec![]);
    }

    let mut groups = group_chains(candidates, read_length, rng);

    // sort groups by representative score, descending
    groups.sort_unstable_by(|a, b| b.chains[0].score.partial_cmp(&a.chains[0].score).unwrap());

    // secondaries are the non-representative chains of the best-scoring group only
    let secondary = groups[0].chains[1..].to_vec();

    // greedy primary-set construction over group representatives
    let mut primary: Vec<(Nam, u8)> = vec![];
    for mut group in groups {
        let rep = &group.chains[0];
        let qualifies = rep.anchors.len() >= 3
            && !primary
                .iter()
                .any(|(p, _)| query_overlap(rep, p, read_length) > 0.01);

        if qualifies {
            let mapq = group.mapq;
            primary.push((group.chains.remove(0), mapq));
        }
    }

    (primary, secondary)
}
