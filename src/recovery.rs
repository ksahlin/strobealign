use crate::chain::Chain;
use crate::chainer::{Anchor, Chainer, filtered_hits_to_anchors};
use crate::hit::Hit;
use crate::index::StrobemerIndex;

const MAX_BAND: usize = 1000;

impl Chainer {
    pub fn recover_anchors(
        &self,
        chains: &mut [Chain],
        hits: &[Vec<Hit>; 2],
        index: &StrobemerIndex,
        read_len: usize,
    ) {
        if chains.is_empty() {
            return;
        }

        let max_dist = index.parameters.randstrobe.max_dist as usize;
        let mut anchors: [Vec<Anchor>; 2] = [vec![], vec![]];
        for is_revcomp in 0..2 {
            let bands = chain_bands(chains, is_revcomp == 1, read_len);
            if bands.is_empty() {
                continue;
            }
            anchors[is_revcomp] =
                filtered_hits_to_anchors(&hits[is_revcomp], index, &bands, max_dist);
            anchors[is_revcomp].sort_unstable_by_key(|a| (a.ref_id, a.ref_start, a.query_start));
            anchors[is_revcomp].dedup();
        }

        for chain in chains {
            let anchors = &anchors[chain.is_revcomp as usize];
            self.recover_into_chain(chain, anchors, index.k(), read_len);
        }
    }

    fn recover_into_chain(&self, chain: &mut Chain, anchors: &[Anchor], k: usize, read_len: usize) {
        let candidates = anchors_near(anchors, chain, read_len);
        if candidates.is_empty() {
            return;
        }

        let mut recovered = Vec::with_capacity(chain.anchors.len() + candidates.len());
        let mut candidate = 0;
        let mut gained = 0f32;

        for anchor in chain.anchors.iter().rev() {
            let gap = candidate;
            while candidate < candidates.len() && candidates[candidate].ref_start < anchor.ref_start
            {
                candidate += 1;
            }
            gained += self.take_best_run(
                &candidates[gap..candidate],
                Some(anchor),
                &mut recovered,
                read_len,
            );
            recovered.push(*anchor);
        }
        gained += self.take_best_run(&candidates[candidate..], None, &mut recovered, read_len);

        if recovered.len() == chain.anchors.len() {
            return;
        }

        recovered.reverse();
        chain.score += gained;
        chain.matching_bases = matching_bases(&recovered, k);
        chain.ref_start = recovered.last().unwrap().ref_start;
        chain.query_start = recovered.last().unwrap().query_start;
        chain.ref_end = recovered[0].ref_start + k;
        chain.query_end = recovered[0].query_start + k;
        chain.anchors = recovered;
    }

    fn link_score(&self, from: &Anchor, to: &Anchor, read_len: usize) -> Option<f32> {
        let dq = to.query_start.checked_sub(from.query_start)?;
        let dr = to.ref_start.checked_sub(from.ref_start)?;
        if dq == 0 || dr == 0 || dr >= self.parameters.max_ref_gap.unwrap_or(read_len) {
            return None;
        }
        if dq.max(dr) as f32 / dq.min(dr) as f32 > self.parameters.max_diagonal_ratio {
            return None;
        }

        Some(self.compute_score_cached(dq, dr))
    }

    fn take_best_run(
        &self,
        gap: &[Anchor],
        next: Option<&Anchor>,
        recovered: &mut Vec<Anchor>,
        read_len: usize,
    ) -> f32 {
        if gap.is_empty() {
            return 0f32;
        }
        const MAX_LOOKBACK: usize = 32;
        let previous = recovered.last().copied();

        let mut replaced_link = 0f32;
        if let (Some(previous), Some(next)) = (previous, next) {
            replaced_link = self.compute_score_cached(
                next.query_start - previous.query_start,
                next.ref_start - previous.ref_start,
            );
        }

        let mut best_ending_at = vec![f32::NEG_INFINITY; gap.len()];
        let mut predecessor = vec![usize::MAX; gap.len()];
        for i in 0..gap.len() {
            let opening_link = if let Some(previous) = previous {
                self.link_score(&previous, &gap[i], read_len)
                    .unwrap_or(f32::NEG_INFINITY)
            } else {
                0f32
            };
            best_ending_at[i] = opening_link + self.parameters.matches_weight;

            for j in i.saturating_sub(MAX_LOOKBACK)..i {
                if best_ending_at[j] == f32::NEG_INFINITY {
                    continue;
                }
                let Some(link) = self.link_score(&gap[j], &gap[i], read_len) else {
                    continue;
                };
                let score = best_ending_at[j] + link + self.parameters.matches_weight;
                if score > best_ending_at[i] {
                    best_ending_at[i] = score;
                    predecessor[i] = j;
                }
            }
        }

        let mut best_total = replaced_link;
        let mut best_end = usize::MAX;
        for i in 0..gap.len() {
            if best_ending_at[i] == f32::NEG_INFINITY {
                continue;
            }
            let mut total = best_ending_at[i];
            if let Some(next) = next {
                let Some(closing_link) = self.link_score(&gap[i], next, read_len) else {
                    continue;
                };
                total += closing_link;
            }
            if total > best_total {
                best_total = total;
                best_end = i;
            }
        }
        if best_end == usize::MAX {
            return 0f32;
        }

        let start = recovered.len();
        let mut anchor = best_end;
        while anchor != usize::MAX {
            recovered.push(gap[anchor]);
            anchor = predecessor[anchor];
        }
        recovered[start..].reverse();

        best_total - replaced_link
    }
}

/// Returns the diagonals each chain spans, widened by the band and sorted.
fn chain_bands(chains: &[Chain], is_revcomp: bool, read_len: usize) -> Vec<(i64, i64)> {
    let band = (read_len / 2).min(MAX_BAND) as i64;
    let mut bands = vec![];
    for chain in chains {
        if chain.is_revcomp == is_revcomp {
            let start_diagonal = chain.ref_start as i64 - chain.query_start as i64;
            let end_diagonal = chain.ref_end as i64 - chain.query_end as i64;
            bands.push((
                start_diagonal.min(end_diagonal) - band,
                start_diagonal.max(end_diagonal) + band,
            ));
        }
    }
    bands.sort_unstable();

    bands
}

/// Returns the set of filtered anchors that are near the chain and the
/// diagonals it spans
fn anchors_near(anchors: &[Anchor], chain: &Chain, read_len: usize) -> Vec<Anchor> {
    let reach = read_len;
    let first = chain.ref_start.saturating_sub(reach);
    let last = chain.ref_end + reach;
    let start = anchors.partition_point(|a| (a.ref_id, a.ref_start) < (chain.ref_id, first));
    let end = anchors.partition_point(|a| (a.ref_id, a.ref_start) <= (chain.ref_id, last));

    // dynamic diagonal band length
    let band = (read_len / 2).min(MAX_BAND) as i64;
    let start_diagonal = chain.ref_start as i64 - chain.query_start as i64;
    let end_diagonal = chain.ref_end as i64 - chain.query_end as i64;
    let lowest = start_diagonal.min(end_diagonal) - band;
    let highest = start_diagonal.max(end_diagonal) + band;

    let mut near = vec![];
    for anchor in &anchors[start..end] {
        let diagonal = anchor.ref_start as i64 - anchor.query_start as i64;
        if diagonal >= lowest && diagonal <= highest {
            near.push(*anchor);
        }
    }

    near
}

fn matching_bases(anchors: &[Anchor], k: usize) -> usize {
    let mut matching_bases = k;
    let mut ref_coverage = anchors[0].ref_start;
    for anchor in &anchors[1..] {
        matching_bases += ref_coverage.saturating_sub(anchor.ref_start).min(k);
        ref_coverage = anchor.ref_start;
    }

    matching_bases
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::chainer::ChainingParameters;

    #[test]
    fn an_anchor_fills_a_gap_in_a_chain() {
        let chainer = Chainer::new(20, ChainingParameters::default());
        #[rustfmt::skip]
        let mut chain = Chain {
            id: 3,
            ref_start: 0,
            ref_end: 140,
            query_start: 0,
            query_end: 140,
            matching_bases: 20 * 4,
            ref_id: 0,
            score: 0.0,
            is_revcomp: false,
            anchors: vec![
                Anchor { ref_id: 0, ref_start: 120, query_start: 120 },
                Anchor { ref_id: 0, ref_start:  60, query_start:  60 },
                Anchor { ref_id: 0, ref_start:  30, query_start:  30 },
                Anchor { ref_id: 0, ref_start:   0, query_start:   0 },
            ],
        };
        #[rustfmt::skip]
        let missing = [Anchor { ref_id: 0, ref_start: 90, query_start: 90 }];

        chainer.recover_into_chain(&mut chain, &missing, 20, 200);

        #[rustfmt::skip]
        assert_eq!(chain.anchors, vec![
            Anchor { ref_id: 0, ref_start: 120, query_start: 120 },
            Anchor { ref_id: 0, ref_start:  90, query_start:  90 },
            Anchor { ref_id: 0, ref_start:  60, query_start:  60 },
            Anchor { ref_id: 0, ref_start:  30, query_start:  30 },
            Anchor { ref_id: 0, ref_start:   0, query_start:   0 },
        ]);
        assert_eq!(chain.matching_bases, 20 * 5);
        assert_eq!((chain.ref_start, chain.ref_end), (0, 140));
        assert_eq!((chain.query_start, chain.query_end), (0, 140));
        assert_eq!(chain.id, 3);
    }

    #[test]
    fn each_end_of_a_chain_takes_more_than_one_anchor() {
        let chainer = Chainer::new(20, ChainingParameters::default());
        #[rustfmt::skip]
        let mut chain = Chain {
            id: 0,
            ref_start: 100,
            ref_end: 150,
            query_start: 100,
            query_end: 150,
            matching_bases: 20 * 2,
            ref_id: 0,
            score: 39.52,
            is_revcomp: false,
            anchors: vec![
                Anchor { ref_id: 0, ref_start: 130, query_start: 130 },
                Anchor { ref_id: 0, ref_start: 100, query_start: 100 },
            ],
        };

        #[rustfmt::skip]
        let outside = [
            Anchor { ref_id: 0, ref_start:  40, query_start:  40 },
            Anchor { ref_id: 0, ref_start:  70, query_start:  70 },
            Anchor { ref_id: 0, ref_start: 160, query_start: 160 },
            Anchor { ref_id: 0, ref_start: 190, query_start: 190 },
        ];

        chainer.recover_into_chain(&mut chain, &outside, 20, 200);

        #[rustfmt::skip]
        assert_eq!(chain.anchors, vec![
            Anchor { ref_id: 0, ref_start: 190, query_start: 190 },
            Anchor { ref_id: 0, ref_start: 160, query_start: 160 },
            Anchor { ref_id: 0, ref_start: 130, query_start: 130 },
            Anchor { ref_id: 0, ref_start: 100, query_start: 100 },
            Anchor { ref_id: 0, ref_start:  70, query_start:  70 },
            Anchor { ref_id: 0, ref_start:  40, query_start:  40 },
        ]);
        assert_eq!(chain.matching_bases, 20 * 6);
        assert_eq!((chain.ref_start, chain.ref_end), (40, 210));
        assert_eq!((chain.query_start, chain.query_end), (40, 210));
        assert!((chain.score - (20.0 + 5.0 * 19.5 + 6.0 * 0.01)).abs() < 0.0001);
    }

    #[test]
    fn the_best_of_the_candidates_for_a_gap_is_the_one_taken() {
        let chainer = Chainer::new(20, ChainingParameters::default());
        #[rustfmt::skip]
        let mut chain = Chain {
            id: 0,
            ref_start: 0,
            ref_end: 220,
            query_start: 0,
            query_end: 220,
            matching_bases: 20 * 2,
            ref_id: 0,
            score: 31.02,
            is_revcomp: false,
            anchors: vec![
                Anchor { ref_id: 0, ref_start: 200, query_start: 200 },
                Anchor { ref_id: 0, ref_start:   0, query_start:   0 },
            ],
        };

        #[rustfmt::skip]
        let competing = [
            Anchor { ref_id: 0, ref_start: 100, query_start: 150 },
            Anchor { ref_id: 0, ref_start: 110, query_start: 110 },
        ];

        chainer.recover_into_chain(&mut chain, &competing, 20, 400);

        #[rustfmt::skip]
        assert_eq!(chain.anchors, vec![
            Anchor { ref_id: 0, ref_start: 200, query_start: 200 },
            Anchor { ref_id: 0, ref_start: 110, query_start: 110 },
            Anchor { ref_id: 0, ref_start:   0, query_start:   0 },
        ]);
        assert!((chain.score - (31.02 + 15.5 + 16.5 - 11.0 + 0.01)).abs() < 0.0001);
    }

    #[test]
    fn a_gap_takes_every_anchor_that_raises_the_score() {
        let chainer = Chainer::new(20, ChainingParameters::default());
        #[rustfmt::skip]
        let mut chain = Chain {
            id: 0,
            ref_start: 0,
            ref_end: 220,
            query_start: 0,
            query_end: 220,
            matching_bases: 20 * 2,
            ref_id: 0,
            score: 31.02,
            is_revcomp: false,
            anchors: vec![
                Anchor { ref_id: 0, ref_start: 200, query_start: 200 },
                Anchor { ref_id: 0, ref_start:   0, query_start:   0 },
            ],
        };

        #[rustfmt::skip]
        let missing = [
            Anchor { ref_id: 0, ref_start:  50, query_start:  50 },
            Anchor { ref_id: 0, ref_start: 100, query_start: 100 },
            Anchor { ref_id: 0, ref_start: 150, query_start: 150 },
        ];

        chainer.recover_into_chain(&mut chain, &missing, 20, 400);

        #[rustfmt::skip]
        assert_eq!(chain.anchors, vec![
            Anchor { ref_id: 0, ref_start: 200, query_start: 200 },
            Anchor { ref_id: 0, ref_start: 150, query_start: 150 },
            Anchor { ref_id: 0, ref_start: 100, query_start: 100 },
            Anchor { ref_id: 0, ref_start:  50, query_start:  50 },
            Anchor { ref_id: 0, ref_start:   0, query_start:   0 },
        ]);
        assert_eq!(chain.matching_bases, 20 * 5);
        assert!((chain.score - (20.0 + 4.0 * 18.5 + 5.0 * 0.01)).abs() < 0.0001);
    }

    #[test]
    fn an_anchor_that_pays_gives_way_to_a_better_one_behind_it() {
        let chainer = Chainer::new(20, ChainingParameters::default());
        #[rustfmt::skip]
        let mut chain = Chain {
            id: 0,
            ref_start: 0,
            ref_end: 220,
            query_start: 0,
            query_end: 220,
            matching_bases: 20 * 2,
            ref_id: 0,
            score: 100.0,
            is_revcomp: false,
            anchors: vec![
                Anchor { ref_id: 0, ref_start: 200, query_start: 200 },
                Anchor { ref_id: 0, ref_start:   0, query_start:   0 },
            ],
        };

        #[rustfmt::skip]
        let gap = [
            Anchor { ref_id: 0, ref_start: 20, query_start: 110 },
            Anchor { ref_id: 0, ref_start: 30, query_start:  30 },
        ];

        chainer.recover_into_chain(&mut chain, &gap, 20, 400);

        #[rustfmt::skip]
        assert_eq!(chain.anchors, vec![
            Anchor { ref_id: 0, ref_start: 200, query_start: 200 },
            Anchor { ref_id: 0, ref_start:  30, query_start:  30 },
            Anchor { ref_id: 0, ref_start:   0, query_start:   0 },
        ]);
        assert!((chain.score - 121.01).abs() < 0.0001);
    }

    #[test]
    fn an_anchor_off_the_diagonal_is_left_out() {
        let chainer = Chainer::new(20, ChainingParameters::default());
        #[rustfmt::skip]
        let mut chain = Chain {
            id: 0,
            ref_start: 0,
            ref_end: 140,
            query_start: 0,
            query_end: 140,
            matching_bases: 20 * 5,
            ref_id: 0,
            score: 50.0,
            is_revcomp: false,
            anchors: vec![
                Anchor { ref_id: 0, ref_start: 120, query_start: 120 },
                Anchor { ref_id: 0, ref_start:  90, query_start:  90 },
                Anchor { ref_id: 0, ref_start:  60, query_start:  60 },
                Anchor { ref_id: 0, ref_start:  30, query_start:  30 },
                Anchor { ref_id: 0, ref_start:   0, query_start:   0 },
            ],
        };

        #[rustfmt::skip]
        let off_diagonal = [Anchor { ref_id: 0, ref_start: 75, query_start: 62 }];

        chainer.recover_into_chain(&mut chain, &off_diagonal, 20, 200);

        assert_eq!(chain.anchors.len(), 5);
        assert_eq!(chain.matching_bases, 20 * 5);
        assert_eq!(chain.score, 50.0);
    }

    #[test]
    fn a_run_drifting_off_the_diagonal_is_cut_where_it_leaves_the_band() {
        let chainer = Chainer::new(20, ChainingParameters::default());
        #[rustfmt::skip]
        let mut chain = Chain {
            id: 0,
            ref_start: 0,
            ref_end: 120,
            query_start: 0,
            query_end: 120,
            matching_bases: 20 * 2,
            ref_id: 0,
            score: 50.0,
            is_revcomp: false,
            anchors: vec![
                Anchor { ref_id: 0, ref_start: 100, query_start: 100 },
                Anchor { ref_id: 0, ref_start:   0, query_start:   0 },
            ],
        };
        #[rustfmt::skip]
        let drifting = [
            Anchor { ref_id: 0, ref_start: 220, query_start: 120 },
            Anchor { ref_id: 0, ref_start: 340, query_start: 140 },
            Anchor { ref_id: 0, ref_start: 460, query_start: 160 },
            Anchor { ref_id: 0, ref_start: 580, query_start: 180 },
            Anchor { ref_id: 0, ref_start: 700, query_start: 200 },
        ];

        chainer.recover_into_chain(&mut chain, &drifting, 20, 600);

        #[rustfmt::skip]
        assert_eq!(chain.anchors, vec![
            Anchor { ref_id: 0, ref_start: 460, query_start: 160 },
            Anchor { ref_id: 0, ref_start: 340, query_start: 140 },
            Anchor { ref_id: 0, ref_start: 220, query_start: 120 },
            Anchor { ref_id: 0, ref_start: 100, query_start: 100 },
            Anchor { ref_id: 0, ref_start:   0, query_start:   0 },
        ]);
        assert!((chain.score - 70.0427).abs() < 0.001);
    }

    #[test]
    fn an_anchor_going_backwards_in_the_query_is_left_out() {
        let chainer = Chainer::new(20, ChainingParameters::default());
        #[rustfmt::skip]
        let mut chain = Chain {
            id: 0,
            ref_start: 0,
            ref_end: 80,
            query_start: 0,
            query_end: 80,
            matching_bases: 20 * 3,
            ref_id: 0,
            score: 50.0,
            is_revcomp: false,
            anchors: vec![
                Anchor { ref_id: 0, ref_start: 60, query_start: 60 },
                Anchor { ref_id: 0, ref_start: 30, query_start: 30 },
                Anchor { ref_id: 0, ref_start:  0, query_start:  0 },
            ],
        };
        #[rustfmt::skip]
        let backwards = [Anchor { ref_id: 0, ref_start: 45, query_start: 10 }];

        chainer.recover_into_chain(&mut chain, &backwards, 20, 200);

        assert_eq!(chain.anchors.len(), 3);
        assert_eq!(chain.score, 50.0);
    }

    #[test]
    fn an_anchor_the_chain_already_has_is_left_out() {
        let chainer = Chainer::new(20, ChainingParameters::default());
        #[rustfmt::skip]
        let mut chain = Chain {
            id: 0,
            ref_start: 0,
            ref_end: 80,
            query_start: 0,
            query_end: 80,
            matching_bases: 20 * 3,
            ref_id: 0,
            score: 50.0,
            is_revcomp: false,
            anchors: vec![
                Anchor { ref_id: 0, ref_start: 60, query_start: 60 },
                Anchor { ref_id: 0, ref_start: 30, query_start: 30 },
                Anchor { ref_id: 0, ref_start:  0, query_start:  0 },
            ],
        };
        #[rustfmt::skip]
        let known = [
            Anchor { ref_id: 0, ref_start:  0, query_start:  0 },
            Anchor { ref_id: 0, ref_start: 30, query_start: 30 },
            Anchor { ref_id: 0, ref_start: 60, query_start: 60 },
        ];

        chainer.recover_into_chain(&mut chain, &known, 20, 200);

        assert_eq!(chain.anchors.len(), 3);
        assert_eq!(chain.matching_bases, 20 * 3);
        assert_eq!(chain.score, 50.0);
    }

    #[test]
    fn an_anchor_on_another_contig_is_left_out() {
        let chainer = Chainer::new(20, ChainingParameters::default());
        #[rustfmt::skip]
        let mut chain = Chain {
            id: 0,
            ref_start: 0,
            ref_end: 80,
            query_start: 0,
            query_end: 80,
            matching_bases: 20 * 3,
            ref_id: 1,
            score: 50.0,
            is_revcomp: false,
            anchors: vec![
                Anchor { ref_id: 1, ref_start: 60, query_start: 60 },
                Anchor { ref_id: 1, ref_start: 30, query_start: 30 },
                Anchor { ref_id: 1, ref_start:  0, query_start:  0 },
            ],
        };
        #[rustfmt::skip]
        let elsewhere = [
            Anchor { ref_id: 0, ref_start: 90, query_start: 90 },
            Anchor { ref_id: 2, ref_start: 90, query_start: 90 },
        ];

        chainer.recover_into_chain(&mut chain, &elsewhere, 20, 200);

        assert_eq!(chain.anchors.len(), 3);
        assert_eq!(chain.score, 50.0);
    }

    #[test]
    fn an_anchor_too_far_from_the_chain_is_left_out() {
        let chainer = Chainer::new(20, ChainingParameters::default());
        #[rustfmt::skip]
        let mut chain = Chain {
            id: 0,
            ref_start: 0,
            ref_end: 80,
            query_start: 0,
            query_end: 80,
            matching_bases: 20 * 3,
            ref_id: 0,
            score: 50.0,
            is_revcomp: false,
            anchors: vec![
                Anchor { ref_id: 0, ref_start: 60, query_start: 60 },
                Anchor { ref_id: 0, ref_start: 30, query_start: 30 },
                Anchor { ref_id: 0, ref_start:  0, query_start:  0 },
            ],
        };
        #[rustfmt::skip]
        let far = [Anchor { ref_id: 0, ref_start: 100000, query_start: 90 }];

        chainer.recover_into_chain(&mut chain, &far, 20, 200);

        assert_eq!(chain.anchors.len(), 3);
        assert_eq!(chain.score, 50.0);
    }
}
