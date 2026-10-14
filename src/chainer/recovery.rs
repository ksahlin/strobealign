use super::{Anchor, Chainer};
use crate::chain::Chain;
use crate::hit::Hit;
use crate::index::{REF_RANDSTROBE_HASH_MASK, RefRandstrobe, StrobemerIndex};
use crate::partition::custom_partition_point;
use crate::refseq::ContigStarts;

const MAX_BAND: usize = 1000;
const MAX_LOOKBACK: usize = 32;
/// Reference positions are binned into bins of 2^BIN_BITS, which wrap around
/// after N_BINS
const BIN_BITS: u32 = 12;
const N_BINS: usize = 1 << 15;
/// Entries skipped at once
const SKIP: usize = 32;
/// Of the filtered hits, at most one per read_len / HIT_SPACING_DIVISOR bases
/// of the read is looked up: in satellite arrays, they are dense
const HIT_SPACING_DIVISOR: usize = 400;

/// Buffers reused for the chains of a read: the candidates of a chain, its
/// anchors with the recovered ones and the DP of a gap (see take_best_run)
type Scratch = (Vec<Anchor>, Vec<Anchor>, Vec<(f32, usize, f32)>);

impl Chainer {
    /// Adds anchors from the filtered (too frequent) hits to the chains (all
    /// of the same orientation as the hits): the hits are looked up again,
    /// restricted to the diagonals near the chains, and the best run of the
    /// resulting anchors is inserted into each gap.
    pub fn recover_anchors(
        &self,
        chains: &mut [Chain],
        hits: &[Hit],
        index: &StrobemerIndex,
        read_len: usize,
    ) {
        if chains.is_empty() {
            return;
        }
        let bands = chains.iter().map(|chain| band(chain, read_len)).collect();
        let spacing = read_len / HIT_SPACING_DIVISOR;
        let mut anchors = filtered_hits_to_anchors(hits, index, bands, spacing);
        anchors.sort_unstable();
        anchors.dedup();

        let (contig_starts, mut scratch) = (&index.contig_starts, Scratch::default());
        for chain in chains {
            self.recover_into_chain(chain, &anchors, contig_starts, read_len, &mut scratch);
        }
    }

    /// Inserts the best run of the (sorted) anchors near the chain into each
    /// of its gaps and before and after it
    fn recover_into_chain(
        &self,
        chain: &mut Chain,
        anchors: &[Anchor],
        contig_starts: &ContigStarts,
        read_len: usize,
        (candidates, recovered, dp): &mut Scratch,
    ) {
        // Anchors on a neighbouring contig are not candidates
        let first = chain
            .ref_start
            .saturating_sub(read_len)
            .max(chain.ref_contig_start);
        let last = (chain.ref_end + read_len).min(contig_starts.contig_end(chain.ref_start) - 1);
        let start = anchors.partition_point(|a| a.ref_start() < first);
        let (low, high) = band(chain, read_len);
        let diagonal = |a: &&Anchor| a.ref_start() as i64 - a.query_start() as i64;
        candidates.clear();
        candidates.extend(
            anchors[start..]
                .iter()
                .take_while(|a| a.ref_start() <= last)
                .filter(|a| (low..=high).contains(&diagonal(a))),
        );
        if candidates.is_empty() {
            return;
        }

        recovered.clear();
        recovered.reserve(chain.anchors.len() + candidates.len());
        let mut gap_start = 0;
        let mut gained = 0f32;
        for next in chain.anchors.iter().rev().copied().map(Some).chain([None]) {
            let gap_end = next.map_or(candidates.len(), |next| {
                gap_start
                    + candidates[gap_start..].partition_point(|a| a.ref_start() < next.ref_start())
            });
            let gap = &candidates[gap_start..gap_end];
            gained += self.take_best_run(gap, next, recovered, dp, read_len);
            recovered.extend(next);
            gap_start = gap_end;
        }
        if recovered.len() == chain.anchors.len() {
            return;
        }

        let k = self.k;
        recovered.reverse();
        chain.score += gained;
        let covered = |p: &[Anchor]| p[0].ref_start().saturating_sub(p[1].ref_start()).min(k);
        chain.matching_bases = k + recovered.windows(2).map(covered).sum::<usize>();
        let (first, last) = (recovered.last().unwrap(), recovered[0]);
        (chain.ref_start, chain.query_start) = (first.ref_start(), first.query_start());
        (chain.ref_end, chain.query_end) = (last.ref_start() + k, last.query_start() + k);
        // The chain's previous anchors become the buffer for the next chain
        std::mem::swap(&mut chain.anchors, recovered);
    }

    fn link_score(&self, from: Anchor, to: Anchor, read_len: usize) -> Option<f32> {
        let dq = to.query_start().checked_sub(from.query_start())?;
        let dr = to.ref_start().checked_sub(from.ref_start())?;
        if dq == 0 || dr == 0 || dr >= self.parameters.max_ref_gap.unwrap_or(read_len) {
            return None;
        }
        if dq.max(dr) as f32 / dq.min(dr) as f32 > self.parameters.max_diagonal_ratio {
            return None;
        }

        Some(self.compute_score_cached(dq, dr))
    }

    /// An upper bound of the score of a link: its score without the
    /// logarithmic penalty, which is non-negative, computed with the same
    /// floating-point operations as compute_score
    fn link_bound(&self, dq: usize, dr: usize) -> f32 {
        let (dd, dg) = (dr.abs_diff(dq), dq.min(dr));
        let penalty = self.parameters.diag_diff_penalty * dd as f32
            + self.parameters.gap_length_penalty * dg.saturating_sub(self.k) as f32;
        self.k.min(dg) as f32 - penalty
    }

    /// Appends the run of anchors from gap that best links the last
    /// recovered anchor to next, if it scores better than linking them
    /// directly, and returns how much it adds to the score. dp holds for each
    /// anchor the best score of a run ending there, its predecessor in that
    /// run and the best score of a run ending there or before.
    fn take_best_run(
        &self,
        gap: &[Anchor],
        next: Option<Anchor>,
        recovered: &mut Vec<Anchor>,
        dp: &mut Vec<(f32, usize, f32)>,
        read_len: usize,
    ) -> f32 {
        if gap.is_empty() {
            return 0f32;
        }
        let previous = recovered.last().copied();
        let replaced_link = match (previous, next) {
            (Some(previous), Some(next)) => self.compute_score_cached(
                next.query_start() - previous.query_start(),
                next.ref_start() - previous.ref_start(),
            ),
            _ => 0f32,
        };

        let parameters = &self.parameters;
        let weight = parameters.matches_weight;
        // With non-negative penalties, a link scores at most k
        let bounded = parameters.diag_diff_penalty >= 0.0 && parameters.gap_length_penalty >= 0.0;
        let max_link = if bounded {
            self.k as f32
        } else {
            f32::INFINITY
        };
        dp.clear();
        let (mut best_total, mut best_end) = (replaced_link, usize::MAX);
        let mut best_so_far = f32::NEG_INFINITY;
        // The first of the anchors at the reference position of the current
        // one, which cannot link to it
        let mut same_ref_start = 0;
        for (i, &anchor) in gap.iter().enumerate() {
            let (q, r) = (anchor.query_start(), anchor.ref_start());
            if gap[same_ref_start].ref_start() != r {
                same_ref_start = i;
            }
            // A run has to advance from previous and end before next, and an
            // anchor out of these bounds could only lead to other such anchors
            if previous.is_some_and(|p| q <= p.query_start() || r <= p.ref_start())
                || next.is_some_and(|n| q >= n.query_start())
            {
                dp.push((f32::NEG_INFINITY, usize::MAX, best_so_far));
                continue;
            }
            let opening_link =
                previous.map_or(Some(0f32), |p| self.link_score(p, anchor, read_len));
            let mut best = opening_link.unwrap_or(f32::NEG_INFINITY) + weight;
            let mut predecessor = usize::MAX;
            // Going backwards, stop when no earlier predecessor can reach best
            // even with the largest link. Of equal scores, the earliest wins.
            for j in (i.saturating_sub(MAX_LOOKBACK)..same_ref_start).rev() {
                let (score, _, best_up_to_j) = dp[j];
                if best_up_to_j + max_link + weight < best {
                    break;
                }
                let from = gap[j];
                let Some(dq) = q.checked_sub(from.query_start()) else {
                    continue;
                };
                if score + self.link_bound(dq, r - from.ref_start()) + weight < best {
                    continue;
                }
                if let Some(link) = self.link_score(from, anchor, read_len) {
                    let score = score + link + weight;
                    if score > best || (score == best && predecessor != usize::MAX) {
                        (best, predecessor) = (score, j);
                    }
                }
            }
            best_so_far = best_so_far.max(best);
            dp.push((best, predecessor, best_so_far));
            let closing_link = next.map_or(Some(0f32), |n| self.link_score(anchor, n, read_len));
            if let Some(closing_link) = closing_link
                && best + closing_link > best_total
            {
                best_total = best + closing_link;
                best_end = i;
            }
        }

        let start = recovered.len();
        let mut i = best_end;
        while i != usize::MAX {
            recovered.push(gap[i]);
            i = dp[i].1;
        }
        recovered[start..].reverse();

        best_total - replaced_link
    }
}

/// Returns the diagonals the chain spans, widened by a band
fn band(chain: &Chain, read_len: usize) -> (i64, i64) {
    let band = (read_len / 2).min(MAX_BAND) as i64;
    let start_diagonal = chain.ref_start as i64 - chain.query_start as i64;
    let end_diagonal = chain.ref_end as i64 - chain.query_end as i64;

    (
        start_diagonal.min(end_diagonal) - band,
        start_diagonal.max(end_diagonal) + band,
    )
}

/// Merges overlapping (closed) intervals
fn merge(intervals: &mut Vec<(i64, i64)>) {
    intervals.sort_unstable();
    intervals.dedup_by(|next, previous| {
        next.0 <= previous.1 && {
            previous.1 = previous.1.max(next.1);
            true
        }
    });
}

/// Looks up the filtered hits again, skipping those closer than spacing to
/// the previous one looked up, and returns their anchors, among them
/// all those on the diagonals of the bands
fn filtered_hits_to_anchors(
    hits: &[Hit],
    index: &StrobemerIndex,
    mut bands: Vec<(i64, i64)>,
    spacing: usize,
) -> Vec<Anchor> {
    let k = index.k();
    let max_dist = index.parameters.randstrobe.max_dist as i64;
    merge(&mut bands);
    // The first index entry of each filtered hit, and the mask of the hash
    // its entries share
    let mut runs = vec![];
    let mut next = 0;
    for hit in hits.iter().filter(|hit| hit.is_filtered) {
        if hit.query_start < next {
            continue;
        }
        next = hit.query_start + spacing;
        if !hit.is_partial {
            runs.push((hit.position, REF_RANDSTROBE_HASH_MASK, hit));
        } else if let Some(entry) = index.get_partial_forward_from(hit.hash, hit.position) {
            let mask = index.parameters.randstrobe.forward_main_hash_mask;
            runs.push((entry.position, mask, hit));
        }
    }
    runs.sort_unstable_by_key(|&(start, mask, _)| (start, mask));

    let mut anchors = vec![];
    let mut windows = vec![];
    // Hits that share their entries (a repeat in the read) scan them once
    for group in runs.chunk_by(|a, b| (a.0, a.1) == (b.0, b.1)) {
        // An entry at r gives an anchor on diagonal r - query_start and, for
        // a full hit, one on diagonal r + offset - (query_end - k), where both
        // offset and query_end - k - query_start are in 0..=max_dist. So each
        // anchor in a band comes from an entry in a window of its hit.
        windows.clear();
        for &(_, _, hit) in group {
            let (low, high) = match hit.is_partial {
                true => (0, 0),
                false => {
                    let query_offset = (hit.query_end - k - hit.query_start) as i64;
                    (query_offset - max_dist, query_offset)
                }
            };
            let q = hit.query_start as i64;
            windows.extend(
                bands
                    .iter()
                    .map(|&(start, end)| (start + q + low, end + q + high)),
            );
        }
        merge(&mut windows);
        // Most entries are rejected by one lookup in a bitmap of the bins the
        // windows touch, and the others by a search of the windows
        let mut bitmap = [0u64; N_BINS / 64];
        for &(start, end) in &windows {
            let bins = start.max(0) as usize >> BIN_BITS..=end.max(0) as usize >> BIN_BITS;
            for bin in bins.take(N_BINS) {
                bitmap[bin % N_BINS / 64] |= 1 << (bin % 64);
            }
        }
        let (start, mask, _) = group[0];
        let entries = index.entries_from(start);
        let key = entries[0].hash() & mask;
        let in_group = |e: &RefRandstrobe| e.hash() & mask == key;
        let marked = |e: &RefRandstrobe| {
            let bin = e.ref_start() >> BIN_BITS;
            bitmap[bin % N_BINS / 64] & (1 << (bin % 64)) != 0
        };
        let mut i = 0;
        while i < entries.len() {
            // The entries are looked up SKIP at a time. Those of the group are
            // contiguous, so if the last of a chunk is in the group, all are.
            let chunk = &entries[i..entries.len().min(i + SKIP)];
            let n = match chunk.last().is_some_and(in_group) {
                true => chunk.len(),
                false => chunk.partition_point(in_group),
            };
            let found = chunk[..n].iter().position(marked);
            if found.is_none() && n < SKIP {
                break;
            }
            i += found.unwrap_or(SKIP - 1);
            let entry = &entries[i];
            let r = entry.ref_start();
            i += 1;
            // Entries of the same hash and offset are sorted by position, so
            // after SKIP rejected ones, those before the next window can be
            // skipped if the run goes on
            let offset = entry.strobe2_offset() as u8;
            let run_bound = |ref_start| RefRandstrobe::new(entry.hash(), ref_start, offset);
            let run_end = run_bound(usize::MAX);
            if found.is_none() && !entries.get(i + SKIP).is_some_and(|e| *e < run_end) {
                continue;
            }
            let w = custom_partition_point(&windows, |&(_, end)| end < r as i64);
            let next_window = windows.get(w).map_or(i64::MAX, |&(start, _)| start);
            if next_window <= r as i64 {
                for &(_, _, hit) in group {
                    anchors.push(Anchor::new(r, hit.query_start));
                    if !hit.is_partial {
                        anchors.push(Anchor::new(r + entry.strobe2_offset(), hit.query_end - k));
                    }
                }
            } else {
                let bound = run_bound(next_window as usize);
                while entries.get(i + SKIP).is_some_and(|e| *e < bound) {
                    i += SKIP;
                }
            }
        }
    }

    anchors
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::anchors;
    use crate::chainer::ChainingParameters;
    use crate::index::RefRandstrobe;
    use crate::seeding::SeedingParameters;

    /// Returns a chain of the given anchors (last first)
    fn chain(anchors: Vec<Anchor>, k: usize, score: f32) -> Chain {
        let (first, last) = (anchors.last().unwrap(), anchors[0]);
        Chain {
            id: 3,
            ref_start: first.ref_start(),
            ref_end: last.ref_start() + k,
            query_start: first.query_start(),
            query_end: last.query_start() + k,
            matching_bases: k * anchors.len(),
            ref_contig_start: 0,
            score,
            is_revcomp: false,
            anchors,
        }
    }

    /// Recovers candidates into a chain of the given anchors (last first)
    fn recover(
        chain_anchors: Vec<Anchor>,
        candidates: Vec<Anchor>,
        contigs: Vec<usize>,
        read_len: usize,
        score: f32,
    ) -> Chain {
        let chainer = Chainer::new(20, ChainingParameters::default());
        let mut chain = chain(chain_anchors, 20, score);
        let contig_starts = ContigStarts::new(contigs, 1_000_000);
        let scratch = &mut Default::default();
        chainer.recover_into_chain(&mut chain, &candidates, &contig_starts, read_len, scratch);

        chain
    }

    #[test]
    fn a_filtered_hit_is_recovered_near_a_chain() {
        let hash = 0x0123_4567_89AB_CD00u64 & REF_RANDSTROBE_HASH_MASK;
        // The entries before 1090 are skipped, and the last one is in the
        // same bin of the bitmap as 1090
        let randstrobes: Vec<_> = (0..100)
            .chain([1090])
            .chain((1..40).map(|i| i * 50_000))
            .chain([1090 + (N_BINS << BIN_BITS)])
            .map(|ref_start| RefRandstrobe::new(hash, ref_start, 10))
            .collect();
        let n = randstrobes.len();
        let contig_starts = ContigStarts::new(vec![0], (2 * N_BINS) << BIN_BITS);
        let parameters = SeedingParameters::new(300);
        let index = StrobemerIndex::new(
            parameters,
            1,
            100,
            randstrobes,
            vec![0, n, n],
            contig_starts,
        );
        let k = index.k();
        let chainer = Chainer::new(k, ChainingParameters::default());
        #[rustfmt::skip]
        let hit = Hit {
            query_start: 90, query_end: 100 + k, position: 0, hash, hash_revcomp: 0,
            is_partial: false, is_filtered: true,
        };
        let mut chains = [chain(anchors![(1180, 180), (1000, 0)], k, 0.0)];

        chainer.recover_anchors(&mut chains, &[hit], &index, 300);

        #[rustfmt::skip]
        assert_eq!(chains[0].anchors, anchors![(1180, 180), (1100, 100), (1090, 90), (1000, 0)]);
    }

    #[test]
    fn anchors_that_do_not_fit_are_left_out() {
        let three = anchors![(60, 60), (30, 30), (0, 0)];
        let five = anchors![(120, 120), (90, 90), (60, 60), (30, 30), (0, 0)];
        for (chain_anchors, candidates, contigs) in [
            // on another contig
            (
                three.clone(),
                anchors![(120, 120), (120, 120)],
                vec![0, 100],
            ),
            // scoring worse than the link it would replace
            (five.clone(), anchors![(75, 62)], vec![0]),
            // going backwards in the query
            (three.clone(), anchors![(45, 10)], vec![0]),
            // already in the chain
            (three.clone(), anchors![(0, 0), (30, 30), (60, 60)], vec![0]),
            // too far from the chain
            (three.clone(), anchors![(100000, 90)], vec![0]),
        ] {
            let chain = recover(chain_anchors.clone(), candidates, contigs, 200, 50.0);

            assert_eq!(chain.anchors, chain_anchors);
            assert_eq!(chain.matching_bases, 20 * chain_anchors.len());
            assert_eq!(chain.score, 50.0);
        }
    }

    #[test]
    fn each_end_of_a_chain_takes_more_than_one_anchor() {
        let outside = anchors![(40, 40), (70, 70), (160, 160), (190, 190)];

        let chain = recover(
            anchors![(130, 130), (100, 100)],
            outside,
            vec![0],
            200,
            39.52,
        );

        #[rustfmt::skip]
        assert_eq!(chain.anchors, anchors![
            (190, 190), (160, 160), (130, 130), (100, 100), (70, 70), (40, 40),
        ]);
        assert_eq!(chain.matching_bases, 20 * 6);
        assert_eq!((chain.ref_start, chain.ref_end), (40, 210));
        assert_eq!((chain.query_start, chain.query_end), (40, 210));
        assert!((chain.score - (20.0 + 5.0 * 19.5 + 6.0 * 0.01)).abs() < 0.0001);
    }

    #[test]
    fn the_best_of_the_candidates_for_a_gap_is_the_one_taken() {
        let competing = anchors![(100, 150), (110, 110)];

        let chain = recover(anchors![(200, 200), (0, 0)], competing, vec![0], 400, 31.02);

        assert_eq!(chain.anchors, anchors![(200, 200), (110, 110), (0, 0)]);
        assert!((chain.score - (31.02 + 15.5 + 16.5 - 11.0 + 0.01)).abs() < 0.0001);
    }

    #[test]
    fn a_gap_takes_every_anchor_that_raises_the_score() {
        let missing = anchors![(50, 50), (100, 100), (150, 150)];

        let chain = recover(anchors![(200, 200), (0, 0)], missing, vec![0], 400, 31.02);

        #[rustfmt::skip]
        assert_eq!(chain.anchors, anchors![
            (200, 200), (150, 150), (100, 100), (50, 50), (0, 0),
        ]);
        assert_eq!(chain.matching_bases, 20 * 5);
        assert!((chain.score - (20.0 + 4.0 * 18.5 + 5.0 * 0.01)).abs() < 0.0001);
    }
}
