use std::time::Instant;

use log::trace;

use crate::details::NamDetails;
use crate::hit::{Hit, HitsDetails, find_hits};
use crate::index::{IndexEntry, StrobemerIndex};
use crate::mcsstrategy::McsStrategy;
use crate::nam::Nam;
use crate::seeding::QueryRandstrobe;

const N_PRECOMPUTED: usize = 1024;

#[derive(Debug, Ord, PartialOrd, Eq, PartialEq, Clone, Copy)]
pub struct Anchor {
    pub ref_id: usize,
    pub ref_start: usize,
    pub query_start: usize,
}

#[derive(Debug, Clone)]
pub struct ChainingParameters {
    pub max_lookback: usize,
    pub diag_diff_penalty: f32,
    pub gap_length_penalty: f32,
    pub valid_score_threshold: f32,
    pub max_ref_gap: usize,
    pub matches_weight: f32,
}

impl Default for ChainingParameters {
    fn default() -> Self {
        ChainingParameters {
            max_lookback: 50,
            diag_diff_penalty: 0.1,
            gap_length_penalty: 0.05,
            valid_score_threshold: 0.7,
            max_ref_gap: 1000,
            matches_weight: 0.01,
        }
    }
}

#[derive(Debug, Default)]
pub struct ChainingResult {
    best_score: f32,
    dp: Vec<f32>,
    predecessors: Vec<usize>,
    anchors: Vec<Anchor>,
    parameters: ChainingParameters,
}

#[derive(Debug)]
pub struct Chainer {
    k: usize,
    parameters: ChainingParameters,
    precomputed_scores: [f32; N_PRECOMPUTED],
}

impl Chainer {
    pub fn new(k: usize, parameters: ChainingParameters) -> Self {
        let mut precomputed_scores = [0f32; N_PRECOMPUTED];
        for i in 0..N_PRECOMPUTED {
            precomputed_scores[i] = compute_score(i, i, k, &parameters);
        }
        Chainer {
            k,
            parameters,
            precomputed_scores,
        }
    }

    fn compute_score_cached(&self, dq: usize, dr: usize) -> f32 {
        // dq == dr is usually the most common case
        if dq == dr && dq < N_PRECOMPUTED {
            self.precomputed_scores[dq]
        } else {
            compute_score(dq, dr, self.k, &self.parameters)
        }
    }

    fn collinear_chaining(&self, anchors: Vec<Anchor>) -> ChainingResult {
        let n = anchors.len();
        if n == 0 {
            return ChainingResult::default();
        }

        let mut dp = vec![self.k as f32; n];
        let mut predecessors = vec![usize::MAX; n];
        let mut best_score = 0.0;
        let mut best_index = usize::MAX;

        for i in 0..n {
            let lookup_end = i.saturating_sub(self.parameters.max_lookback);
            let ai = &anchors[i];
            for j in (lookup_end..i).rev() {
                let aj = &anchors[j];

                if ai.ref_id != aj.ref_id {
                    break;
                }

                let Some(dq) = ai.query_start.checked_sub(aj.query_start) else {
                    // Not collinear
                    continue;
                };

                debug_assert!(
                    ai.ref_start >= aj.ref_start,
                    "anchors must be sorted by reference start position"
                );

                let dr = ai.ref_start - aj.ref_start;

                if dr >= self.parameters.max_ref_gap {
                    break;
                }

                let score = self.compute_score_cached(dq, dr);
                let new_score = dp[j] + score;
                if new_score >= dp[i] {
                    dp[i] = new_score;
                    predecessors[i] = j;
                    // Runtime heuristic: If the predecessor is on the same diagonal,
                    // assume that it is the best one and skip the remaining ones.
                    if dq == dr {
                        break;
                    }
                }
            }

            // The above runtime heuristic sometimes skips the anchor that
            // represents the so-far optimal chain and can then give suboptimal
            // results. To mitigate the issue, we explicitly check that anchor.
            if best_index != usize::MAX && ai.ref_id == anchors[best_index].ref_id {
                let aj = &anchors[best_index];
                if ai.query_start > aj.query_start {
                    let dq = ai.query_start - aj.query_start;
                    let dr = ai.ref_start - aj.ref_start;

                    if dr < self.parameters.max_ref_gap && dr > 0 {
                        let score = self.compute_score_cached(dq, dr);
                        let new_score = dp[best_index] + score;
                        if new_score > dp[i] {
                            dp[i] = new_score;
                            predecessors[i] = best_index;
                        }
                    }
                }
            }

            if dp[i] > best_score {
                best_score = dp[i];
                best_index = i;
            }
        }

        ChainingResult {
            best_score,
            dp,
            predecessors,
            anchors,
            parameters: self.parameters.clone(),
        }
    }

    pub fn get_chains(
        &self,
        query_randstrobes: &[Vec<QueryRandstrobe>; 2],
        index: &StrobemerIndex,
        mcs_strategy: McsStrategy,
    ) -> (NamDetails, Vec<Nam>) {
        let hits_timer = Instant::now();

        let mut hits = [vec![], vec![]];
        let mut hits_details = [HitsDetails::default(), HitsDetails::default()];
        for is_revcomp in 0..2 {
            (hits_details[is_revcomp], hits[is_revcomp]) = find_hits(
                &query_randstrobes[is_revcomp],
                index,
                mcs_strategy,
                index.filter_cutoff(),
            );
        }
        let mut time_find_hits = hits_timer.elapsed().as_secs_f64();

        // Runtime heuristic: If one orientation appears to have many more hits
        // than the other, we assume it is the correct one and do not check the
        // other.
        let mut orientations = vec![];
        if hits_details[0].is_better_than(&hits_details[1]) {
            orientations.push(0);
        } else if hits_details[1].is_better_than(&hits_details[0]) {
            orientations.push(1);
        } else {
            orientations.push(0);
            orientations.push(1);
        }

        let mut n_anchors = 0;
        let mut time_chaining = 0.0;
        let mut chains = vec![];
        // Allocate once and reuse across both orientations to avoid per-call heap alloc.
        let mut anchor_diags: Vec<(usize, isize)> = Vec::new();
        for &is_revcomp in &orientations {
            let hits_timer = Instant::now();
            let mut anchors = hits_to_anchors(&hits[is_revcomp], index);
            add_repetitive_via_sweep(&mut anchors, &hits[is_revcomp], index, &mut anchor_diags);
            time_find_hits += hits_timer.elapsed().as_secs_f64();
            n_anchors += anchors.len();
            let chaining_timer = Instant::now();
            trace!("Chaining {} anchors", anchors.len());
            anchors.sort_unstable_by_key(|a| (a.ref_id, a.ref_start, a.query_start));
            anchors.dedup();

            // trace!(
            //     "Anchors for {} strand [{}]",
            //     if is_revcomp == 0 {
            //         "forward"
            //     } else {
            //         "reverse"
            //     },
            //     anchors
            //         .iter()
            //         .map(|anchor| format!("{:?}", anchor))
            //         .collect::<Vec<_>>()
            //         .join(",")
            // );

            let chaining_result = self.collinear_chaining(anchors);

            chaining_result.extract_chains(index.k(), is_revcomp == 1, &mut chains);
            time_chaining += chaining_timer.elapsed().as_secs_f64();
        }
        let mut hits_details12 = hits_details[0].clone();
        hits_details12 += hits_details[1].clone();
        let details = NamDetails {
            hits: hits_details12,
            n_reads: 1,
            n_randstrobes: query_randstrobes[0].len() + query_randstrobes[1].len(),
            n_anchors,
            n_nams: chains.len(),
            time_randstrobes: 0.0,
            time_find_hits,
            time_chaining,
            time_rescue: 0.0,
            time_sort_nams: 0f64,
            both_orientations: orientations.len() > 1,
        };

        (details, chains)
    }
}

/// Returns the chaining score between two anchors based on their distance.
///
/// This function calculates the score contribution when chaining two anchors,
/// considering their relative positions on the query and reference sequences.
/// It penalizes large gaps and diagonal differences to encourage collinear chains.
///
/// # Parameters
///
/// - `dq`: Difference in query start positions between anchor i and anchor j
///   (i.e., ai.query_start - aj.query_start). Must be > 0.
/// - `dr`: Difference in reference start positions between anchor i and anchor j
///   (i.e., ai.ref_start - aj.ref_start). Must be > 0.
/// - `k`: Length of the k-mer used to form the anchor.
/// - `chaining_params`: Parameters controlling penalties:
///   - `diag_diff_penalty`: Multiplier for the absolute difference |dr - dq|, penalizing non-diagonal moves.
///   - `gap_length_penalty`: Multiplier for min(dq, dr), penalizing longer gaps.
fn compute_score(dq: usize, dr: usize, k: usize, parameters: &ChainingParameters) -> f32 {
    let dd = dr.abs_diff(dq);
    let dg = dq.min(dr);
    let mut score = k.min(dg) as f32;
    let dg = dg.saturating_sub(k);

    let lin_penalty =
        parameters.diag_diff_penalty * dd as f32 + parameters.gap_length_penalty * dg as f32;
    let log_penalty = if dd >= 1 {
        (dd as f32 + 1f32).log2()
    } else {
        0.0
    };
    score -= lin_penalty + 0.5 * log_penalty;

    score
}

fn add_to_anchors_full(
    anchors: &mut Vec<Anchor>,
    query_start: usize,
    query_end: usize,
    index: &StrobemerIndex,
    entry: IndexEntry,
) {
    let mut min_length_diff = usize::MAX;
    let forward_hash = entry.hash();
    for i in entry.position..index.len() {
        let entry = index.entry(i);
        if entry.hash() != forward_hash {
            break;
        }
        let ref_start = entry.position();
        let ref_end = ref_start + entry.strobe2_offset() + index.k();
        let length_diff = (query_end - query_start).abs_diff(ref_end - ref_start);
        if length_diff <= min_length_diff {
            let ref_id = entry.reference_index();
            anchors.push(Anchor {
                ref_id,
                ref_start,
                query_start,
            });
            anchors.push(Anchor {
                ref_id,
                ref_start: ref_end - index.k(),
                query_start: query_end - index.k(),
            });
            min_length_diff = length_diff;
        }
    }
}

fn add_to_anchors_partial(
    anchors: &mut Vec<Anchor>,
    query_start: usize,
    index: &StrobemerIndex,
    entry: IndexEntry,
) {
    let forward_hash = entry.get_hash_partial_forward();
    for i in entry.position..index.len() {
        let entry = index.entry(i);
        // Filter out partial lookups with different orientation
        if entry.get_hash_partial_forward() != forward_hash {
            break;
        }
        let ref_id = entry.reference_index();
        let ref_start = entry.strobe_extent_partial().0;

        anchors.push(Anchor {
            ref_id,
            ref_start,
            query_start,
        });
    }
}

/// For each filtered (repetitive) hit, add its reference positions that land on
/// the same diagonal as an existing anchor (within `k` bp tolerance).
///
/// This lets repetitive seeds contribute to chaining without the N×R candidate
/// explosion: only positions co-linear with at least one non-repetitive anchor
/// are admitted.
///
/// # Sort-order note
/// `RefRandstrobe` is now sorted by `(hash, ref_index, position)`, so within a
/// hash group entries are grouped by chromosome and then sorted by position.
/// `anchor_diags` is sorted by `(ref_id, diagonal)`.  Both sequences increase
/// in the same (ref_id, value) order, enabling a single two-pointer sweep with
/// no per-ref_id bookkeeping.
fn add_repetitive_via_sweep(
    anchors: &mut Vec<Anchor>,
    hits: &[Hit],
    index: &StrobemerIndex,
    anchor_diags: &mut Vec<(usize, isize)>,
) {
    if anchors.is_empty() {
        return;
    }

    // Select a minimum-cost non-overlapping tiling of filtered hits via
    // Weighted Interval Scheduling (WIS) where weight = CAP − count + 1.
    // Two hits may overlap by at most OVERLAP_THRESHOLD bp on the query.
    // The cap is applied here — hits above it are excluded before WIS,
    // avoiding a redundant get_count_full_forward call in the sweep loop.
    //
    // DP (O(n log n)):
    //   Sort by query_end.  For each hit i, p[i] = number of hits whose
    //   query_end ≤ query_start[i] + OVERLAP_THRESHOLD (the compatible
    //   prefix).  Then:
    //     dp[i] = max(dp[i−1],  weight[i] + dp[p[i]])
    //   where dp[0] = 0 and weight[i] = CAP − count[i] + 1 > 0.
    const CAP: usize = 10_000;
    let mut candidates: Vec<(&Hit, usize)> = hits
        .iter()
        .filter(|h| h.is_filtered && !h.is_partial)
        .filter_map(|h| {
            let count = index.entry(h.position).get_count_full_forward();
            if count <= CAP { Some((h, count)) } else { None }
        })
        .collect();

    // Fast path: nothing to do if no eligible filtered hits remain.
    if candidates.is_empty() {
        return;
    }
    // Sort by query_end for the WIS DP; ties broken by count.
    candidates.sort_unstable_by_key(|&(h, count)| (h.query_end, count));

    const OVERLAP_THRESHOLD: usize = 10;
    let n = candidates.len();

    // p[i] = number of hits compatible with (i.e. not overlapping) hit i.
    let p: Vec<usize> = (0..n)
        .map(|i| {
            let cutoff = candidates[i].0.query_start + OVERLAP_THRESHOLD;
            candidates[..i].partition_point(|&(h, _)| h.query_end <= cutoff)
        })
        .collect();

    // weight[i] = CAP - count[i] + 1: positive for all hits (count ≤ CAP),
    // higher for less-repetitive seeds.  Ensures the DP always selects
    // at least one hit (empty set has value 0 < any single hit's weight),
    // while two non-overlapping low-count hits beat one high-count hit
    // when their combined weight exceeds it.
    let weight = |i: usize| (CAP + 1 - candidates[i].1.min(CAP)) as isize;

    // dp[i] = maximum total weight over a valid tiling of hits 0..i.
    let mut dp = vec![0isize; n + 1];
    for i in 1..=n {
        dp[i] = dp[i - 1].max(weight(i - 1) + dp[p[i - 1]]);
    }

    // Traceback to recover selected hits.
    let mut selected: Vec<&Hit> = Vec::new();
    let mut i = n;
    while i > 0 {
        if weight(i - 1) + dp[p[i - 1]] >= dp[i - 1] {
            selected.push(candidates[i - 1].0);
            i = p[i - 1];
        } else {
            i -= 1;
        }
    }

    // Build (ref_id, diagonal) sorted by (ref_id, diagonal).
    // diagonal = ref_start - query_start; co-linear hits share the same value.
    // Reuse the caller-provided buffer to avoid per-call heap allocation.
    anchor_diags.clear();
    anchor_diags.extend(
        anchors
            .iter()
            .map(|a| (a.ref_id, a.ref_start as isize - a.query_start as isize)),
    );
    anchor_diags.sort_unstable();
    anchor_diags.dedup();

    let tolerance = index.k() as isize;

    for hit in &selected {
        let q = hit.query_start as isize;
        let entry = index.entry(hit.position);
        let forward_hash = entry.hash();

        // Skip seeds so repetitive that even co-linear filtering would flood the chainer.
        if entry.get_count_full_forward() > 10_000 {
            continue;
        }

        // Single anchor pointer — never goes backwards because:
        // - index positions are now sorted by (ref_id, ref_start) within a hash group
        // - anchor_diags is sorted by (ref_id, diagonal)
        // - both sequences increase in the same (ref_id, value) direction
        let mut ptr_anc = 0usize;

        for pos in entry.position..index.len() {
            let entry = index.entry(pos);
            if entry.hash() != forward_hash {
                break;
            }
            let ref_id = entry.reference_index();
            let ref_start = entry.position();
            let diag = ref_start as isize - q;
            let lo = (ref_id, diag - tolerance);
            let hi = (ref_id, diag + tolerance);

            // Advance past anchor entries that are before the window.
            // Tuple comparison handles both ref_id changes and diagonal range.
            while ptr_anc < anchor_diags.len() && anchor_diags[ptr_anc] < lo {
                ptr_anc += 1;
            }

            if ptr_anc < anchor_diags.len() && anchor_diags[ptr_anc] <= hi {
                let ref_end = ref_start + entry.strobe2_offset() + index.k();
                anchors.push(Anchor {
                    ref_id,
                    ref_start,
                    query_start: hit.query_start,
                });
                anchors.push(Anchor {
                    ref_id,
                    ref_start: ref_end - index.k(),
                    query_start: hit.query_end - index.k(),
                });
            }
        }
    }
}

fn hits_to_anchors(hits: &Vec<Hit>, index: &StrobemerIndex) -> Vec<Anchor> {
    let mut anchors = vec![];
    for hit in hits {
        if hit.is_filtered {
            continue;
        }
        if hit.is_partial {
            if let Some(forward_entry) = index.get_partial_forward_from(hit.hash, hit.position) {
                add_to_anchors_partial(&mut anchors, hit.query_start, index, forward_entry);
            }
        } else {
            add_to_anchors_full(
                &mut anchors,
                hit.query_start,
                hit.query_end,
                index,
                index.entry(hit.position),
            );
        }
    }

    anchors
}

impl ChainingResult {
    fn extract_chains(&self, k: usize, is_revcomp: bool, chains: &mut Vec<Nam>) {
        let n = self.anchors.len();
        let valid_score = self.best_score * self.parameters.valid_score_threshold;

        let mut candidates = vec![];
        for i in 0..n {
            if self.dp[i] >= valid_score {
                candidates.push((i, self.dp[i]));
            }
        }

        candidates.sort_by(|a, b| b.1.total_cmp(&a.1));

        let mut used = vec![false; n];
        for (i, score) in candidates {
            if used[i] {
                continue;
            }

            let mut j = i;
            let mut overlaps = false;
            let mut chain_anchors = vec![self.anchors[i]];

            let mut matching_bases = k;
            let mut ref_coverage = self.anchors[i].ref_start;

            while self.predecessors[j] != usize::MAX {
                j = self.predecessors[j];
                if used[j] {
                    overlaps = true;
                    break;
                }
                chain_anchors.push(self.anchors[j]);
                used[j] = true;

                matching_bases += ref_coverage
                    .saturating_sub(self.anchors[j].ref_start)
                    .min(k);
                ref_coverage = self.anchors[j].ref_start;
            }

            if overlaps {
                continue;
            }

            let first = &self.anchors[j];
            let last = &self.anchors[i];

            chains.push(Nam {
                nam_id: chains.len(),
                query_start: first.query_start,
                query_end: last.query_start + k,
                ref_start: first.ref_start,
                ref_end: last.ref_start + k,
                matching_bases,
                ref_id: last.ref_id,
                score: score + chain_anchors.len() as f32 * self.parameters.matches_weight,
                is_revcomp,
                anchors: chain_anchors,
            });
        }
    }
}

#[cfg(test)]
mod test {
    use super::{Anchor, Chainer, ChainingParameters};

    #[test]
    fn test_chainer_early_break() {
        let chainer = Chainer::new(20, ChainingParameters::default());
        #[rustfmt::skip]
        let anchors = vec![
            Anchor { ref_id: 0, ref_start:  0, query_start:  0, },
            Anchor { ref_id: 0, ref_start: 30, query_start: 20, },
            Anchor { ref_id: 0, ref_start: 60, query_start:  0, },
            Anchor { ref_id: 0, ref_start: 95, query_start: 35, },
        ];

        let chaining_result = chainer.collinear_chaining(anchors);

        // The best chain has score 42.342842 and uses anchors 0, 1, 3.
        // When using the heuristic that breaks early if the predecessor is on the
        // same diagonal, a suboptimal chain is found that has score 38.2 and
        // consists of anchors 2 and 3.
        assert!(chaining_result.best_score > 42.0);
    }

    #[test]
    fn test_linear_score_adjacent_anchors() {
        let chainer = Chainer::new(20, ChainingParameters::default());
        #[rustfmt::skip]
        let anchors = [
            Anchor { ref_id: 0, ref_start:   0, query_start:  0, },
            Anchor { ref_id: 0, ref_start:  20, query_start: 20, },
            Anchor { ref_id: 0, ref_start:  40, query_start: 40, },
        ];
        let chaining_result = chainer.collinear_chaining(anchors[0..1].to_vec());
        let score1 = chaining_result.best_score;
        assert_eq!(
            chainer
                .collinear_chaining(anchors[0..2].to_vec())
                .best_score,
            score1 * 2.0
        );
        assert_eq!(
            chainer
                .collinear_chaining(anchors[0..3].to_vec())
                .best_score,
            score1 * 3.0
        );
    }
}
