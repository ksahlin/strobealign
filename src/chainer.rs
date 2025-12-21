use std::time::Instant;
use log::{log_enabled, Level, trace};

use crate::details::NamDetails;
use crate::hit::{find_hits, Hit, HitsDetails};
use crate::index::StrobemerIndex;
use crate::mapper::QueryRandstrobe;
use crate::mcsstrategy::McsStrategy;
use crate::nam::Nam;

const N_PRECOMPUTED: usize = 1024;

#[derive(Debug, Ord, PartialOrd, Eq, PartialEq)]
struct Anchor {
    ref_id: usize,
    ref_start: usize,
    query_start: usize,
}

#[derive(Debug)]
pub struct ChainingParameters {
    pub max_lookback: usize,
    pub diag_diff_penalty: f32,
    pub gap_length_penalty: f32,
    pub valid_score_threshold: f32,
    pub max_ref_gap: usize,
    pub matches_weight: f32,
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
        Chainer { k, parameters, precomputed_scores }
    }

    fn compute_score_cached(&self, dq: usize, dr: usize) -> f32 {
        // dq == dr is usually the most common case
        if dq == dr && (dq as usize) < N_PRECOMPUTED {
            self.precomputed_scores[dq]
        } else {
            compute_score(dq, dr, self.k, &self.parameters)
        }
    }

    fn collinear_chaining(
        &self,
        anchors: &[Anchor],
        // dp: &[f32],
        // predecessors: &[usize],
    ) -> (f32, Vec<f32>, Vec<usize>) {
        let n = anchors.len();
        if n == 0 {
            return (0.0, vec![], vec![]);
        }

        let mut dp = vec![self.k as f32; n];
        let mut predecessors = vec![usize::MAX; n];
        let mut best_score = 0.0;

        for i in 0..n {
            let lookup_end = i.saturating_sub(self.parameters.max_lookback);
            for j in (lookup_end..i).rev() {
                let ai = &anchors[i];
                let aj = &anchors[j];

                if ai.ref_id != aj.ref_id {
                    break;
                }

                if ai.query_start < aj.query_start {
                    // Not collinear
                    continue;
                }
                let dq = ai.query_start - aj.query_start;
                let dr = ai.ref_start - aj.ref_start;

                if dr >= self.parameters.max_ref_gap {
                    break;
                }
                if /*dq <= 0 ||*/ dr <= 0 {
                    // Not collinear
                    continue;
                }

                let score = self.compute_score_cached(dq, dr);

                let new_score = dp[j] + score;
                if new_score > dp[i] {
                    dp[i] = new_score;
                    predecessors[i] = j;
                    // Runtime heuristic: If the predecessor is on the same diagonal,
                    // assume that it is the best one and skip the remaining ones.
                    if dq == dr {
                        break;
                    }
                }
            }
            if dp[i] > best_score {
                best_score = dp[i];
            }
        }

        (best_score, dp, predecessors)
    }

    pub fn get_chains(
        &self,
        query_randstrobes: &[Vec<QueryRandstrobe>; 2],
        index: &StrobemerIndex,
        rescue_distance: usize,
        mcs_strategy: McsStrategy,
    ) -> (NamDetails, Vec<Nam>) {
        let hits_timer = Instant::now();

        let mut total_hits = 0;
        let mut n_rescue_hits = 0;
        let mut n_rescue_partial_hits = 0;
        let mut hits = [vec![], vec![]];
        let mut hits_details = [HitsDetails::default(), HitsDetails::default()];
        for is_revcomp in 0..2 {
            let sorting_needed1;

            (hits_details[is_revcomp], sorting_needed1, hits[is_revcomp]) =
                find_hits(&query_randstrobes[is_revcomp], index, mcs_strategy, index.filter_cutoff, rescue_distance);

            if log_enabled!(Level::Trace) {
                trace!("Found {} hits", hits_details[is_revcomp].total_hits());
                for hit in &hits[is_revcomp] {
                    trace!("Hit: {:?}", hit);
                }
            }
        }
        let mut time_find_hits = hits_timer.elapsed().as_secs_f64();

        let mut n_anchors = 0;
        let mut n_rescue_hits = 0;
        let mut n_rescue_nams = 0;
        let mut rescue_done = false;
        let mut time_rescue = 0.0;
        let mut time_chaining = 0.0;

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

        let mut chains = vec![];
        for is_revcomp in orientations {
            let hits_timer = Instant::now();
            let mut anchors = hits_to_anchors(&hits[is_revcomp], index);
            time_find_hits += hits_timer.elapsed().as_secs_f64();
            /*for anchor in anchors.iter().take(100) {
                trace!("{:?}", anchor);
            }*/
            n_anchors += anchors.len();
            let chaining_timer = Instant::now();
            trace!("Chaining {} anchors", anchors.len());
            // TODO this used to be pdqsort
            anchors.sort_by_key(|a| (a.ref_id, a.ref_start, a.query_start));
            anchors.dedup();
            let (best_score, dp, predecessors) = self.collinear_chaining(&anchors);

            extract_chains_from_dp(
                &anchors, &dp, &predecessors, best_score,
                index.k(), is_revcomp == 1, &mut chains, &self.parameters
            );
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
            n_rescue_nams,
            nam_rescue: rescue_done as usize,
            n_rescue_hits,
            n_rescue_partial_hits,
            time_randstrobes: 0f64,
            time_find_hits,
            time_chaining,
            time_rescue,
            time_sort_nams: 0f64,
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

    let lin_penalty = parameters.diag_diff_penalty * dd as f32 + parameters.gap_length_penalty * dg as f32;
    let log_penalty = if dd >= 1 { (dd as f32 + 1f32).log2() } else { 0.0 };
    score -= lin_penalty + 0.5 * log_penalty;

    score
}

fn add_to_anchors_full(
    anchors: &mut Vec<Anchor>,
    query_start: usize,
    query_end: usize,
    index: &StrobemerIndex,
    mut position: usize,
) {
    let mut min_length_diff = usize::MAX;
    let hash = index.randstrobes[position].hash();
    for randstrobe in &index.randstrobes[position..] {
        if randstrobe.hash() != hash {
            break;
        }
        let ref_start = randstrobe.position();
        let ref_end = ref_start + randstrobe.strobe2_offset() + index.k();
        let length_diff = (query_end - query_start).abs_diff(ref_end - ref_start);
        if length_diff <= min_length_diff {
            let ref_id = randstrobe.reference_index();
            anchors.push(Anchor {ref_id, ref_start, query_start});
            anchors.push(Anchor {ref_id, ref_start: ref_end - index.k(), query_start: query_end - index.k()});
            min_length_diff = length_diff;
        }
    }
}

fn add_to_anchors_partial(
    anchors: &mut Vec<Anchor>,
    query_start: usize,
    index: &StrobemerIndex,
    position: usize,
) {
    let hash = index.get_hash_partial(position);
    for pos in position..index.randstrobes.len() {
        if index.get_hash_partial(pos) != hash {
            break;
        }
        let randstrobe = &index.randstrobes[pos];
        let ref_id = randstrobe.reference_index();
        let (ref_start, ref_end) = index.strobe_extent_partial(pos);

        anchors.push(Anchor { ref_id, ref_start, query_start });
    }
}

fn hits_to_anchors(
    hits: &Vec<Hit>,
    index: &StrobemerIndex,
) -> Vec<Anchor> {
    let mut anchors = vec![];
    for hit in hits {
        if hit.is_filtered {
            continue;
        }
        if hit.is_partial {
            add_to_anchors_partial(&mut anchors, hit.query_start, index, hit.position);
        } else {
            add_to_anchors_full(&mut anchors, hit.query_start, hit.query_end, index, hit.position);
        }
    }

    anchors
}

fn extract_chains_from_dp(
    anchors: &[Anchor],
    dp: &[f32],
    predecessors: &[usize],
    best_score: f32,
    k: usize,
    is_revcomp: bool,
    chains: &mut Vec<Nam>,
    chaining_parameters: &ChainingParameters,
) {
    let n = anchors.len();
    let valid_score = best_score * chaining_parameters.valid_score_threshold;

    let mut candidates = vec![];
    for i in 0..n {
        if dp[i] >= valid_score {
            candidates.push((i, dp[i]));
        }
    }

    candidates.sort_by(|a, b| b.1.total_cmp(&a.1));

    let mut used = vec![false; n];
    for (i, score) in candidates {
        if used[i] {
            continue;
        }

        let mut j = i;
        let mut c = 1;
        let mut overlaps = false;

        while predecessors[j] != usize::MAX {
            j = predecessors[j];
            if used[j] {
                overlaps = true;
                break;
            }
            used[j] = true;
            c += 1;
        }

        if overlaps {
            continue;
        }

        let first = &anchors[j];
        let last = &anchors[i];
        chains.push(
            Nam {
                nam_id: chains.len(),
                query_start: first.query_start,
                query_end: last.query_start + k,
                query_prev_match_startpos: usize::MAX,
                ref_start: first.ref_start,
                ref_end: last.ref_start + k,
                ref_prev_match_startpos: usize::MAX,
                n_matches: c,
                ref_id: last.ref_id,
                score: score + c as f32 * chaining_parameters.matches_weight,
                is_revcomp,
            }
        );
    }
}
