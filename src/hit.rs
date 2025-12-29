use std::ops;

use log;
use log::trace;
use log::Level;

use crate::index::StrobemerIndex;
use crate::mapper::QueryRandstrobe;
use crate::mcsstrategy::McsStrategy;

#[derive(Debug)]
pub struct Hit {
    pub query_start: usize,
    pub query_end: usize,
    pub position: usize,
    pub hash_revcomp: u64,
    pub is_partial: bool,
    pub is_filtered: bool,
}

/// Aggregate statistics resulting from looking up all strobemers of a single
/// query
#[derive(Clone, Debug)]
pub struct HitsDetails {
    /// full hit not found
    pub full_not_found: usize,

    /// found hit found but filtered
    pub full_filtered: usize,

    /// found hit found and not filtered
    pub full_found: usize,

    pub partial_not_found: usize,
    pub partial_filtered: usize,
    pub partial_found: usize,

    pub rescued: usize,
}

impl Default for HitsDetails {
    fn default() -> Self {
        HitsDetails {
            full_not_found: 0,
            full_filtered: 0,
            full_found: 0,
            partial_not_found: 0,
            partial_filtered: 0,
            partial_found: 0,
            rescued: 0,
        }
    }
}

impl ops::AddAssign<HitsDetails> for HitsDetails {
    fn add_assign(&mut self, rhs: HitsDetails) {
        self.full_not_found += rhs.full_not_found;
        self.full_filtered += rhs.full_filtered;
        self.full_found += rhs.full_found;
        self.partial_not_found += rhs.partial_not_found;
        self.partial_filtered += rhs.partial_filtered;
        self.partial_found += rhs.partial_found;
        self.rescued += rhs.rescued;
    }
}

impl HitsDetails {
    pub fn total_hits(&self) -> usize {
        self.partial_filtered + self.partial_found + self.full_filtered + self.full_found
    }

    pub fn total_filtered(&self) -> usize {
        self.partial_filtered + self.full_filtered
    }

    pub fn total_found(&self) -> usize {
        self.full_found + self.partial_found
    }

    /// Used as a heuristic to compare the two orientations of a query
    pub fn is_better_than(&self, other: &HitsDetails) -> bool {
        let total = self.full_found + self.full_filtered;
        let other_total = other.full_found + other.full_filtered;

        total >= other_total * 2 && total > other_total + 5
    }

}

/// Finds least frequent hits in a portion of the hits vector and set their
/// 'is_filtered' attribute to false (thus "rescuing" them).
///
/// If `threshold` is set to a nonzero value, hits with a larger count are
/// ignored.
///
/// Returns the number of hits that were rescued.
fn rescue_least_frequent(
    index: &StrobemerIndex,
    hits: &mut [Hit],
    start: usize,
    end: usize,
    to_rescue: usize,
    rescue_threshold: Option<usize>,
) -> usize {
    let mut rescued = 0;

    // Index and hit count
    let mut hit_counts = vec![];
    for i in start..end {
        let cnt =
            if hits[i].is_partial {
                index.get_count_partial(hits[i].position)
            } else {
                index.get_count_full(hits[i].position, hits[i].hash_revcomp)
            };
        if rescue_threshold.is_none() || cnt <= rescue_threshold.unwrap() {
            hit_counts.push((i, cnt));
        }
    }

    // Sort by count ascending
    hit_counts.sort_by_key(|hc| hc.1);

    // Take up to num_to_rescue lowest count
    for &(hit_index, _cnt) in hit_counts.iter().take(to_rescue) {
        rescued += hits[hit_index].is_filtered as usize;
        hits[hit_index].is_filtered = false;
    }

    rescued
}

/// Find all hits for a query, using the requested MCS strategy.
/// Repetitive hits are included. `is_filtered` is set to false for all hits.
fn find_all_hits(
    query_randstrobes: &[QueryRandstrobe], index: &StrobemerIndex, filter_cutoff: usize, mcs_strategy: McsStrategy,
) -> (HitsDetails, bool, Vec<Hit>) {

    let mut hits = Vec::with_capacity(query_randstrobes.len());
    let mut sorting_needed = mcs_strategy == McsStrategy::Always || mcs_strategy == McsStrategy::FirstStrobe;
    let mut hits_details = HitsDetails::default();

    if mcs_strategy != McsStrategy::FirstStrobe {
        for randstrobe in query_randstrobes {
            if let Some(position) = index.get_full(randstrobe.hash) {
                let is_filtered = index.is_too_frequent(position, filter_cutoff, randstrobe.hash_revcomp);
                if is_filtered {
                    hits_details.full_filtered += 1;
                } else {
                    hits_details.full_found += 1;
                }

                let hit = Hit { position, query_start: randstrobe.start, query_end: randstrobe.end, is_partial: false, is_filtered, hash_revcomp: randstrobe.hash_revcomp };
                hits.push(hit);
            } else {
                hits_details.full_not_found += 1;
                if mcs_strategy == McsStrategy::Always {
                    if let Some(position) = index.get_partial(randstrobe.hash) {
                        let is_filtered = index.is_too_frequent_partial(position, filter_cutoff, randstrobe.hash_revcomp);
                        if is_filtered {
                            hits_details.partial_filtered += 1;
                        } else {
                            hits_details.partial_found += 1;
                        }
                        let hit = Hit { position, query_start: randstrobe.start, query_end: randstrobe.start + index.k(), is_partial: true, is_filtered, hash_revcomp: randstrobe.hash_revcomp };
                        hits.push(hit);
                    } else {
                        hits_details.partial_not_found += 1;
                    }
                }
            }
        }
    }
    if mcs_strategy == McsStrategy::Always {
        debug_assert!(hits_details.full_not_found == hits_details.partial_not_found + hits_details.partial_filtered + hits_details.partial_found);
    }

    // Rescue using partial hits even in non-MCS mode
    if mcs_strategy == McsStrategy::FirstStrobe || (mcs_strategy == McsStrategy::Rescue && hits_details.full_filtered + hits_details.full_found == 0) {
        for randstrobe in query_randstrobes {
            if let Some(position) = index.get_partial(randstrobe.hash) {
                let is_filtered = index.is_too_frequent_partial(position, filter_cutoff, randstrobe.hash_revcomp);
                if is_filtered {
                    hits_details.partial_filtered += 1;
                } else {
                    hits_details.partial_found += 1;
                }
                let hit = Hit { position, query_start: randstrobe.start, query_end: randstrobe.start + index.k(), is_partial: true, is_filtered, hash_revcomp: randstrobe.hash_revcomp };
                hits.push(hit);
            } else {
                hits_details.partial_not_found += 1;
            }
        }
        sorting_needed = true;
    }

    (hits_details, sorting_needed, hits)
}

/// Rescue seeds from filtered regions that have a given minimum length (in
/// nucleotides)
///
/// All stretches of consecutive filtered hits are considered. For any stretch
/// that is longer than rescue_distance, rescue_least_frequent is called.
fn rescue_all_least_frequent(
    index: &StrobemerIndex,
    hits: &mut [Hit],
    rescue_distance: usize,
) -> usize {
    let mut last_unfiltered_start = 0;
    let mut first_filtered = 0;
    let mut intervals = vec![];

    for (i, hit) in hits.iter().enumerate() {
        if hit.is_filtered {
            continue;
        }
        if hit.query_start > last_unfiltered_start + rescue_distance {
            let to_rescue = (hit.query_start - last_unfiltered_start) / rescue_distance;
            intervals.push((first_filtered, i, to_rescue));
        }
        last_unfiltered_start = hit.query_start;
        first_filtered = i + 1;
    }
    // Ensure we consider the end as well
    if let Some(last_hit) = hits.last() {
        if last_hit.query_start - last_unfiltered_start > rescue_distance {
            let to_rescue = (last_hit.query_start - last_unfiltered_start) / rescue_distance;
            intervals.push((first_filtered, hits.len(), to_rescue));
        }
    }

    let mut rescued = 0;
    for (start, end, to_rescue) in intervals {
        rescued += rescue_least_frequent(index, hits, start, end, to_rescue, None);
    }

    rescued
}

/// Find a query’s hits
pub fn find_hits(
    query_randstrobes: &[QueryRandstrobe],
    index: &StrobemerIndex,
    mcs_strategy: McsStrategy,
    filter_cutoff: usize,
    rescue_distance: usize,
) -> (HitsDetails, bool, Vec<Hit>) {
    let (mut details, sorting_needed, mut hits) = find_all_hits(query_randstrobes, index, filter_cutoff, mcs_strategy);

    let total_hits = details.total_hits();
    let nonrepetitive_hits = details.total_found();
    let nonrepetitive_fraction = if total_hits > 0 { nonrepetitive_hits as f32 / total_hits as f32 } else { 1.0 };

    // rescue distance 0 disables both global and local rescue
    if rescue_distance > 0 {
        if nonrepetitive_fraction < 0.7 {
            // "global" rescue
            let n = hits.len();
            details.rescued += rescue_least_frequent(index, &mut hits, 0, n, 5, Some(1000));
        }
        // "local" rescue
        details.rescued += rescue_all_least_frequent(index, &mut hits, rescue_distance);
    }

    if log::log_enabled!(Level::Trace) {
        trace!("Found {} hits ({} of those rescued):", hits.len(), details.rescued);
        trace!("querypos count (p=partial, F=filtered)");
        for hit in &hits {
            let cnt =
                if hit.is_partial {
                    index.get_count_partial(hit.position)
                } else {
                    index.get_count_full(hit.position, hit.hash_revcomp)
                };
            trace!("{:6} {}{:6} {}", hit.query_start, if hit.is_partial { "p" } else { " " }, cnt, if hit.is_filtered { "F" } else { " " });
        }
    }

    (details, sorting_needed, hits)
}
