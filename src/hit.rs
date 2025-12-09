use std::ops;
use crate::index::StrobemerIndex;
use crate::mapper::QueryRandstrobe;
use crate::mcsstrategy::McsStrategy;

#[derive(Debug)]
pub struct Hit {
    pub query_start: usize,
    pub query_end: usize,
    pub position: usize,
    pub is_partial: bool,
}

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
    }
}

impl HitsDetails {
    pub fn total_hits(&self) -> usize {
        self.partial_filtered + self.partial_found + self.full_filtered + self.full_found
    }
}

/// Find a query’s hits, ignoring randstrobes that occur too often in the
/// reference (have a count above filter_cutoff).
pub fn find_hits(
    query_randstrobes: &[QueryRandstrobe], index: &StrobemerIndex, filter_cutoff: usize, mcs_strategy: McsStrategy
) -> (HitsDetails, bool, Vec<Hit>) {

    let mut hits = vec![];
    let mut sorting_needed = mcs_strategy == McsStrategy::Always || mcs_strategy == McsStrategy::FirstStrobe;
    let mut hits_details = HitsDetails::default();

    if mcs_strategy != McsStrategy::FirstStrobe {
        for randstrobe in query_randstrobes {
            if let Some(position) = index.get_full(randstrobe.hash) {
                if index.is_too_frequent(position, filter_cutoff, randstrobe.hash_revcomp) {
                    hits_details.full_filtered += 1;
                    continue;
                }
                hits_details.full_found += 1;

                let hit = Hit { position, query_start: randstrobe.start, query_end: randstrobe.end, is_partial: false };
                hits.push(hit);
            } else {
                hits_details.full_not_found += 1;
                if mcs_strategy == McsStrategy::Always {
                    if let Some(position) = index.get_partial(randstrobe.hash) {
                        if index.is_too_frequent_partial(position, filter_cutoff, randstrobe.hash_revcomp) {
                            hits_details.partial_filtered += 1;
                            continue;
                        }
                        hits_details.partial_found += 1;
                        let hit = Hit { position, query_start: randstrobe.start, query_end: randstrobe.start + index.k(), is_partial: true };
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
                if index.is_too_frequent_partial(position, filter_cutoff, randstrobe.hash_revcomp) {
                    hits_details.partial_filtered += 1;
                    continue;
                }
                hits_details.partial_found += 1;
                let hit = Hit { position, query_start: randstrobe.start, query_end: randstrobe.start + index.k(), is_partial: true };
                hits.push(hit);
            } else {
                hits_details.partial_not_found += 1;
            }
        }
        sorting_needed = true;
    }

    debug_assert!(hits_details.full_found + hits_details.partial_found == hits.len());

    (hits_details, sorting_needed, hits)
}
