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

/// Find a query’s hits, ignoring randstrobes that occur too often in the
/// reference (have a count above filter_cutoff).
///
/// Return the fraction of nonrepetitive hits (those not above the filter_cutoff threshold)
///
pub fn find_hits(
    query_randstrobes: &[QueryRandstrobe], index: &StrobemerIndex, filter_cutoff: usize, mcs_strategy: McsStrategy
) -> (usize, usize, bool, Vec<Hit>) {

    let mut hits = vec![];
    let mut sorting_needed = mcs_strategy == McsStrategy::Always || mcs_strategy == McsStrategy::FirstStrobe;
    let mut total_hits = 0;
    let mut partial_hits = 0;

    if mcs_strategy == McsStrategy::FirstStrobe {
        for randstrobe in query_randstrobes {
            if let Some(position) = index.get_partial(randstrobe.hash) {
                partial_hits += 1;
                if index.is_too_frequent_partial(position, filter_cutoff, randstrobe.hash_revcomp) {
                    continue;
                }

                let hit = Hit { position, query_start: randstrobe.start, query_end: randstrobe.start + index.k(), is_partial: true };
                hits.push(hit);
            }
        }

        return (partial_hits, partial_hits, sorting_needed, hits);
    }

    for randstrobe in query_randstrobes {
        if let Some(position) = index.get_full(randstrobe.hash) {
            total_hits += 1;
            if index.is_too_frequent(position, filter_cutoff, randstrobe.hash_revcomp) {
                continue;
            }

            let hit = Hit { position, query_start: randstrobe.start, query_end: randstrobe.end, is_partial: false };
            hits.push(hit);
        } else if mcs_strategy == McsStrategy::Always {
            if let Some(position) = index.get_partial(randstrobe.hash) {
                total_hits += 1;
                if index.is_too_frequent_partial(position, filter_cutoff, randstrobe.hash_revcomp) {
                    continue;
                }
                partial_hits += 1;
                let hit = Hit { position, query_start: randstrobe.start, query_end: randstrobe.start + index.k(), is_partial: true };
                hits.push(hit);
            }
        }
    }

    // Rescue using partial hits even in non-MCS mode
    if total_hits == 0 && mcs_strategy == McsStrategy::Rescue {
        for randstrobe in query_randstrobes {
            if let Some(position) = index.get_partial(randstrobe.hash) {
                total_hits += 1;
                if index.is_too_frequent_partial(position, filter_cutoff, randstrobe.hash_revcomp) {
                    continue;
                }
                partial_hits += 1;
                let hit = Hit { position, query_start: randstrobe.start, query_end: randstrobe.start + index.k(), is_partial: true };
                hits.push(hit);
            }
        }
        sorting_needed = true;
    }

    (total_hits, partial_hits, sorting_needed, hits)
}
