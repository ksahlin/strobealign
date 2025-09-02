#include "hits.hpp"

/*
 * Find a query’s hits, ignoring randstrobes that occur too often in the
 * reference (have a count above filter_cutoff).
 *
 * Return the fraction of nonrepetitive hits (those not above the filter_cutoff threshold)
 */
std::tuple<int, int, bool, std::vector<Hit>> find_hits(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    McsStrategy mcs_strategy
) {
    // If we produce matches in sorted order, then merge_matches_into_nams()
    // does not have to re-sort
    bool sorting_needed{mcs_strategy == McsStrategy::Always || mcs_strategy == McsStrategy::FirstStrobe};
    std::vector<Hit> hits;
    int total_hits = 0;
    int partial_hits = 0;

    if (mcs_strategy == McsStrategy::FirstStrobe) {
        for (const auto &q : query_randstrobes) {
            size_t partial_position = index.find_partial(q.hash);
            if (partial_position != index.end()) {
                partial_hits++;
                if (index.is_partial_filtered(partial_position, q.hash_revcomp)) {
                    continue;
                }
                hits.push_back(Hit{partial_position, q.start, q.start + index.k(), true});
            }
        }

        return {total_hits, partial_hits, sorting_needed, hits};
    }

    for (const auto &q : query_randstrobes) {
        size_t position = index.find_full(q.hash);
        if (position != index.end()) {
            total_hits++;
            if (index.is_filtered(position, q.hash_revcomp)) {
                continue;
            }
            hits.push_back(Hit{position, q.start, q.end, false});
        } else if (mcs_strategy == McsStrategy::Always) {
            size_t partial_pos = index.find_partial(q.hash);
            if (partial_pos != index.end()) {
                total_hits++;
                if (index.is_partial_filtered(partial_pos, q.hash_revcomp)) {
                    continue;
                }
                partial_hits++;
                hits.push_back(Hit{partial_pos, q.start, q.start + index.k(), true});
            }
        }
    }

    // Rescue using partial hits
    if (mcs_strategy == McsStrategy::Rescue && total_hits == 0) {
        for (const auto &q : query_randstrobes) {
            size_t partial_pos = index.find_partial(q.hash);
            if (partial_pos != index.end()) {
                total_hits++;
                if (index.is_partial_filtered(partial_pos, q.hash_revcomp)) {
                    continue;
                }
                partial_hits++;
                hits.push_back(Hit{partial_pos, q.start, q.start + index.k(), true});
            }
        }
        sorting_needed = true;
    }

    return {total_hits, partial_hits, sorting_needed, hits};
}
