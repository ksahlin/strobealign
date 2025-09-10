#include "hits.hpp"

/*
 * Find a query’s hits, ignoring randstrobes that occur too often in the
 * reference (have a count above filter_cutoff).
 */
std::tuple<HitsDetails, int, bool, std::vector<Hit>> find_hits(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    McsStrategy mcs_strategy
) {
    // If we produce matches in sorted order, then merge_matches_into_nams()
    // does not have to re-sort
    bool sorting_needed{mcs_strategy == McsStrategy::Always || mcs_strategy == McsStrategy::FirstStrobe};
    std::vector<Hit> hits;
    HitsDetails details;

    if (mcs_strategy == McsStrategy::FirstStrobe) {
        for (const auto &q : query_randstrobes) {
            size_t partial_position = index.find_partial(q.hash);
            if (partial_position != index.end()) {
                if (index.is_partial_filtered(partial_position, q.hash_revcomp)) {
                    details.partial_filtered++;
                    continue;
                }
                details.partial_found++;
                hits.push_back(Hit{partial_position, q.start, q.start + index.k(), true});
            } else {
                details.partial_not_found++;
            }
        }
        uint partial_hits = details.partial_filtered + details.partial_found;

        return {details, partial_hits, sorting_needed, hits};
    }

    for (const auto &q : query_randstrobes) {
        size_t position = index.find_full(q.hash);
        if (position != index.end()) {
            if (index.is_filtered(position, q.hash_revcomp)) {
                details.full_filtered++;
                continue;
            }
            details.full_found++;
            hits.push_back(Hit{position, q.start, q.end, false});
        } else {
            details.full_not_found++;
            if (mcs_strategy == McsStrategy::Always) {
                size_t partial_pos = index.find_partial(q.hash);
                if (partial_pos != index.end()) {
                    if (index.is_partial_filtered(partial_pos, q.hash_revcomp)) {
                        details.partial_filtered++;
                        continue;
                    }
                    details.partial_found++;
                    hits.push_back(Hit{partial_pos, q.start, q.start + index.k(), true});
                } else {
                    details.partial_not_found++;
                }
            }
        }
    }
    if (mcs_strategy == McsStrategy::Always) {
        assert(details.full_not_found == details.partial_not_found + details.partial_filtered + details.partial_found);
    }

    // Rescue using partial hits
    if (mcs_strategy == McsStrategy::Rescue && details.full_filtered + details.full_found == 0) {
        for (const auto &q : query_randstrobes) {
            size_t partial_pos = index.find_partial(q.hash);
            if (partial_pos != index.end()) {
                if (index.is_partial_filtered(partial_pos, q.hash_revcomp)) {
                    details.partial_filtered++;
                    continue;
                }
                details.partial_found++;
                hits.push_back(Hit{partial_pos, q.start, q.start + index.k(), true});
            } else {
                details.partial_not_found++;
            }
        }
        sorting_needed = true;
    }

    uint total_hits = details.partial_filtered + details.partial_found + details.full_filtered + details.full_found;

    return {details, total_hits, sorting_needed, hits};
}
