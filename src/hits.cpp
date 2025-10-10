#include "hits.hpp"

/*
 * Go back over randstrobes that were previously filtered and add the least
 * frequent ones to the hits vector.
 *
 * Return the number of hits that were rescued
 */
uint rescue_least_frequent(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    size_t start,
    size_t end,
    std::vector<Hit>& hits,
    int distance,
    int L
) {
    struct Candidate {
        int count;
        const QueryRandstrobe* query_randstrobe;
        size_t position;
        bool is_partial;
    };
    uint rescued{0};
    size_t num_to_rescue = distance / L;
    std::vector<Candidate> candidates;

    for (size_t i = start; i < end; ++i) { // seeds between that first and last unfiltered seeds
        const auto &q = query_randstrobes[i];
        size_t position = index.find_full(q.hash);
        if (position != index.end()) {
            int cnt = index.get_count_full(position);
            candidates.emplace_back(Candidate{cnt, &q, position, false});
        } else {
            size_t partial_position = index.find_partial(q.hash);
            if (partial_position != index.end()) {
                int cnt = index.get_count_partial(partial_position);
                candidates.emplace_back(Candidate{cnt, &q, partial_position, true});
            }
        }
    }

    // Sort by count ascending
    std::sort(candidates.begin(), candidates.end(), [](auto &a, auto &b) { return a.count < b.count; });

    // Take up to num_to_rescue lowest count
    for (size_t i = 0; i < std::min(num_to_rescue, candidates.size()); ++i) {
        auto [cnt, q, pos, is_partial] = candidates[i];
        if (is_partial)
            hits.push_back(Hit{pos, q->start, q->start + index.k(), is_partial});
        else{
            hits.push_back(Hit{pos, q->start, q->end, is_partial});
        }
        rescued++;
    }

    return rescued;
}

/*
 * Find a query’s hits, ignoring randstrobes that occur too often in the
 * reference (have a count above filter_cutoff).
 */
std::tuple<HitsDetails, bool, std::vector<Hit>> find_hits(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    McsStrategy mcs_strategy,
    int rescue_threshold
) {
    // If we produce matches in sorted order, then merge_matches_into_nams()
    // does not have to re-sort
    bool sorting_needed{mcs_strategy == McsStrategy::Always || mcs_strategy == McsStrategy::FirstStrobe};
    std::vector<Hit> hits;
    HitsDetails details;

    // Rescue threshold: If all hits over a region of this length (in nucleotides)
    // are initially filtered out, we go back and add the least frequent of them
    const int L = rescue_threshold; // threshold for rescue
    int n_filtered = 0; // For checking that there are eligible seeds in window L
    int last_unfiltered_start = 0;
    size_t first_filtered = 0;

    if (mcs_strategy != McsStrategy::FirstStrobe) {
        for (size_t i = 0; i < query_randstrobes.size(); ++i) {
            const auto &q = query_randstrobes[i];
            size_t position = index.find_full(q.hash);
            if (position != index.end()) {
                if (index.is_filtered(position, q.hash_revcomp)) {
                    details.full_filtered++;
                    n_filtered++;
                    continue;
                }
                details.full_found++;
                if ((q.start - last_unfiltered_start > L) && n_filtered > 0) {
                    details.rescued += rescue_least_frequent(query_randstrobes, index, first_filtered, i, hits, q.start - last_unfiltered_start, L);
                }
                last_unfiltered_start = q.start;
                first_filtered = i + 1;
                n_filtered = 0;
                hits.push_back(Hit{position, q.start, q.end, false});
            } else {
                details.full_not_found++;
                if (mcs_strategy == McsStrategy::Always) {
                    size_t partial_pos = index.find_partial(q.hash);
                    if (partial_pos != index.end()) {
                        if (index.is_partial_filtered(partial_pos, q.hash_revcomp)) {
                            details.partial_filtered++;
                            n_filtered++;
                            continue;
                        }
                        details.partial_found++;
                        if ((q.start - last_unfiltered_start > L) && n_filtered > 0) {
                            details.rescued += rescue_least_frequent(query_randstrobes, index, first_filtered, i, hits, q.start - last_unfiltered_start, L);
                        }
                        last_unfiltered_start = q.start;
                        first_filtered = i + 1;
                        n_filtered = 0;
                        hits.push_back(Hit{partial_pos, q.start, q.start + index.k(), true});
                    } else {
                        details.partial_not_found++;
                    }
                }
            }
        }
        if (!query_randstrobes.empty() && query_randstrobes.back().start - last_unfiltered_start > L && n_filtered > 0) { // End case we have not sampled the end
            details.rescued += rescue_least_frequent(query_randstrobes, index, first_filtered, query_randstrobes.size(), hits, query_randstrobes.back().start - last_unfiltered_start, L);
        }
    }

    // Only partial lookups
    if (
        mcs_strategy == McsStrategy::FirstStrobe
        || (mcs_strategy == McsStrategy::Rescue && details.full_filtered + details.full_found == 0)
    ) {
        n_filtered = 0;
        last_unfiltered_start = 0;
        first_filtered = 0;
        for (size_t i = 0; i < query_randstrobes.size(); ++i) {
            const auto &q = query_randstrobes[i];
            size_t partial_pos = index.find_partial(q.hash);
            if (partial_pos != index.end()) {
                if (index.is_partial_filtered(partial_pos, q.hash_revcomp)) {
                    details.partial_filtered++;
                    n_filtered++;
                    continue;
                }
                details.partial_found++;
                if ((q.start - last_unfiltered_start > L) && n_filtered > 0) {
                    details.rescued += rescue_least_frequent(query_randstrobes, index, first_filtered, i, hits, q.start - last_unfiltered_start, L);
                }
                last_unfiltered_start = q.start;
                first_filtered = i + 1;
                n_filtered = 0;
                hits.push_back(Hit{partial_pos, q.start, q.start + index.k(), true});
            } else {
                details.partial_not_found++;
            }
        }
        sorting_needed = true;
    }

    if (mcs_strategy == McsStrategy::Always) {
        assert(details.full_not_found == details.partial_not_found + details.partial_filtered + details.partial_found);
    }


    return {details, sorting_needed, hits};
}
