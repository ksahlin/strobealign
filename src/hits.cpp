#include "hits.hpp"

void add_seeds(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    size_t start,
    size_t end,
    std::vector<Hit>& hits,
    int distance,
    int L
) {
    int num_to_rescue = distance / L;
    std::vector<std::tuple<int, const QueryRandstrobe*, size_t, bool>> candidates;

    for (size_t i = start+1; i < end-1; ++i) { // seeds between that first and last unfiltered seeds
        const auto &q = query_randstrobes[i];
        size_t position = index.find_full(q.hash);
        if (position != index.end()) {
            int cnt = index.get_count_full(position);
            candidates.emplace_back(cnt, &q, position, false);
        } else {
            size_t partial_position = index.find_partial(q.hash);
            if (partial_position != index.end()) {
                int cnt = index.get_count_partial(partial_position);
                candidates.emplace_back(cnt, &q, partial_position, true);
            }
    }
    }


    // Sort by count ascending
    std::sort(candidates.begin(), candidates.end(),
                [](auto &a, auto &b) { return std::get<0>(a) < std::get<0>(b); });

    // Take up to num_to_rescue lowest count
    for (int i = 0; i < num_to_rescue && i < (int)candidates.size(); ++i) {
        auto [cnt, q, pos, is_partial] = candidates[i];
        if (is_partial)
            hits.push_back(Hit{pos, q->start, q->start + index.k(), is_partial});
        else{
            hits.push_back(Hit{pos, q->start, q->end, is_partial});
        }
    }

}

/*
 * Find a query’s hits, ignoring randstrobes that occur too often in the
 * reference (have a count above filter_cutoff).
 */
std::tuple<HitsDetails, bool, std::vector<Hit>> find_hits(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    McsStrategy mcs_strategy
) {
    // If we produce matches in sorted order, then merge_matches_into_nams()
    // does not have to re-sort
    bool sorting_needed{mcs_strategy == McsStrategy::Always || mcs_strategy == McsStrategy::FirstStrobe};
    std::vector<Hit> hits;
    HitsDetails details;

    const int L = 100; // threshold for rescue
    int last_unfiltered_start = 0;
    int filter_start = 0;

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

        return {details, sorting_needed, hits};
    }

    for (size_t i = 0; i < query_randstrobes.size(); ++i) {
        const auto &q = query_randstrobes[i];
        size_t position = index.find_full(q.hash);
        if (position != index.end()) {
            if (index.is_filtered(position, q.hash_revcomp)) {
                details.full_filtered++;
                continue;
            }
            details.full_found++;
            if ((q.start - last_unfiltered_start) > L) { // Rescue
                add_seeds(query_randstrobes, index, filter_start, i, hits, q.start - last_unfiltered_start, L);
            }
            last_unfiltered_start = q.start;
            filter_start = i;
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
                    if ((q.start - last_unfiltered_start) > L) { // Rescue
                        add_seeds(query_randstrobes, index, filter_start, i, hits, q.start - last_unfiltered_start, L);
                    }
                    last_unfiltered_start = q.start;
                    filter_start = i;
                    hits.push_back(Hit{partial_pos, q.start, q.start + index.k(), true});                 
                } else {
                    details.partial_not_found++;
                }
            }
        }
    }

    if (!query_randstrobes.empty() && query_randstrobes.back().start - last_unfiltered_start > L) { // End case we have not sampled the end
        add_seeds(query_randstrobes, index, filter_start, query_randstrobes.size(), hits, query_randstrobes.back().start - last_unfiltered_start, L);
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

    assert(details.full_found + details.partial_found == hits.size());

    return {details, sorting_needed, hits};
}
