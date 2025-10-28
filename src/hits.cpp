#include "hits.hpp"

struct CandidateHit {
    size_t position;
    randstrobe_hash_t hash_revcomp;
    size_t query_start;
    size_t query_end;
    bool is_partial;
    bool is_filtered;
};

/*
 * Find least frequent hits in a portion of the hits vector and set their
 * 'is_filtered' attribute to false (thus "rescuing" them).
 *
 * Return the number of hits that were rescued.
 */
uint rescue_least_frequent(
    const StrobemerIndex& index,
    std::vector<CandidateHit>& hits,
    size_t start,
    size_t end,
    int distance,
    int L
) {
    uint rescued{0};
    size_t num_to_rescue = distance / L;

    // Index and hit count
    std::vector<std::pair<size_t, uint>> hit_counts;
    for (size_t i = start; i < end; ++i) {
        uint cnt;
        if (hits[i].is_partial) {
            cnt = index.get_count_partial(hits[i].position);
        } else {
            cnt = index.get_count_full(hits[i].position, hits[i].hash_revcomp);
        }
        hit_counts.push_back({i, cnt});
    }

    // Sort by count ascending
    std::sort(hit_counts.begin(), hit_counts.end(), [](auto &a, auto &b) { return a.second < b.second; });

    // Take up to num_to_rescue lowest count
    for (size_t i = 0; i < std::min(num_to_rescue, hit_counts.size()); ++i) {
        auto hit_index = hit_counts[i].first;
        rescued += hits[hit_index].is_filtered;
        hits[hit_index].is_filtered = false;
    }

    return rescued;
}

/*
 * Find all hits using the requested MCS strategy. Keep repetitive hits.
 */
std::tuple<std::vector<CandidateHit>, HitsDetails, bool> find_candidate_hits(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    McsStrategy mcs_strategy
) {
    // If we produce matches in sorted order, then merge_matches_into_nams()
    // does not have to re-sort
    bool sorting_needed{mcs_strategy == McsStrategy::Always || mcs_strategy == McsStrategy::FirstStrobe};
    HitsDetails details;
    std::vector<CandidateHit> hits;
    if (mcs_strategy != McsStrategy::FirstStrobe) {
        for (size_t i = 0; i < query_randstrobes.size(); ++i) {
            const auto &q = query_randstrobes[i];
            size_t position = index.find_full(q.hash);
            if (position != index.end()) {
                bool is_filtered = index.is_filtered(position, q.hash_revcomp);
                if (is_filtered) {
                    details.full_filtered++;
                } else {
                    details.full_found++;
                }
                hits.push_back(CandidateHit{position, q.hash_revcomp, q.start, q.end, false, is_filtered});
            } else {
                details.full_not_found++;
                if (mcs_strategy == McsStrategy::Always) {
                    size_t partial_pos = index.find_partial(q.hash);
                    if (partial_pos != index.end()) {
                        bool is_filtered = index.is_partial_filtered(partial_pos, q.hash_revcomp);
                        if (is_filtered) {
                            details.partial_filtered++;
                        } else {
                            details.partial_found++;
                        }
                        hits.push_back(CandidateHit{partial_pos, q.hash_revcomp, q.start, q.start + index.k(), true, is_filtered});
                    } else {
                        details.partial_not_found++;
                    }
                }
            }
        }
    }

    if (
        mcs_strategy == McsStrategy::FirstStrobe
        || (mcs_strategy == McsStrategy::Rescue && hits.size() == 0)
    ) {
        // Only partial lookups
        for (size_t i = 0; i < query_randstrobes.size(); ++i) {
            const auto &q = query_randstrobes[i];
            size_t partial_pos = index.find_partial(q.hash);
            if (partial_pos != index.end()) {
                bool is_filtered = index.is_partial_filtered(partial_pos, q.hash_revcomp);
                if (is_filtered) {
                    details.partial_filtered++;
                } else {
                    details.partial_found++;
                }
                hits.push_back(CandidateHit{partial_pos, q.hash_revcomp, q.start, q.start + index.k(), true, is_filtered});
            } else {
                details.partial_not_found++;
            }
        }
        sorting_needed = true;
    }

    if (mcs_strategy == McsStrategy::Always) {
        assert(details.full_not_found == details.partial_not_found + details.partial_filtered + details.partial_found);
    }

    return {hits, details, sorting_needed};
}


uint rescue_all_least_frequent(
    const StrobemerIndex& index,
    std::vector<CandidateHit>& candidates,
    uint rescue_threshold
) {
    // Rescue threshold: If all hits over a region of this length (in nucleotides)
    // are initially filtered out, we go back and add the least frequent of them
    const uint L = rescue_threshold; // threshold for rescue
    int last_unfiltered_start = 0;
    size_t first_filtered = 0;
    uint rescued = 0;

    for (size_t i = 0; i < candidates.size(); ++i) {
        const auto &hit = candidates[i];

        if (hit.is_filtered) {
            continue;
        }
        if (hit.query_start > last_unfiltered_start + L) {
            rescued += rescue_least_frequent(index, candidates, first_filtered, i, hit.query_start - last_unfiltered_start, L);
        }
        last_unfiltered_start = hit.query_start;
        first_filtered = i + 1;
    }
    if (!candidates.empty() && candidates.back().query_start - last_unfiltered_start > L) { // End case we have not sampled the end
        rescued += rescue_least_frequent(index, candidates, first_filtered, candidates.size(), candidates.back().query_start - last_unfiltered_start, L);
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
    auto [candidates, details, sorting_needed] = find_candidate_hits(query_randstrobes, index, mcs_strategy);
    if (rescue_threshold != -1) {
        details.rescued += rescue_all_least_frequent(index, candidates, rescue_threshold);
    }
    std::vector<Hit> hits;
    for (const auto& candidate : candidates) {
        if (!candidate.is_filtered) {
            hits.push_back(Hit{candidate.position, candidate.query_start, candidate.query_end, candidate.is_partial});
        }
    }

    return {details, sorting_needed, hits};
}
