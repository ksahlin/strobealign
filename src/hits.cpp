#include "hits.hpp"
#include "logger.hpp"

#include <iomanip>

static Logger& logger = Logger::get();

/*
 * Find least frequent hits in a portion of the hits vector and set their
 * 'is_filtered' attribute to false (thus "rescuing" them).
 *
 * Return the number of hits that were rescued.
 */
uint rescue_least_frequent(
    const StrobemerIndex& index,
    std::vector<Hit>& hits,
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
 * Find all hits, including repetitive ones, using the requested MCS strategy.
 */
std::tuple<std::vector<Hit>, HitsDetails, bool> find_all_hits(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    McsStrategy mcs_strategy
) {
    // If we produce matches in sorted order, then merge_matches_into_nams()
    // does not have to re-sort
    bool sorting_needed{mcs_strategy == McsStrategy::Always || mcs_strategy == McsStrategy::FirstStrobe};
    HitsDetails details;
    std::vector<Hit> hits;
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
                hits.push_back(Hit{position, q.hash_revcomp, q.start, q.end, false, is_filtered});
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
                        hits.push_back(Hit{partial_pos, q.hash_revcomp, q.start, q.start + index.k(), true, is_filtered});
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
                hits.push_back(Hit{partial_pos, q.hash_revcomp, q.start, q.start + index.k(), true, is_filtered});
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
    std::vector<Hit>& hits,
    uint rescue_threshold
) {
    // Rescue threshold: If all hits over a region of this length (in nucleotides)
    // are initially filtered out, we go back and add the least frequent of them
    const uint L = rescue_threshold; // threshold for rescue
    int last_unfiltered_start = 0;
    size_t first_filtered = 0;
    uint rescued = 0;

    for (size_t i = 0; i < hits.size(); ++i) {
        const auto &hit = hits[i];

        if (hit.is_filtered) {
            continue;
        }
        if (hit.query_start > last_unfiltered_start + L) {
            rescued += rescue_least_frequent(index, hits, first_filtered, i, hit.query_start - last_unfiltered_start, L);
        }
        last_unfiltered_start = hit.query_start;
        first_filtered = i + 1;
    }
    if (!hits.empty() && hits.back().query_start - last_unfiltered_start > L) { // End case we have not sampled the end
        rescued += rescue_least_frequent(index, hits, first_filtered, hits.size(), hits.back().query_start - last_unfiltered_start, L);
    }

    return rescued;
}

/*
 * Find a query’s hits
 */
std::tuple<HitsDetails, bool, std::vector<Hit>> find_hits(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    McsStrategy mcs_strategy,
    int rescue_threshold
) {
    auto [hits, details, sorting_needed] = find_all_hits(query_randstrobes, index, mcs_strategy);
    details.rescued += rescue_all_least_frequent(index, hits, rescue_threshold);

    if (logger.level() <= LOG_TRACE) {
        logger.trace() << "Found " << hits.size() << " hits (" << details.rescued << " of those rescued):\n";
        logger.trace() << "querypos count (p=partial, F=filtered)\n";
        for (const auto& hit : hits) {
            int cnt;
            if (hit.is_partial) {
                cnt = index.get_count_partial(hit.position);
            } else {
                cnt = index.get_count_full(hit.position, hit.hash_revcomp);
            }
            logger.trace()
                << std::setw(6) << hit.query_start
                << (hit.is_partial ? " p" : "  ")
                << std::setw(6) << cnt
                << (hit.is_filtered ? " F" : "  ")
                << "\n";
        }
    }

    return {details, sorting_needed, hits};
}


std::ostream& operator<<(std::ostream& os, const HitsDetails& details) {
    os  << "HitsDetails("
        << "full_not_found=" << details.full_not_found
        << ", full_filtered=" << details.full_filtered
        << ", full_found=" << details.full_found
        << ", partial_not_found=" << details.partial_not_found
        << ", partial_filtered=" << details.partial_filtered
        << ", partial_found=" << details.partial_found
        << ", rescued=" << details.rescued
        << ", filtered_nucleotides=" << details.filtered_nucleotides
        << ")";
    return os;
}
