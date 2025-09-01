#include "nam.hpp"

bool operator==(const Match& lhs, const Match& rhs) {
    return (lhs.query_start == rhs.query_start) && (lhs.query_end == rhs.query_end) && (lhs.ref_start == rhs.ref_start) && (lhs.ref_end == rhs.ref_end);
}


namespace {

inline void add_to_matches_map_full(
    robin_hood::unordered_map<unsigned int, std::vector<Match>>& matches_map,
    int query_start,
    int query_end,
    const StrobemerIndex& index,
    size_t position
) {
    int min_diff = std::numeric_limits<int>::max();
    for (const auto hash = index.get_hash(position); index.get_hash(position) == hash; ++position) {
        int ref_start = index.get_strobe1_position(position);
        int ref_end = ref_start + index.strobe2_offset(position) + index.strobe3_offset(position) + index.k();
        int diff = std::abs((query_end - query_start) - (ref_end - ref_start));
        if (diff <= min_diff) {
            matches_map[index.reference_index(position)].push_back(
                Match{query_start, query_end, ref_start, ref_end}
            );
            min_diff = diff;
        }
    }
}

/*
This function often produces the same Match multiple times.
This happens when multiple randstrobes involve a syncmer that has a low hash
value (and which therefore gets chosen as main hash multiple times).

For example, assume we have syncmers starting at positions 14, 21, 39, 51
and randstrobes like this:

14, 21, 39, 51
14------39
    21--39
        39--51

If 39 now has the lowest hash value, it will be chosen as the main hash for all
three randstrobes. Thus, when looking up that main hash, all three randstrobes
are found. After the start coordinate of the randstrobe is converted to the
start coordinate of the single syncmer (using strobe_extent_partial()), they all
result in a hit that starts at 39 (and has length k).
*/
inline void add_to_matches_map_partial(
    robin_hood::unordered_map<unsigned int, std::vector<Match>>& matches_map,
    int query_start,
    int query_end,
    const StrobemerIndex& index,
    size_t position
) {
    for (const auto hash = index.get_main_hash(position);
         index.get_main_hash(position) == hash;
         ++position
    ) {
        auto [ref_start, ref_end] = index.strobe_extent_partial(position);
        matches_map[index.reference_index(position)].push_back(
            Match{query_start, query_end, ref_start, ref_end}
        );
    }
}

} // namespace

void merge_matches_into_nams(
    robin_hood::unordered_map<unsigned int, std::vector<Match>>& matches_map,
    int k,
    bool sort,
    bool is_revcomp,
    std::vector<Nam>& nams  // inout
) {
    for (auto &[ref_id, matches] : matches_map) {
        if (sort) {
            std::sort(matches.begin(), matches.end(), [](const Match& a, const Match& b) -> bool {
                    // first sort on query starts, then on reference starts, finally prefer full matches over partial
                    return (a.query_start < b.query_start) ||
                            ( (a.query_start == b.query_start) && (a.ref_start < b.ref_start)) ||
                            ( (a.query_start == b.query_start) && (a.ref_start == b.ref_start) && (a.query_end > b.query_end)  );
                }
            );
        }

        std::vector<Nam> open_nams;
        int prev_q_start = 0;
        auto prev_match = Match{0,0,0,0};
        for (auto &m : matches) {
            if ( (prev_match == m) || ( ((m.query_end - m.query_start) == k) && (prev_match.query_start == m.query_start) && (prev_match.ref_start == m.ref_start)) )  { // is a redundant partial match
                continue;
            }
            bool is_added = false;
            for (auto & o : open_nams) {

                // Extend NAM
                if ((o.query_prev_match_startpos <= m.query_start) && (m.query_start <= o.query_end ) && (o.ref_prev_match_startpos <= m.ref_start) && (m.ref_start <= o.ref_end) ){
                    if ( (m.query_end > o.query_end) && (m.ref_end > o.ref_end) ) {
                        o.query_end = m.query_end;
                        o.ref_end = m.ref_end;
//                        o.previous_query_start = h.query_s;
//                        o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
                        o.query_prev_match_startpos = m.query_start;
                        o.ref_prev_match_startpos = m.ref_start;
                        o.n_matches ++;
//                        o.score += (float)1/ (float)h.count;
                        is_added = true;
                        break;
                    }
                    else if ((m.query_end <= o.query_end) && (m.ref_end <= o.ref_end)) {
//                        o.previous_query_start = h.query_s;
//                        o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
                        o.query_prev_match_startpos = m.query_start;
                        o.ref_prev_match_startpos = m.ref_start;
                        o.n_matches ++;
//                        o.score += (float)1/ (float)h.count;
                        is_added = true;
                        break;
                    }
                }

            }
            // Add to open matches
            if (!is_added){
                Nam n;
                n.query_start = m.query_start;
                n.query_end = m.query_end;
                n.ref_start = m.ref_start;
                n.ref_end = m.ref_end;
                n.ref_id = ref_id;
//                n.previous_query_start = h.query_s;
//                n.previous_ref_start = h.ref_s;
                n.query_prev_match_startpos = m.query_start;
                n.ref_prev_match_startpos = m.ref_start;
                n.n_matches = 1;
                n.is_revcomp = is_revcomp;
//                n.score += (float)1 / (float)h.count;
                open_nams.push_back(n);
            }

            // Only filter if we have advanced at least k nucleotides
            if (m.query_start > prev_q_start + k) {

                // Output all NAMs from open_matches to final_nams that the current match have passed
                for (auto &n : open_nams) {
                    if (n.query_end < m.query_start) {
                        int n_max_span = std::max(n.query_span(), n.ref_span());
                        int n_min_span = std::min(n.query_span(), n.ref_span());
                        float n_score;
                        n_score = ( 2*n_min_span -  n_max_span) > 0 ? (float) (n.n_matches * ( 2*n_min_span -  n_max_span) ) : 1;   // this is really just n_matches * ( min_span - (offset_in_span) ) );
//                        n_score = n.n_hits * n.query_span();
                        n.score = n_score;
                        n.nam_id = nams.size();
                        nams.push_back(n);
                    }
                }

                // Remove all NAMs from open_matches that the current match have passed
                auto c = m.query_start;
                auto predicate = [c](decltype(open_nams)::value_type const &nam) { return nam.query_end < c; };
                open_nams.erase(std::remove_if(open_nams.begin(), open_nams.end(), predicate), open_nams.end());
                prev_q_start = m.query_start;
            }
            prev_match = m;
        }

        // Add all current open_matches to final NAMs
        for (auto &n : open_nams) {
            int n_max_span = std::max(n.query_span(), n.ref_span());
            int n_min_span = std::min(n.query_span(), n.ref_span());
            float n_score;
            n_score = ( 2*n_min_span -  n_max_span) > 0 ? (float) (n.n_matches * ( 2*n_min_span -  n_max_span) ) : 1;   // this is really just n_matches * ( min_span - (offset_in_span) ) );
//            n_score = n.n_matches * n.query_span();
            n.score = n_score;
            n.nam_id = nams.size();
            nams.push_back(n);
        }
    }
}

robin_hood::unordered_map<unsigned int, std::vector<Match>> hits_to_matches(
    const std::vector<Hit>& hits,
    const StrobemerIndex& index
) {
    robin_hood::unordered_map<unsigned int, std::vector<Match>> matches_map;
    matches_map.reserve(100);

    for (const auto& hit : hits) {
        if (hit.is_partial) {
            add_to_matches_map_partial(matches_map, hit.query_start, hit.query_end, index, hit.position);
        } else {
            add_to_matches_map_full(matches_map, hit.query_start, hit.query_end, index, hit.position);
        }
    }

    return matches_map;
}

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
    int nonrepetitive_hits = 0;
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

/*
 * Find a query’s NAMs, using also some of the randstrobes that occur more often
 * than filter_cutoff.
 *
 * Return the number of hits and the vector of NAMs.
 */
std::tuple<int, int, robin_hood::unordered_map<unsigned int, std::vector<Match>>> find_matches_rescue(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    unsigned int rescue_cutoff,
    McsStrategy mcs_strategy
) {
    struct RescueHit {
        size_t position;
        unsigned int count;
        unsigned int query_start;
        unsigned int query_end;
        bool is_partial;

        bool operator< (const RescueHit& rhs) const {
            return std::tie(count, query_start, query_end)
                < std::tie(rhs.count, rhs.query_start, rhs.query_end);
        }
    };
    robin_hood::unordered_map<unsigned int, std::vector<Match>> matches_map;
    matches_map.reserve(100);
    int n_hits = 0;
    int partial_hits = 0;
    std::vector<RescueHit> rescue_hits;
    rescue_hits.reserve(5000);
    for (auto &qr : query_randstrobes) {
        size_t position = index.find_full(qr.hash);
        if (position != index.end()) {
            unsigned int count = index.get_count_full(position);

            size_t position_revcomp = index.find_full(qr.hash_revcomp);
            if (position_revcomp != index.end()) {
                count += index.get_count_full(position_revcomp);
            }
            RescueHit rh{position, count, qr.start, qr.end, false};
            rescue_hits.push_back(rh);
        }
        else if (mcs_strategy == McsStrategy::Always) {
            size_t partial_pos = index.find_partial(qr.hash);
            if (partial_pos != index.end()) {
                unsigned int partial_count = index.get_count_partial(partial_pos);
                size_t position_revcomp = index.find_partial(qr.hash_revcomp);
                if (position_revcomp != index.end()) {
                    partial_count += index.get_count_partial(position_revcomp);
                }
                RescueHit rh{partial_pos, partial_count, qr.start, qr.start + index.k(), true};
                rescue_hits.push_back(rh);
                partial_hits++;
            }
        }
    }
    std::sort(rescue_hits.begin(), rescue_hits.end());

    int cnt = 0;
    for (auto &rh : rescue_hits) {
        if ((rh.count > rescue_cutoff && cnt >= 5) || rh.count > 1000) {
            break;
        }
        if (rh.is_partial){
            partial_hits++;
            add_to_matches_map_partial(matches_map, rh.query_start, rh.query_end, index, rh.position);
        } else{
            add_to_matches_map_full(matches_map, rh.query_start, rh.query_end, index, rh.position);
        }
        cnt++;
        n_hits++;
    }

    return {n_hits, partial_hits, matches_map};
}

std::ostream& operator<<(std::ostream& os, const Hit& hit) {
    os << "Hit(query_start=" << hit.query_start << ", query_end=" << hit.query_end << ", position=" << hit.position << ", is_partial=" << hit.is_partial << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const Match& match) {
    os << "Match(query_start=" << match.query_start << ", query_end=" << match.query_end << ", ref_start=" << match.ref_start << ", ref_end=" << match.ref_end << ")";
    return os;
}


std::ostream& operator<<(std::ostream& os, const Nam& n) {
    os << "Nam(ref_id=" << n.ref_id << ", query: " << n.query_start << ".." << n.query_end << ", ref: " << n.ref_start << ".." << n.ref_end << ", rc=" << static_cast<int>(n.is_revcomp) << ", score=" << n.score << ")";
    return os;
}
