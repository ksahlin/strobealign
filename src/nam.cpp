#include "nam.hpp"

bool operator==(const Match& lhs, const Match& rhs) {
    return (lhs.query_start == rhs.query_start) && (lhs.query_end == rhs.query_end) && (lhs.ref_start == rhs.ref_start) && (lhs.ref_end == rhs.ref_end);
}

namespace {

/*
 * A partial hit is a hit where not the full randstrobe hash could be found in
 * the index but only the "main" hash (only the first aux_len bits).
 */
struct PartialHit {
    randstrobe_hash_t hash;
    unsigned int start;  // position in strobemer index
    bool is_reverse;
    bool operator==(const PartialHit& rhs) const {
        return (hash == rhs.hash) && (start == rhs.start) && (is_reverse == rhs.is_reverse);
    }
};

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
        int ref_end = ref_start + index.strobe2_offset(position) + index.k();
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

} // namespace

std::vector<Nam> merge_matches_into_nams_forward_and_reverse(
    std::array<robin_hood::unordered_map<unsigned int, std::vector<Match>>, 2>& matches_map,
    int k,
    bool sort
) {
    std::vector<Nam> nams;
    for (size_t is_revcomp = 0; is_revcomp < 2; ++is_revcomp) {
        auto& matches_oriented = matches_map[is_revcomp];
        merge_matches_into_nams(matches_oriented, k, sort, is_revcomp, nams);
    }
    return nams;
}


/*
 * Find a query’s NAMs, ignoring randstrobes that occur too often in the
 * reference (have a count above filter_cutoff).
 *
 * Return the fraction of nonrepetitive hits (those not above the filter_cutoff threshold)
 */
std::tuple<float, int, int, bool, std::array<robin_hood::unordered_map<unsigned int, std::vector<Match>>, 2>> find_matches(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    bool use_mcs
) {
    // If we produce matches in sorted order, then merge_matches_into_nams()
    // does not have to re-sort
    bool sorting_needed{use_mcs};
    std::vector<PartialHit> partial_queried; // TODO: is a small set more efficient than linear search in a small vector?
    if (use_mcs) {
        partial_queried.reserve(10);
    }
    std::array<robin_hood::unordered_map<unsigned int, std::vector<Match>>, 2> matches_map;
    matches_map[0].reserve(100);
    matches_map[1].reserve(100);
    int nonrepetitive_hits = 0;
    int total_hits = 0;
    int partial_hits = 0;
    for (const auto &q : query_randstrobes) {
        size_t position = index.find_full(q.hash);
        if (position != index.end()) {
            total_hits++;
            if (index.is_filtered(position)) {
                continue;
            }
            nonrepetitive_hits++;
            add_to_matches_map_full(matches_map[q.is_revcomp], q.start, q.end, index, position);
        }
        else if (use_mcs) {
            PartialHit ph{q.hash & index.get_main_hash_mask(), q.partial_start, q.is_revcomp};
            if (std::find(partial_queried.begin(), partial_queried.end(), ph) != partial_queried.end()) {
                // already queried
                continue;
            }
            size_t partial_pos = index.find_partial(q.hash);
            if (partial_pos != index.end()) {
                total_hits++;
                if (index.is_partial_filtered(partial_pos)) {
                    partial_queried.push_back(ph);
                    continue;
                }
                nonrepetitive_hits++;
                partial_hits++;
                add_to_matches_map_partial(matches_map[q.is_revcomp], q.partial_start, q.partial_end, index, partial_pos);
            }
            partial_queried.push_back(ph);
        }
    }

    // Rescue using partial hits, even in non-MCS mode
    if (total_hits == 0 && !use_mcs) {
        for (const auto &q : query_randstrobes) {
            size_t partial_pos = index.find_partial(q.hash);
            if (partial_pos != index.end()) {
                total_hits++;
                if (index.is_partial_filtered(partial_pos)) {
                    continue;
                }
                nonrepetitive_hits++;
                partial_hits++;
                add_to_matches_map_partial(matches_map[q.is_revcomp], q.partial_start, q.partial_end, index, partial_pos);
            }
        }
        sorting_needed = true;
    }

    float nonrepetitive_fraction = total_hits > 0 ? ((float) nonrepetitive_hits) / ((float) total_hits) : 1.0;
    return {nonrepetitive_fraction, nonrepetitive_hits, partial_hits, sorting_needed, matches_map};
}

/*
 * Find a query’s NAMs, using also some of the randstrobes that occur more often
 * than filter_cutoff.
 *
 * Return the number of hits and the vector of NAMs.
 */
std::tuple<int, int, std::vector<Nam>> find_nams_rescue(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    unsigned int rescue_cutoff,
    bool use_mcs
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
    std::vector<PartialHit> partial_queried;  // TODO: is a small set more efficient than linear search in a small vector?
    partial_queried.reserve(10);
    std::array<robin_hood::unordered_map<unsigned int, std::vector<Match>>, 2> matches_map;
    std::array<std::vector<RescueHit>, 2> hits;
    matches_map[0].reserve(100);
    matches_map[1].reserve(100);
    hits[0].reserve(5000);
    hits[1].reserve(5000);

    for (auto &qr : query_randstrobes) {
        size_t position = index.find_full(qr.hash);
        if (position != index.end()) {
            unsigned int count = index.get_count_full(position);
            RescueHit rh{position, count, qr.start, qr.end, false};
            hits[qr.is_revcomp].push_back(rh);
        }
        else if (use_mcs) {
            PartialHit ph = {qr.hash & index.get_main_hash_mask(), qr.partial_start, qr.is_revcomp};
            if (std::find(partial_queried.begin(), partial_queried.end(), ph) != partial_queried.end()) {
                // already queried
                continue;
            }
            size_t partial_pos = index.find_partial(qr.hash);
            if (partial_pos != index.end()) {
                unsigned int partial_count = index.get_count_partial(partial_pos);
                RescueHit rh{partial_pos, partial_count, qr.partial_start, qr.partial_end, true};
                hits[qr.is_revcomp].push_back(rh);
            }
            partial_queried.push_back(ph);
        }
    }

    std::sort(hits[0].begin(), hits[0].end());
    std::sort(hits[1].begin(), hits[1].end());
    int n_hits = 0;
    int partial_hits = 0;
    size_t is_revcomp = 0;
    for (auto& rescue_hits : hits) {
        int cnt = 0;
        for (auto &rh : rescue_hits) {
            if ((rh.count > rescue_cutoff && cnt >= 5) || rh.count > 1000) {
                break;
            }
            if (rh.is_partial){
                partial_hits++;
                add_to_matches_map_partial(matches_map[is_revcomp], rh.query_start, rh.query_end, index, rh.position);
            } else{
                add_to_matches_map_full(matches_map[is_revcomp], rh.query_start, rh.query_end, index, rh.position);
            }
            cnt++;
            n_hits++;
        }
        is_revcomp++;
    }

    return {n_hits, partial_hits, merge_matches_into_nams_forward_and_reverse(matches_map, index.k(), true)};
}

std::ostream& operator<<(std::ostream& os, const Nam& n) {
    os << "Nam(ref_id=" << n.ref_id << ", query: " << n.query_start << ".." << n.query_end << ", ref: " << n.ref_start << ".." << n.ref_end << ", rc=" << static_cast<int>(n.is_revcomp) << ", score=" << n.score << ")";
    return os;
}
