#include "nam.hpp"

namespace {

struct Hit {
    int query_start;
    int query_end;
    int ref_start;
    int ref_end;
};

struct PartialSeed {
    size_t hash;
    unsigned int start;
    bool is_reverse;
    bool operator==(const PartialSeed& rhs) const {
        return (hash == rhs.hash) && (start == rhs.start) && (is_reverse == rhs.is_reverse);
    }
};

inline void add_to_hits_per_ref_full(
    robin_hood::unordered_map<unsigned int, std::vector<Hit>>& hits_per_ref,
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
            hits_per_ref[index.reference_index(position)].push_back(
                Hit{query_start, query_end, ref_start, ref_end}
            );
            min_diff = diff;
        }
    }
}

inline void add_to_hits_per_ref_partial(
    robin_hood::unordered_map<unsigned int, std::vector<Hit>>& hits_per_ref,
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
        hits_per_ref[index.reference_index(position)].push_back(
            Hit{query_start, query_end, ref_start, ref_end}
        );
    }
}

bool operator==(const Hit& lhs, const Hit& rhs) {
    return (lhs.query_start == rhs.query_start) && (lhs.query_end == rhs.query_end) && (lhs.ref_start == rhs.ref_start) && (lhs.ref_end == rhs.ref_end);
}

void merge_hits_into_nams(
    robin_hood::unordered_map<unsigned int, std::vector<Hit>>& hits_per_ref,
    int k,
    bool sort,
    bool is_revcomp,
    std::vector<Nam>& nams  // inout
) {
    for (auto &[ref_id, hits] : hits_per_ref) {
        if (sort) {
            std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) -> bool {
                    // first sort on query starts, then on reference starts, finally prefer full matches over partial
                    return (a.query_start < b.query_start) ||
                            ( (a.query_start == b.query_start) && (a.ref_start < b.ref_start)) ||
                            ( (a.query_start == b.query_start) && (a.ref_start == b.ref_start) && (a.query_end > b.query_end)  );
                }
            );
        }

        std::vector<Nam> open_nams;
        int prev_q_start = 0;
        auto prev_hit = Hit{0,0,0,0};
        for (auto &h : hits) {
            if ( (prev_hit == h) || ( ((h.query_end - h.query_start) == k) && (prev_hit.query_start == h.query_start) && (prev_hit.ref_start == h.ref_start)) )  { // is a redundant partial hit
                continue;
            }
            bool is_added = false;
            for (auto & o : open_nams) {

                // Extend NAM
                if ((o.query_prev_hit_startpos <= h.query_start) && (h.query_start <= o.query_end ) && (o.ref_prev_hit_startpos <= h.ref_start) && (h.ref_start <= o.ref_end) ){
                    if ( (h.query_end > o.query_end) && (h.ref_end > o.ref_end) ) {
                        o.query_end = h.query_end;
                        o.ref_end = h.ref_end;
//                        o.previous_query_start = h.query_s;
//                        o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
                        o.query_prev_hit_startpos = h.query_start;
                        o.ref_prev_hit_startpos = h.ref_start;
                        o.n_hits ++;
//                        o.score += (float)1/ (float)h.count;
                        is_added = true;
                        break;
                    }
                    else if ((h.query_end <= o.query_end) && (h.ref_end <= o.ref_end)) {
//                        o.previous_query_start = h.query_s;
//                        o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
                        o.query_prev_hit_startpos = h.query_start;
                        o.ref_prev_hit_startpos = h.ref_start;
                        o.n_hits ++;
//                        o.score += (float)1/ (float)h.count;
                        is_added = true;
                        break;
                    }
                }

            }
            // Add the hit to open matches
            if (!is_added){
                Nam n;
                n.query_start = h.query_start;
                n.query_end = h.query_end;
                n.ref_start = h.ref_start;
                n.ref_end = h.ref_end;
                n.ref_id = ref_id;
//                n.previous_query_start = h.query_s;
//                n.previous_ref_start = h.ref_s;
                n.query_prev_hit_startpos = h.query_start;
                n.ref_prev_hit_startpos = h.ref_start;
                n.n_hits = 1;
                n.is_rc = is_revcomp;
//                n.score += (float)1 / (float)h.count;
                open_nams.push_back(n);
            }

            // Only filter if we have advanced at least k nucleotides
            if (h.query_start > prev_q_start + k) {

                // Output all NAMs from open_matches to final_nams that the current hit have passed
                for (auto &n : open_nams) {
                    if (n.query_end < h.query_start) {
                        int n_max_span = std::max(n.query_span(), n.ref_span());
                        int n_min_span = std::min(n.query_span(), n.ref_span());
                        float n_score;
                        n_score = ( 2*n_min_span -  n_max_span) > 0 ? (float) (n.n_hits * ( 2*n_min_span -  n_max_span) ) : 1;   // this is really just n_hits * ( min_span - (offset_in_span) ) );
//                        n_score = n.n_hits * n.query_span();
                        n.score = n_score;
                        n.nam_id = nams.size();
                        nams.push_back(n);
                    }
                }

                // Remove all NAMs from open_matches that the current hit have passed
                auto c = h.query_start;
                auto predicate = [c](decltype(open_nams)::value_type const &nam) { return nam.query_end < c; };
                open_nams.erase(std::remove_if(open_nams.begin(), open_nams.end(), predicate), open_nams.end());
                prev_q_start = h.query_start;
            }
            prev_hit = h;
        }

        // Add all current open_matches to final NAMs
        for (auto &n : open_nams) {
            int n_max_span = std::max(n.query_span(), n.ref_span());
            int n_min_span = std::min(n.query_span(), n.ref_span());
            float n_score;
            n_score = ( 2*n_min_span -  n_max_span) > 0 ? (float) (n.n_hits * ( 2*n_min_span -  n_max_span) ) : 1;   // this is really just n_hits * ( min_span - (offset_in_span) ) );
//            n_score = n.n_hits * n.query_span();
            n.score = n_score;
            n.nam_id = nams.size();
            nams.push_back(n);
        }
    }
}

std::vector<Nam> merge_hits_into_nams_forward_and_reverse(
    std::array<robin_hood::unordered_map<unsigned int, std::vector<Hit>>, 2>& hits_per_ref,
    int k,
    bool sort
) {
    std::vector<Nam> nams;
    for (size_t is_revcomp = 0; is_revcomp < 2; ++is_revcomp) {
        auto& hits_oriented = hits_per_ref[is_revcomp];
        merge_hits_into_nams(hits_oriented, k, sort, is_revcomp, nams);
    }
    return nams;
}

} // namespace

/*
 * Find a query’s NAMs, ignoring randstrobes that occur too often in the
 * reference (have a count above filter_cutoff).
 *
 * Return the fraction of nonrepetitive hits (those not above the filter_cutoff threshold)
 */
std::tuple<float, int, std::vector<Nam>> find_nams(
    const QueryRandstrobeVector &query_randstrobes,
    const StrobemerIndex& index,
    bool use_mcs
) {
    std::vector<PartialSeed> partial_queried; // TODO: is a small set more efficient than linear search in a small vector?
    partial_queried.reserve(10);
    std::array<robin_hood::unordered_map<unsigned int, std::vector<Hit>>, 2> hits_per_ref;
    hits_per_ref[0].reserve(100);
    hits_per_ref[1].reserve(100);
    int nr_good_hits = 0;
    int total_hits = 0;
    for (const auto &q : query_randstrobes) {
        size_t position = index.find(q.hash);
        if (position != index.end()){
            total_hits++;
            if (index.is_filtered(position)) {
                continue;
            }
            nr_good_hits++;
            add_to_hits_per_ref_full(hits_per_ref[q.is_reverse], q.start, q.end, index, position);
        }
        else if (use_mcs) {
            PartialSeed ph{q.hash >> index.get_aux_len(), q.partial_start, q.is_reverse};
            if (std::find(partial_queried.begin(), partial_queried.end(), ph) != partial_queried.end()) {
                // already queried
                continue;
            }
            size_t partial_pos = index.partial_find(q.hash);
            if (partial_pos != index.end()) {
                total_hits++;
                if (index.is_partial_filtered(partial_pos)) {
                    partial_queried.push_back(ph);
                    continue;
                }
                nr_good_hits++;
                add_to_hits_per_ref_partial(hits_per_ref[q.is_reverse], q.partial_start, q.partial_end, index, partial_pos);
            }
            partial_queried.push_back(ph);
        }
    }
    float nonrepetitive_fraction = total_hits > 0 ? ((float) nr_good_hits) / ((float) total_hits) : 1.0;
    auto nams = merge_hits_into_nams_forward_and_reverse(hits_per_ref, index.k(), true);
    return {nonrepetitive_fraction, nr_good_hits, nams};
}

/*
 * Find a query’s NAMs, using also some of the randstrobes that occur more often
 * than filter_cutoff.
 *
 * Return the number of hits and the vector of NAMs.
 */
std::pair<int, std::vector<Nam>> find_nams_rescue(
    const QueryRandstrobeVector &query_randstrobes,
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
    std::vector<PartialSeed> partial_queried;  // TODO: is a small set more efficient than linear search in a small vector?
    partial_queried.reserve(10);
    std::array<robin_hood::unordered_map<unsigned int, std::vector<Hit>>, 2> hits_per_ref;
    std::vector<RescueHit> hits_fw;
    std::vector<RescueHit> hits_rc;
    hits_per_ref[0].reserve(100);
    hits_per_ref[1].reserve(100);
    hits_fw.reserve(5000);
    hits_rc.reserve(5000);

    for (auto &qr : query_randstrobes) {
        size_t position = index.find(qr.hash);
        if (position != index.end()) {
            unsigned int count = index.get_count(position);
            RescueHit rh{position, count, qr.start, qr.end, false};
            if (qr.is_reverse){
                hits_rc.push_back(rh);
            } else {
                hits_fw.push_back(rh);
            }
        }
        else if (use_mcs) {
            PartialSeed ph = {qr.hash >> index.get_aux_len(), qr.partial_start, qr.is_reverse};
            if (std::find(partial_queried.begin(), partial_queried.end(), ph) != partial_queried.end()) {
                // already queried
                continue;
            }
            size_t partial_pos = index.partial_find(qr.hash);
            if (partial_pos != index.end()) {
                unsigned int partial_count = index.get_partial_count(partial_pos);
                RescueHit rh{partial_pos, partial_count, qr.partial_start, qr.partial_end, true};
                if (qr.is_reverse){
                    hits_rc.push_back(rh);
                } else {
                    hits_fw.push_back(rh);
                }
            }
            partial_queried.push_back(ph);
        }
    }

    std::sort(hits_fw.begin(), hits_fw.end());
    std::sort(hits_rc.begin(), hits_rc.end());
    int n_hits = 0;
    size_t is_revcomp = 0;
    for (auto& rescue_hits : {hits_fw, hits_rc}) {
        int cnt = 0;
        for (auto &rh : rescue_hits) {
            if ((rh.count > rescue_cutoff && cnt >= 5) || rh.count > 1000) {
                break;
            }
            if (rh.is_partial){
                add_to_hits_per_ref_partial(hits_per_ref[is_revcomp], rh.query_start, rh.query_end, index, rh.position);
            } else{
                add_to_hits_per_ref_full(hits_per_ref[is_revcomp], rh.query_start, rh.query_end, index, rh.position);
            }
            cnt++;
            n_hits++;
        }
        is_revcomp++;
    }

    return {n_hits, merge_hits_into_nams_forward_and_reverse(hits_per_ref, index.k(), true)};
}

std::ostream& operator<<(std::ostream& os, const Nam& n) {
    os << "Nam(ref_id=" << n.ref_id << ", query: " << n.query_start << ".." << n.query_end << ", ref: " << n.ref_start << ".." << n.ref_end << ", score=" << n.score << ")";
    return os;
}
