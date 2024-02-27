#include "nam.hpp"

namespace {

struct Hit {
    int query_start;
    int query_end;
    int ref_start;
    int ref_end;
};

inline void add_to_hits_per_ref(
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
            hits_per_ref[index.reference_index(position)].push_back(Hit{query_start, query_end, ref_start, ref_end});
            min_diff = diff;
        }
    }
}

void merge_hits_into_nams_one_ref(
    std::vector<Hit>& hits,
    size_t ref_id,
    int k,
    bool is_revcomp,
    std::vector<Nam>& nams  // inout
) {
    std::vector<Nam> open_nams;
    unsigned int prev_q_start = 0;
    for (auto &h : hits) {
        bool is_added = false;
        for (auto & o : open_nams) {

            // Extend NAM
            if ((o.query_prev_hit_startpos < h.query_start) && (h.query_start <= o.query_end ) && (o.ref_prev_hit_startpos < h.ref_start) && (h.ref_start <= o.ref_end) ){
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

void merge_nearby_nams(std::vector<Nam>& nams) {
    if (nams.empty()) {
        return;
    }
    std::sort(nams.begin(), nams.end(), [](const Nam& a, const Nam& b) -> bool {
        // by diagonal
        // TODO make branchless
        return (a.diagonal() < b.diagonal()) || ((a.diagonal() == b.diagonal()) && (a.ref_start < b.ref_start));
    });
    if (nams.empty()) {
        return;
    }
    size_t j = 0;
    for (size_t i = 1; i < nams.size(); i++) {
        if (nams[j].ref_id == nams[i].ref_id && nams[j].diagonal() == nams[i].diagonal() && nams[j].is_rc == nams[i].is_rc) {
            assert(nams[j].is_rc == nams[i].is_rc);
            // merge // TODO turn into method
            nams[j].query_end = nams[i].query_end;
            nams[j].ref_end = nams[i].ref_end;
            nams[j].n_hits += nams[i].n_hits;
            int n_max_span = std::max(nams[j].query_span(), nams[j].ref_span());
            int n_min_span = std::min(nams[j].query_span(), nams[j].ref_span());
            float n_score;
            n_score = ( 2*n_min_span -  n_max_span) > 0 ? (float) (nams[j].n_hits * ( 2*n_min_span -  n_max_span) ) : 1;   // this is really just n_hits * ( min_span - (offset_in_span) ) );
//            n_score = n.n_hits * n.query_span();
            nams[j].score = n_score;
        } else {
            j++;
        }
    }
    nams.resize(j + 1);
}

std::vector<Nam> merge_hits_into_nams(
    std::array<robin_hood::unordered_map<unsigned int, std::vector<Hit>>, 2>& hits_per_ref,
    int k,
    bool sort
) {
    std::vector<Nam> nams;
    for (size_t is_revcomp = 0; is_revcomp < 2; ++is_revcomp) {
        auto& hits_oriented = hits_per_ref[is_revcomp];
        for (auto &[ref_id, hits] : hits_oriented) {
            if (sort) {
                std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) -> bool {
                        // first sort on query starts, then on reference starts
                        return (a.query_start < b.query_start) || ( (a.query_start == b.query_start) && (a.ref_start < b.ref_start) );
                    }
                );
            }
            merge_hits_into_nams_one_ref(hits, ref_id, k, is_revcomp, nams);
            merge_nearby_nams(nams);
        }
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
std::pair<float, std::vector<Nam>> find_nams(
    const QueryRandstrobeVector &query_randstrobes,
    const StrobemerIndex& index
) {
    std::array<robin_hood::unordered_map<unsigned int, std::vector<Hit>>, 2> hits_per_ref;
    hits_per_ref[0].reserve(100);
    hits_per_ref[1].reserve(100);
    int nr_good_hits = 0, total_hits = 0;
    for (const auto &q : query_randstrobes) {
        size_t position = index.find(q.hash);
        if (position != index.end()){
            total_hits++;
            if (index.is_filtered(position)) {
                continue;
            }
            nr_good_hits++;
            add_to_hits_per_ref(hits_per_ref[q.is_reverse], q.start, q.end, index, position);
        }
    }
    float nonrepetitive_fraction = total_hits > 0 ? ((float) nr_good_hits) / ((float) total_hits) : 1.0;
    auto nams = merge_hits_into_nams(hits_per_ref, index.k(), false);
    return make_pair(nonrepetitive_fraction, nams);
}

/*
 * Find a query’s NAMs, using also some of the randstrobes that occur more often
 * than filter_cutoff.
 *
 */
std::vector<Nam> find_nams_rescue(
    const QueryRandstrobeVector &query_randstrobes,
    const StrobemerIndex& index,
    unsigned int rescue_cutoff
) {
    struct RescueHit {
        size_t position;
        unsigned int count;
        unsigned int query_start;
        unsigned int query_end;

        bool operator< (const RescueHit& rhs) const {
            return std::tie(count, query_start, query_end)
                < std::tie(rhs.count, rhs.query_start, rhs.query_end);
        }
    };

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
            RescueHit rh{position, count, qr.start, qr.end};
            if (qr.is_reverse){
                hits_rc.push_back(rh);
            } else {
                hits_fw.push_back(rh);
            }
        }
    }

    std::sort(hits_fw.begin(), hits_fw.end());
    std::sort(hits_rc.begin(), hits_rc.end());
    size_t is_revcomp = 0;
    for (auto& rescue_hits : {hits_fw, hits_rc}) {
        int cnt = 0;
        for (auto &rh : rescue_hits) {
            if ((rh.count > rescue_cutoff && cnt >= 5) || rh.count > 1000) {
                break;
            }
            add_to_hits_per_ref(hits_per_ref[is_revcomp], rh.query_start, rh.query_end, index, rh.position);
            cnt++;
        }
        is_revcomp++;
    }

    return merge_hits_into_nams(hits_per_ref, index.k(), true);
}

std::ostream& operator<<(std::ostream& os, const Nam& n) {
    os << "Nam(ref_id=" << n.ref_id << ", query: " << n.query_start << ".." << n.query_end << ", ref: " << n.ref_start << ".." << n.ref_end << ", score=" << n.score << ")";
    return os;
}
