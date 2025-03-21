#ifndef STROBEALIGN_NAM_HPP
#define STROBEALIGN_NAM_HPP

#include <vector>
#include <array>
#include "index.hpp"
#include "randstrobes.hpp"

// Non-overlapping approximate match
struct Nam {
    int nam_id;
    int query_start;
    int query_end;
    int query_prev_match_startpos;
    int ref_start;
    int ref_end;
    int ref_prev_match_startpos;
    int n_matches = 0;
    int ref_id;
    float score;
//    unsigned int previous_query_start;
//    unsigned int previous_ref_start;
    bool is_revcomp = false;

    int ref_span() const {
        return ref_end - ref_start;
    }

    int query_span() const {
        return query_end - query_start;
    }

    int projected_ref_start() const {
        return std::max(0, ref_start - query_start);
    }
};

std::tuple<float, int, int, std::vector<Nam>> find_nams(
    const std::vector<QueryRandstrobe> &query_randstrobes,
    const StrobemerIndex& index,
    bool use_mcs
);

std::tuple<int, int, std::vector<Nam>> find_nams_rescue(
    const std::vector<QueryRandstrobe> &query_randstrobes,
    const StrobemerIndex& index,
    unsigned int rescue_cutoff,
    bool use_mcs
);

std::ostream& operator<<(std::ostream& os, const Nam& nam);

#endif
