#ifndef STROBEALIGN_NAM_HPP
#define STROBEALIGN_NAM_HPP

#include <vector>
#include "index.hpp"
#include "randstrobes.hpp"

// Non-overlapping approximate match
struct nam {
    int nam_id;
    int query_s;
    int query_e;
    int query_prev_hit_startpos;
    int ref_s;
    int ref_e;
    int ref_prev_hit_startpos;
    int n_hits = 0;
    int ref_id;
    float score;
//    unsigned int previous_query_start;
//    unsigned int previous_ref_start;
    bool is_rc = false;

    int ref_span() const {
        return ref_e - ref_s;
    }

    int query_span() const {
        return query_e - query_s;
    }
};

std::pair<float,int> find_nams(
    std::vector<nam> &final_nams,
    const QueryRandstrobeVector &query_randstrobes,
    const StrobemerIndex& index,
    int k
);

void find_nams_rescue(
    std::vector<nam> &final_nams,
    const QueryRandstrobeVector &query_randstrobes,
    const StrobemerIndex& index,
    int k,
    unsigned int filter_cutoff
);

#endif
