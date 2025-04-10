#ifndef STROBEALIGN_NAM_HPP
#define STROBEALIGN_NAM_HPP

#include <vector>
#include "index.hpp"
#include "randstrobes.hpp"

struct Hit {
    size_t position;
    size_t query_start;
    size_t query_end;
    bool is_partial;
};

std::ostream& operator<<(std::ostream& os, const Hit& hit);

struct Match {
    int query_start;
    int query_end;
    int ref_start;
    int ref_end;
};

bool operator==(const Match& lhs, const Match& rhs);
std::ostream& operator<<(std::ostream& os, const Match& match);


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

std::ostream& operator<<(std::ostream& os, const Nam& nam);

std::tuple<int, int, bool, std::vector<Hit>> find_hits(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    bool use_mcs
);

std::tuple<int, int, robin_hood::unordered_map<unsigned int, std::vector<Match>>> find_matches_rescue(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    unsigned int rescue_cutoff,
    bool use_mcs
);

void merge_matches_into_nams(
    robin_hood::unordered_map<unsigned int, std::vector<Match>>& matches_map,
    int k,
    bool sort,
    bool is_revcomp,
    std::vector<Nam>& nams  // inout
);

robin_hood::unordered_map<unsigned int, std::vector<Match>> hits_to_matches(
    const std::vector<Hit>& hits,
    const StrobemerIndex& index
);

#endif
