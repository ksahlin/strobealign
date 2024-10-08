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
    bool is_rc = false;

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

struct NamFinder {
    NamFinder(const StrobemerIndex& index, bool use_mcs)
    : index(index)
    , use_mcs(use_mcs) {
    }

    std::tuple<float, int, std::vector<Nam>> find(const QueryRandstrobeVector &query_randstrobes) const;
    std::pair<int, std::vector<Nam>> find_rescue(
        const QueryRandstrobeVector &query_randstrobes,
        unsigned int rescue_cutoff
    ) const;

private:
    const StrobemerIndex& index;
    bool use_mcs;

    mutable std::vector<PartialHit> partial_queried;
};

std::ostream& operator<<(std::ostream& os, const Nam& nam);

#endif
