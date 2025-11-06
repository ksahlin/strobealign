#ifndef STROBEALIGN_HITS_HPP
#define STROBEALIGN_HITS_HPP

#include <iostream>
#include <tuple>
#include <vector>

#include "randstrobes.hpp"
#include "index.hpp"
#include "mcsstrategy.hpp"

// A Hit is the result of successfully looking up a strobemer in the index
struct Hit {
    size_t position;
    size_t query_start;
    size_t query_end;
    bool is_partial;
};

// Aggregate statistics resulting from looking up all strobemers of a single
// query.
struct HitsDetails {
    uint full_not_found{0};
    uint full_filtered{0};  // found but filtered
    uint full_found{0};  // found and not filtered

    uint partial_not_found{0};
    uint partial_filtered{0};
    uint partial_found{0};

    uint rescued{0};  // first filtered but then rescued

    uint filtered_nucleotides{0};

    uint total_hits() const {
        return partial_filtered + partial_found + full_filtered + full_found;
    }

    // Used as a heuristic to compare the two orientations of a query
    bool is_better_than(HitsDetails& other) const {
        uint total = full_found + full_filtered;
        uint other_total = other.full_found + other.full_filtered;
        return total >= other_total * 2 && total > other_total + 5;
    }

    uint total_found() const {
        return full_filtered + full_found + partial_filtered + partial_found;
    }

    HitsDetails& operator+=(const HitsDetails& other) {
        full_not_found += other.full_not_found;
        full_filtered += other.full_filtered;
        full_found += other.full_found;
        partial_not_found += other.partial_not_found;
        partial_filtered += other.partial_filtered;
        partial_found += other.partial_found;
        rescued += other.rescued;

        return *this;
    }
};

std::ostream& operator<<(std::ostream& os, const Hit& hit);

std::tuple<HitsDetails, bool, std::vector<Hit>> find_hits(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    McsStrategy mcs_strategy,
    int rescue_threshold
);

std::ostream& operator<<(std::ostream& os, const HitsDetails& details);

#endif
