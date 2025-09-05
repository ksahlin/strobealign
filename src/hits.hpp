#ifndef STROBEALIGN_HITS_HPP
#define STROBEALIGN_HITS_HPP

#include <iostream>
#include <tuple>
#include <vector>

#include "randstrobes.hpp"
#include "index.hpp"
#include "mcsstrategy.hpp"

/* A Hit is the result of successfully looking up a strobemer in the index
 *
 */
struct Hit {
    size_t position;
    size_t query_start;
    size_t query_end;
    bool is_partial;
};

std::ostream& operator<<(std::ostream& os, const Hit& hit);

std::tuple<int, int, bool, std::vector<Hit>> find_hits(
    const std::vector<QueryRandstrobe>& query_randstrobes,
    const StrobemerIndex& index,
    McsStrategy mcs_strategy
);

#endif
