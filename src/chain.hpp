#ifndef STROBEALIGN_CHAIN_HPP
#define STROBEALIGN_CHAIN_HPP

#include <vector>

#include "nam.hpp"
#include "index.hpp"
#include "mappingparameters.hpp"
#include "indexparameters.hpp"

struct Anchor {
    uint query_start;
    uint ref_start;
    uint ref_id;

    bool operator<(const Anchor& other) const {
        return (ref_id < other.ref_id) || (ref_id == other.ref_id && ref_start < other.ref_start) || (ref_id == other.ref_id && ref_start == other.ref_start && query_start < other.query_start);
    }

    bool operator==(const Anchor& other) const {
        return (ref_id == other.ref_id) && ref_start == other.ref_start && query_start == other.query_start;
    }
};

std::vector<Nam> get_chains(
    const std::string& seq,
    const StrobemerIndex& index,
    const MappingParameters& map_param,
    const IndexParameters& index_parameters
);

/**
 * Fast log2 function for x >= 1.
 * Copied from https://github.com/lh3/minimap2/blob/master/mmpriv.h
 */
static inline float mg_log2(float x) {
    union { float f; uint32_t i; } z = { x };
    float log_2 = ((z.i >> 23) & 255) - 128;
    z.i &= ~(255 << 23);
    z.i += 127 << 23;
    log_2 += (-0.34484843f * z.f + 2.02466578f) * z.f - 0.67487759f;
    return log_2;
}

#endif
