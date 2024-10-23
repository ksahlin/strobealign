#include "doctest.h"
#include "randstrobes.hpp"


TEST_CASE("test RefRandstrobe constructor") {
    randstrobe_hash_t hash = 0x1234567890ABCDEF & RANDSTROBE_HASH_MASK;
    uint32_t position = ~0u;
    uint32_t ref_index = (1u << 23) - 1;
    uint8_t offset = 255;
    bool first_strobe_is_main = true;
    RefRandstrobe rr{hash, position, ref_index, offset, first_strobe_is_main};

    CHECK(rr.hash() == hash);
    CHECK(rr.position() == position);
    CHECK(rr.reference_index() == ref_index);
    CHECK(rr.strobe2_offset() == offset);
    CHECK(rr.first_strobe_is_main() == first_strobe_is_main);
}
