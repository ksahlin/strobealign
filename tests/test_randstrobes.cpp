#include "doctest.h"
#include "randstrobes.hpp"


TEST_CASE("RefRandstrobe constructor") {
    randstrobe_hash_t hash = 0x1234567890ABCDEF & RANDSTROBE_HASH_MASK;
    uint32_t position = ~0u;
    uint32_t ref_index = RefRandstrobe::max_number_of_references - 1;
    uint8_t offset = 255;
    bool first_strobe_is_main = true;
    RefRandstrobe rr{hash, position, ref_index, offset, first_strobe_is_main};

    CHECK(rr.hash() == hash);
    CHECK(rr.position() == position);
    CHECK(rr.reference_index() == ref_index);
    CHECK(rr.strobe2_offset() == offset);
    CHECK(rr.first_strobe_is_main() == first_strobe_is_main);
}

TEST_CASE("SyncmerIterator") {
    std::string seq{"AAAAAAAAAAAAAAAAAAAA"};
    CHECK(seq.size() == 20);

    SyncmerParameters parameters{8, 4};
    CHECK(parameters.t_syncmer == 3);
    SyncmerIterator si{seq, parameters};

    Syncmer syncmer;
    syncmer = si.next();
    CHECK(!syncmer.is_end());
    CHECK(syncmer.position == 0ul);
}
