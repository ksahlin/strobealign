#include "doctest.h"
#include "randstrobes.hpp"
#include "revcomp.hpp"
#include "refs.hpp"


std::vector<Syncmer> syncmers_of(std::string& seq, SyncmerParameters parameters) {
    SyncmerIterator iterator{seq, parameters};
    std::vector<Syncmer> syncmers;
    Syncmer syncmer;
    while (!(syncmer = iterator.next()).is_end()) {
        syncmers.push_back(syncmer);
    }

    return syncmers;
}

TEST_CASE("RefRandstrobe constructor") {
    randstrobe_hash_t hash = 0x1234567890ABCDEF & RANDSTROBE_HASH_MASK;
    uint32_t position = ~0u;
    uint32_t ref_index = RefRandstrobe::max_number_of_references - 1;
    SUBCASE("one") {
        uint8_t offset = 255;
        RefRandstrobe rr{hash, position, ref_index, offset};

        CHECK(rr.hash() == hash);
        CHECK(rr.position() == position);
        CHECK(rr.reference_index() == ref_index);
        CHECK(rr.strobe2_offset() == offset);
    }

    SUBCASE("two") {
        uint8_t offset = 0;
        RefRandstrobe rr{hash, position, ref_index, offset};

        CHECK(rr.hash() == hash);
        CHECK(rr.position() == position);
        CHECK(rr.reference_index() == ref_index);
        CHECK(rr.strobe2_offset() == offset);
    }


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
