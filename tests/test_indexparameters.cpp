#include "doctest.h"
#include "indexparameters.hpp"


TEST_CASE("test IndexParameters constructor") {
    size_t canonical_read_length = 250;

    int k = 22;
    int s = 18;

    int l = 2;
    int u = 12;
    unsigned w_min = 6;
    unsigned w_max = 16;
    int max_dist = 180;
    uint64_t q = 255;
    uint64_t main_hash_mask = 0xfffffffffc000000;
    int aux_len = 17;

    SyncmerParameters sp{k, s};
    RandstrobeParameters rp{q, max_dist, w_min, w_max, main_hash_mask};
    IndexParameters ip = IndexParameters{
        canonical_read_length,
        k,
        s,
        l,
        u,
        q,
        max_dist,
        aux_len
    };

    CHECK(ip.canonical_read_length == canonical_read_length);
    CHECK(ip.randstrobe == rp);
    CHECK(ip.syncmer == sp);

    auto def = IndexParameters::DEFAULT;
    IndexParameters ip2 = IndexParameters::from_read_length(canonical_read_length + 1, def, def, def, def, def, def);

    CHECK(ip2.canonical_read_length == canonical_read_length);
    CHECK(ip2.randstrobe == rp);
    CHECK(ip2.syncmer == sp);
}
