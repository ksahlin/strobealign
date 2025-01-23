#ifndef STROBEALIGN_INDEXPARAMETERS_HPP
#define STROBEALIGN_INDEXPARAMETERS_HPP

#include <cstdint>
#include <algorithm>
#include <iostream>
#include <limits>
#include "exceptions.hpp"


struct SyncmerParameters {
    const int k;
    const int s;
    const int t_syncmer;

    SyncmerParameters(int k, int s)
        : k(k)
        , s(s)
        , t_syncmer((k - s) / 2 + 1)
    {
        verify();
    }

    void verify() const {
        if (k <= 7 || k > 32) {
            throw BadParameter("k not in [8,32]");
        }
        if (s > k) {
            throw BadParameter("s is larger than k");
        }
        if ((k - s) % 2 != 0) {
            throw BadParameter("(k - s) must be an even number to create canonical syncmers. Please set s to e.g. k-2, k-4, k-6, ...");
        }
    }

    bool operator==(const SyncmerParameters& other) const;
};

struct RandstrobeParameters {
    const uint64_t q;
    const int max_dist;
    const unsigned w_min;
    const unsigned w_max;
    const uint64_t main_hash_mask;
    const uint64_t strobe2_mask;

    RandstrobeParameters(
        uint64_t q,
        int max_dist,
        unsigned w_min,
        unsigned w_max,
        uint64_t main_hash_mask,
        uint64_t strobe2_mask
    )
        : q(q)
        , max_dist(max_dist)
        , w_min(w_min)
        , w_max(w_max)
        , main_hash_mask(main_hash_mask)
        , strobe2_mask(strobe2_mask)
    {
    }

    bool operator==(const RandstrobeParameters& other) const;

private:
    void verify(unsigned aux_len) const {
        if (max_dist > 255) {
            throw BadParameter("maximum seed length (-m <max_dist>) is larger than 255");
        }
        if (w_min > w_max) {
            throw BadParameter("l is greater than u");
        }
        if (aux_len > 63) {
            throw BadParameter("aux_len is larger than 63");
        }
    }
};

/* Settings that influence index creation */
class IndexParameters {
public:
    const size_t canonical_read_length;
    const SyncmerParameters syncmer;
    const RandstrobeParameters randstrobe;

    static const int DEFAULT = std::numeric_limits<int>::min();

    IndexParameters(size_t canonical_read_length, int k, int s, int w_min, int w_max, uint64_t q, int max_dist, int aux_len)
        : canonical_read_length(canonical_read_length)
        , syncmer(k, s)
        , randstrobe(q, max_dist, w_min, w_max, ~0ul << (16 + 2 * aux_len))
    {
        if (w_min < 0 || w_max < 0) {
            throw BadParameter("Neither l nor u can be negative");
        }
        verify(aux_len);
    }

    IndexParameters(size_t canonical_read_length, SyncmerParameters syncmer, RandstrobeParameters randstrobe)
        : canonical_read_length(canonical_read_length)
        , syncmer(syncmer)
        , randstrobe(randstrobe)
    {
    }

    void verify(unsigned aux_len) const {
        if (aux_len < 6 || aux_len > 27) {
            throw BadParameter("aux_len must be between 6 and 27");
        }
    }

    static IndexParameters from_read_length(
        int read_length, int k = DEFAULT, int s = DEFAULT, int w_min = DEFAULT, int w_max = DEFAULT, int c = DEFAULT, int max_seed_len = DEFAULT, int aux_len = DEFAULT);
    static IndexParameters read(std::istream& os);
    std::string filename_extension() const;
    void write(std::ostream& os) const;
    bool operator==(const IndexParameters& other) const;
    bool operator!=(const IndexParameters& other) const { return !(*this == other); }
};

std::ostream& operator<<(std::ostream& os, const SyncmerParameters& parameters);
std::ostream& operator<<(std::ostream& os, const RandstrobeParameters& parameters);
std::ostream& operator<<(std::ostream& os, const IndexParameters& parameters);

#endif
