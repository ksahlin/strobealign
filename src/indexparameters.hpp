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
    const int l;
    const int u;
    const uint64_t q;
    const int max_dist;
    const unsigned w_min;
    const unsigned w_max;

    RandstrobeParameters(int l, int u, uint64_t q, int max_dist, unsigned w_min, unsigned w_max)
        : l(l)
        , u(u)
        , q(q)
        , max_dist(max_dist)
        , w_min(w_min)
        , w_max(w_max)
    {
        verify();
    }

    bool operator==(const RandstrobeParameters& other) const;

private:
    void verify() const {
        if (max_dist > 255) {
            throw BadParameter("maximum seed length (-m <max_dist>) is larger than 255");
        }
        if (w_min > w_max) {
            throw BadParameter("w_min is greater than w_max (choose different -l/-u parameters)");
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

    IndexParameters(size_t canonical_read_length, int k, int s, int l, int u, int q, int max_dist)
        : canonical_read_length(canonical_read_length)
        , syncmer(k, s)
        , randstrobe(l, u, q, max_dist, std::max(0, k / (k - s + 1) + l), k / (k - s + 1) + u)
    {
    }

    static IndexParameters from_read_length(
        int read_length, int k = DEFAULT, int s = DEFAULT, int l = DEFAULT, int u = DEFAULT, int c = DEFAULT, int max_seed_len = DEFAULT
    );
    static IndexParameters read(std::istream& os);
    std::string filename_extension() const;
    void write(std::ostream& os) const;
    bool operator==(const IndexParameters& other) const;
    bool operator!=(const IndexParameters& other) const { return !(*this == other); }
};

std::ostream& operator<<(std::ostream& os, const IndexParameters& parameters);

#endif
