#ifndef INDEXPARAMETERS_HPP
#define INDEXPARAMETERS_HPP

#include <cstdint>
#include <algorithm>
#include <iostream>
#include <limits>
#include "exceptions.hpp"

/* Settings that influence index creation */
class IndexParameters {
public:
    const size_t canonical_read_length;
    const int k;
    const int s;
    const int l;
    const int u;
    const uint64_t q;
    const int max_dist;
    const int t_syncmer;
    const unsigned w_min;
    const unsigned w_max;

    static const int DEFAULT = std::numeric_limits<int>::min();
    IndexParameters(size_t canonical_read_length, int k, int s, int l, int u, int q, int max_dist)
        : canonical_read_length(canonical_read_length)
        , k(k)
        , s(s)
        , l(l)
        , u(u)
        , q(q)
        , max_dist(max_dist)
        , t_syncmer((k - s) / 2 + 1)
        , w_min(std::max(1, k / (k - s + 1) + l))
        , w_max(k / (k - s + 1) + u)
    {
        verify();
    }

    static IndexParameters from_read_length(int read_length, int c = DEFAULT, int k = DEFAULT, int s = DEFAULT, int l = DEFAULT, int u = DEFAULT, int max_seed_len = DEFAULT);
    static IndexParameters read(std::istream& os);
    std::string filename_extension() const;
    void write(std::ostream& os) const;
    bool operator==(const IndexParameters& other) const;
    bool operator!=(const IndexParameters& other) const { return !(*this == other); }

private:
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
        if (max_dist > 255) {
            throw BadParameter("maximum seed length (-m <max_dist>) is larger than 255");
        }
    }
};

std::ostream& operator<<(std::ostream& os, const IndexParameters& parameters);

#endif
