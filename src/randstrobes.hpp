#ifndef STROBEALIGN_RANDSTROBES_HPP
#define STROBEALIGN_RANDSTROBES_HPP

#include <vector>
#include <string>
#include <tuple>
#include <deque>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <inttypes.h>

#include "indexparameters.hpp"


using syncmer_hash_t = uint64_t;
using randstrobe_hash_t = uint64_t;

struct RefRandstrobe {
    using packed_t = uint32_t;
    randstrobe_hash_t hash;
    uint32_t position;

    RefRandstrobe() { }

    RefRandstrobe(randstrobe_hash_t hash, uint32_t position, uint32_t packed)
        : hash(hash)
        , position(position)
        , m_packed(packed) { }

    bool operator< (const RefRandstrobe& other) const {
        return hash < other.hash;
    }

    int reference_index() const {
        return m_packed >> bit_alloc;
    }

    int strobe2_offset() const {
        return m_packed & mask;
    }

private:
    static constexpr int bit_alloc = 8;
    static constexpr int mask = (1 << bit_alloc) - 1;
    packed_t m_packed; // packed representation of ref_index and strobe offset
};

struct QueryRandstrobe {
    randstrobe_hash_t hash;
    unsigned int start;
    unsigned int end;
    bool is_reverse;
};

std::ostream& operator<<(std::ostream& os, const QueryRandstrobe& randstrobe);

using QueryRandstrobeVector = std::vector<QueryRandstrobe>;

QueryRandstrobeVector randstrobes_query(const std::string_view seq, const IndexParameters& parameters);

struct Randstrobe {
    randstrobe_hash_t hash;
    unsigned int strobe1_pos;
    unsigned int strobe2_pos;

    bool operator==(const Randstrobe& other) const {
        return hash == other.hash && strobe1_pos == other.strobe1_pos && strobe2_pos == other.strobe2_pos;
    }

    bool operator!=(const Randstrobe& other) const {
        return !(*this == other);
    }
};

std::ostream& operator<<(std::ostream& os, const Randstrobe& randstrobe);

class RandstrobeIterator {
public:
    RandstrobeIterator(
        const std::vector<uint64_t>& syncmer_hashes,
        const std::vector<unsigned int>& index_to_coordinate,
        unsigned w_min,
        unsigned w_max,
        uint64_t q,
        int max_dist
    ) : syncmer_hashes(syncmer_hashes)
      , index_to_coordinate(index_to_coordinate)
      , w_min(w_min)
      , w_max(w_max)
      , q(q)
      , max_dist(max_dist)
    {
        if (w_min > w_max) {
            throw std::invalid_argument("w_min is greater than w_max");
        }
    }

    Randstrobe next() {
        return get(strobe1_index++);
    }

    bool has_next() {
        return strobe1_index + w_min < syncmer_hashes.size();
    }

private:
    Randstrobe get(unsigned int strobe1_index) const;
    const std::vector<uint64_t>& syncmer_hashes;
    const std::vector<unsigned int>& index_to_coordinate;
    const unsigned w_min;
    const unsigned w_max;
    const uint64_t q;
    const unsigned int max_dist;
    unsigned int strobe1_index = 0;
};

struct Syncmer {
    syncmer_hash_t hash;
    size_t position;
    bool is_end() const {
        return hash == 0 && position == 0;
    }
};

std::ostream& operator<<(std::ostream& os, const Syncmer& syncmer);

class SyncmerIterator {
public:
    SyncmerIterator(const std::string_view seq, size_t k, size_t s, size_t t)
        : seq(seq), k(k), s(s), t(t) { }

    Syncmer next();

private:
    const std::string_view seq;
    const size_t k;
    const size_t s;
    const size_t t;

    const uint64_t kmask = (1ULL << 2*k) - 1;
    const uint64_t smask = (1ULL << 2*s) - 1;
    const uint64_t kshift = (k - 1) * 2;
    const uint64_t sshift = (s - 1) * 2;
    std::deque<uint64_t> qs;  // s-mer hashes
    uint64_t qs_min_val = UINT64_MAX;
    size_t qs_min_pos = -1;
    size_t l = 0;
    uint64_t xk[2] = {0, 0};
    uint64_t xs[2] = {0, 0};
    size_t i = 0;
};

class RandstrobeIterator2 {
public:
    RandstrobeIterator2(
        const std::string& seq, size_t k, size_t s, size_t t,
        unsigned w_min,
        unsigned w_max,
        uint64_t q,
        int max_dist
    ) : syncmer_iterator(SyncmerIterator(seq, k, s, t))
      , w_min(w_min)
      , w_max(w_max)
      , q(q)
      , max_dist(max_dist)
    { }

    Randstrobe next();
    Randstrobe end() const { return Randstrobe{0, 0, 0}; }

private:
    SyncmerIterator syncmer_iterator;
    const unsigned w_min;
    const unsigned w_max;
    const uint64_t q;
    const unsigned int max_dist;
    std::deque<Syncmer> syncmers;
};


std::pair<std::vector<syncmer_hash_t>, std::vector<unsigned int>> make_string_to_hashvalues_open_syncmers_canonical(
    const std::string_view seq,
    const size_t k,
    const size_t s,
    const size_t t
);

#endif
