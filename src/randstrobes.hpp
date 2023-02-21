#ifndef RANDSTROBES_H
#define RANDSTROBES_H

#include <vector>
#include <string>
#include <tuple>
#include <deque>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <inttypes.h>

using syncmer_hash_t = uint64_t;
using randstrobe_hash_t = uint64_t;

// only used during index generation
struct RefRandstrobeWithHash {
    using packed_t = uint32_t;
    randstrobe_hash_t hash;
    uint32_t position;
    packed_t packed; // packed representation of ref_index and strobe offset

    bool operator< (const RefRandstrobeWithHash& other) const {
        return hash < other.hash;
    }
};

struct QueryRandstrobe {
    randstrobe_hash_t hash;
    unsigned int start;
    unsigned int end;
    bool is_reverse;
};

using QueryRandstrobeVector = std::vector<QueryRandstrobe>;

QueryRandstrobeVector randstrobes_query(int k, unsigned w_min, unsigned w_max, const std::string &seq, int s, int t, uint64_t q, int max_dist);

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
        const std::vector<uint64_t> &string_hashes,
        const std::vector<unsigned int> &pos_to_seq_coordinate,
        unsigned w_min,
        unsigned w_max,
        uint64_t q,
        int max_dist
    ) : string_hashes(string_hashes)
      , pos_to_seq_coordinate(pos_to_seq_coordinate)
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
        return get(strobe1_start++);
    }

    bool has_next() {
        return strobe1_start + w_min < string_hashes.size();
    }

private:
    Randstrobe get(unsigned int strobe1_start) const;
    const std::vector<uint64_t> &string_hashes;
    const std::vector<unsigned int> &pos_to_seq_coordinate;
    const unsigned w_min;
    const unsigned w_max;
    const uint64_t q;
    const unsigned int max_dist;
    unsigned int strobe1_start = 0;
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
    SyncmerIterator(const std::string& seq, size_t k, size_t s, size_t t)
        : seq(seq), k(k), s(s), t(t) { }

    Syncmer next();

private:
    const std::string& seq;
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
    const std::string &seq,
    const size_t k,
    const size_t s,
    const size_t t
);

#endif
