#ifndef RANDSTROBES_H
#define RANDSTROBES_H

#include <vector>
#include <string>
#include <tuple>
#include <deque>
#include <algorithm>
#include <iostream>
#include <inttypes.h>

// only used during index generation
struct MersIndexEntry {

    using packed_t = uint32_t;
    uint64_t hash;
    uint32_t position;
    packed_t packed; // packed representation of ref_index and strobe offset

    bool operator< (const MersIndexEntry& other) const {
        return hash < other.hash;
    }
};

typedef std::vector<MersIndexEntry> ind_mers_vector;


struct QueryMer {
    uint64_t hash;
    unsigned int start;
    unsigned int end;
    bool is_reverse;
};


typedef std::vector<QueryMer> mers_vector_read;

void randstrobes_reference(ind_mers_vector& flat_vector, int k, int w_min, int w_max, const std::string &seq, int ref_index, int s, int t, uint64_t q, int max_dist);
mers_vector_read randstrobes_query(int k, int w_min, int w_max, const std::string &seq, int s, int t, uint64_t q, int max_dist);

struct Randstrobe {
    uint64_t hash;
    unsigned int strobe1_pos;
    unsigned int strobe2_pos;
};

std::ostream& operator<<(std::ostream& os, const Randstrobe& randstrobe);

class RandstrobeIterator {
public:
    RandstrobeIterator(
        const std::vector<uint64_t> &string_hashes,
        const std::vector<unsigned int> &pos_to_seq_coordinate,
        int w_min,
        int w_max,
        uint64_t q,
        int max_dist
    ) : string_hashes(string_hashes)
      , pos_to_seq_coordinate(pos_to_seq_coordinate)
      , w_min(w_min)
      , w_max(w_max)
      , q(q)
      , max_dist(max_dist)
    {
    }

    Randstrobe next() {
        return get(strobe1_start++);
    }

    bool has_next() {
        return (strobe1_start + w_max < string_hashes.size())
            || (strobe1_start + w_min + 1 < string_hashes.size());
    }

private:
    Randstrobe get(unsigned int strobe1_start) const;
    const std::vector<uint64_t> &string_hashes;
    const std::vector<unsigned int> &pos_to_seq_coordinate;
    const int w_min;
    const int w_max;
    const uint64_t q;
    const unsigned int max_dist;
    unsigned int strobe1_start = 0;
};

struct Syncmer {
    uint64_t hash;
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

std::pair<std::vector<uint64_t>, std::vector<unsigned int>> make_string_to_hashvalues_open_syncmers_canonical(
    const std::string &seq,
    const size_t k,
    const size_t s,
    const size_t t
);

#endif
