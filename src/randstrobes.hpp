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
private:
    randstrobe_hash_t m_hash;
public:
    uint32_t position;
private:
    uint32_t m_packed; // packed representation of ref_index and strobe offset
public:

    RefRandstrobe() { }

    RefRandstrobe(randstrobe_hash_t hash, uint32_t position, uint32_t ref_index, uint8_t offset, bool first_strobe_is_main)
        : m_hash(hash)
        , position(position)
        , m_packed((ref_index << 9) | (first_strobe_is_main << 8) | offset) { }

    bool operator< (const RefRandstrobe& other) const {
        // Compare both hash and position to ensure that the order of the
        // RefRandstrobes in the index is reproducible no matter which sorting
        // function is used. This branchless comparison is faster than the
        // equivalent one using std::tie.
        __uint128_t lhs = (static_cast<__uint128_t>(m_hash) << 64) | ((static_cast<uint64_t>(position) << 32) | m_packed);
        __uint128_t rhs = (static_cast<__uint128_t>(other.m_hash) << 64) | ((static_cast<uint64_t>(other.position) << 32) | m_packed);
        return lhs < rhs;
    }

    bool first_strobe_is_main() const {
        return (m_packed >> bit_alloc) & 1;
    }

    int reference_index() const {
        return m_packed >> (bit_alloc + 1);
    }

    int strobe2_offset() const {
        return m_packed & mask;
    }

    randstrobe_hash_t hash() const {
        return m_hash;
    }

private:
    static constexpr int bit_alloc = 8;
    static constexpr int mask = (1 << bit_alloc) - 1;

public:
    static constexpr uint32_t max_number_of_references = (1 << (32 - bit_alloc - 1)) - 1; // bit_alloc - 1 because 1 bit to first_strobe_is_main()
};

struct QueryRandstrobe {
    randstrobe_hash_t hash;
    unsigned int start;
    unsigned int end;
    /* Start and end of the main syncmer (relevant if the randstrobe couldn’t
     * be found in the index and we fall back to a partial hit)
     */
    unsigned int partial_start;
    unsigned int partial_end;
    bool is_reverse;
};

std::ostream& operator<<(std::ostream& os, const QueryRandstrobe& randstrobe);

using QueryRandstrobeVector = std::vector<QueryRandstrobe>;

QueryRandstrobeVector randstrobes_query(const std::string_view seq, const IndexParameters& parameters);

struct Randstrobe {
    randstrobe_hash_t hash;
    unsigned int strobe1_pos;
    unsigned int strobe2_pos;
    bool first_strobe_is_main;

    bool operator==(const Randstrobe& other) const {
        return hash == other.hash && strobe1_pos == other.strobe1_pos && strobe2_pos == other.strobe2_pos;
    }

    bool operator!=(const Randstrobe& other) const {
        return !(*this == other);
    }
};

std::ostream& operator<<(std::ostream& os, const Randstrobe& randstrobe);

struct Syncmer {
    syncmer_hash_t hash;
    size_t position;
    bool is_end() const {
        return hash == 0 && position == 0;
    }
};

/*
 * Iterate over randstrobes using a pre-computed vector of syncmers
 */
class RandstrobeIterator {
public:
    RandstrobeIterator(
        const std::vector<Syncmer>& syncmers,
        RandstrobeParameters parameters
    ) : syncmers(syncmers)
      , w_min(parameters.w_min)
      , w_max(parameters.w_max)
      , q(parameters.q)
      , max_dist(parameters.max_dist)
      , aux_len(parameters.aux_len)
    {
        if (w_min > w_max) {
            throw std::invalid_argument("w_min is greater than w_max");
        }
    }

    Randstrobe next() {
        return get(strobe1_index++);
    }

    bool has_next() {
        return strobe1_index + w_min < syncmers.size();
    }

private:
    Randstrobe get(unsigned int strobe1_index) const;
    const std::vector<Syncmer>& syncmers;
    const unsigned w_min;
    const unsigned w_max;
    const uint64_t q;
    const unsigned int max_dist;
    const unsigned int aux_len;
    unsigned strobe1_index = 0;
};

std::ostream& operator<<(std::ostream& os, const Syncmer& syncmer);

class SyncmerIterator {
public:
    SyncmerIterator(const std::string_view seq, SyncmerParameters parameters)
        : seq(seq), k(parameters.k), s(parameters.s), t(parameters.t_syncmer) { }

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

/*
 * Iterate over randstrobes while generating syncmers on the fly
 *
 * Unlike RandstrobeIterator, this does not need a pre-computed vector
 * of syncmers and therefore uses less memory.
 */
class RandstrobeGenerator {
public:
    RandstrobeGenerator(
        const std::string& seq,
        SyncmerParameters syncmer_parameters,
        RandstrobeParameters randstrobe_parameters
    ) : syncmer_iterator(SyncmerIterator(seq, syncmer_parameters))
      , w_min(randstrobe_parameters.w_min)
      , w_max(randstrobe_parameters.w_max)
      , q(randstrobe_parameters.q)
      , max_dist(randstrobe_parameters.max_dist)
      , aux_len(randstrobe_parameters.aux_len)
    { }

    Randstrobe next();
    Randstrobe end() const { return Randstrobe{0, 0, 0, false}; }

private:
    SyncmerIterator syncmer_iterator;
    const unsigned w_min;
    const unsigned w_max;
    const uint64_t q;
    const unsigned int max_dist;
    const unsigned int aux_len;
    std::deque<Syncmer> syncmers;
};


std::vector<Syncmer> canonical_syncmers(const std::string_view seq, SyncmerParameters parameters);

#endif
