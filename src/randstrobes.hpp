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

static constexpr uint64_t RANDSTROBE_HASH_MASK = 0xFFFFFFFFFFFFFE00;

struct RefRandstrobe {
private:
    // packed representation of hash, offset and first_strobe_is_main
    randstrobe_hash_t m_hash_offset_flag;
    uint32_t m_position;
    uint32_t m_ref_index;

public:
    RefRandstrobe() : m_hash_offset_flag(0), m_position(0), m_ref_index(0) { }

    RefRandstrobe(randstrobe_hash_t hash, uint32_t position, uint32_t ref_index, uint8_t offset, bool first_strobe_is_main)
        : m_hash_offset_flag((hash & RANDSTROBE_HASH_MASK) | (offset << 1) | first_strobe_is_main)
        , m_position(position)
        , m_ref_index(ref_index)
    { }

    bool operator< (const RefRandstrobe& other) const {
        // Compare both hash and position to ensure that the order of the
        // RefRandstrobes in the index is reproducible no matter which sorting
        // function is used. This branchless comparison is faster than the
        // equivalent one using std::tie.
        __uint128_t lhs = (static_cast<__uint128_t>(m_hash_offset_flag) << 64) | ((static_cast<uint64_t>(m_position) << 32) | m_ref_index);
        __uint128_t rhs = (static_cast<__uint128_t>(other.m_hash_offset_flag) << 64) | ((static_cast<uint64_t>(other.m_position) << 32) | m_ref_index);
        return lhs < rhs;
    }

    bool first_strobe_is_main() const {
        return m_hash_offset_flag & 1;
    }

    unsigned reference_index() const {
        return m_ref_index;
    }

    unsigned strobe2_offset() const {
        return (m_hash_offset_flag >> 1) & 0xff;
    }

    randstrobe_hash_t hash() const {
        return m_hash_offset_flag & RANDSTROBE_HASH_MASK;
    }

    uint32_t position() const {
        return m_position;
    }

    static constexpr size_t max_number_of_references = (1ul << 32) - 1;
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
    bool is_revcomp;
};

std::ostream& operator<<(std::ostream& os, const QueryRandstrobe& randstrobe);

std::vector<QueryRandstrobe> randstrobes_query(const std::string_view seq, const IndexParameters& parameters);

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

    bool operator==(const Syncmer& rhs) const {
        return this->hash == rhs.hash && this->position == rhs.position;
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
    ) : syncmers(syncmers), parameters(parameters)
    {
        if (parameters.w_min > parameters.w_max) {
            throw std::invalid_argument("w_min is greater than w_max");
        }
    }

    Randstrobe next() {
        return get(strobe1_index++);
    }

    bool has_next() {
        return strobe1_index + parameters.w_min < syncmers.size();
    }

private:
    Randstrobe get(unsigned int strobe1_index) const;
    const std::vector<Syncmer>& syncmers;
    const RandstrobeParameters parameters;

    unsigned strobe1_index = 0;
};

std::ostream& operator<<(std::ostream& os, const Syncmer& syncmer);

class SyncmerIterator {
public:
    SyncmerIterator(const std::string_view seq, SyncmerParameters parameters)
        : seq(seq), parameters(parameters) { }

    Syncmer next();

private:
    const std::string_view seq;
    const SyncmerParameters parameters;

    const uint64_t kmask = (1ULL << 2*parameters.k) - 1;
    const uint64_t smask = (1ULL << 2*parameters.s) - 1;
    const uint64_t kshift = (parameters.k - 1) * 2;
    const uint64_t sshift = (parameters.s - 1) * 2;
    std::deque<uint64_t> qs;  // s-mer hashes
    uint64_t qs_min_val = UINT64_MAX;
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
      , parameters(randstrobe_parameters)
    { }

    Randstrobe next();
    Randstrobe end() const { return Randstrobe{0, 0, 0, false}; }

private:
    SyncmerIterator syncmer_iterator;
    const RandstrobeParameters parameters;
    std::deque<Syncmer> syncmers;
};


std::vector<Syncmer> canonical_syncmers(const std::string_view seq, SyncmerParameters parameters);

#endif
