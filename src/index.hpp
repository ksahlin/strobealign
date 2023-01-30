//
//  index.hpp
//  cpp_course
//
//  Created by Kristoffer Sahlin on 4/21/21.
//
#ifndef index_hpp
#define index_hpp

#include <chrono>
#include <stdio.h>
#include <string>
#include <vector>
#include <deque>
#include <tuple>
#include <cmath>
#include <iostream>
#include <cassert>
#include "robin_hood.h"
#include "exceptions.hpp"
#include "refs.hpp"
#include "randstrobes.hpp"
#include "indexparameters.hpp"

class ReferenceMer {
public:
    ReferenceMer() { }  // TODO should not be needed
    ReferenceMer(uint32_t position, uint32_t packed) : position(position), m_packed(packed) {
    }
    uint32_t position;

    int reference_index() const {
        return m_packed >> bit_alloc;
    }

    int strobe2_offset() const {
        return m_packed & mask;
    }

    MersIndexEntry::packed_t packed() const {
        return m_packed;
    }

private:
    static const int bit_alloc = 8;
    static const int mask = (1 << bit_alloc) - 1;
    MersIndexEntry::packed_t m_packed;
};

typedef std::vector<ReferenceMer> mers_vector;

class KmerLookupEntry {
public:
    KmerLookupEntry() { }
    KmerLookupEntry(unsigned int offset, unsigned int count) : m_offset(offset), m_count(count) { }

    unsigned int count() const {
        if (is_reference_mer()) {
            return 1;
        } else {
            return m_count;
        }
    }

    unsigned int offset() const{
        assert(!is_reference_mer());
        return m_offset;
    }

    bool is_reference_mer() const {
        return m_count & 0x8000'0000;
    }

    ReferenceMer as_reference_mer() const {
        assert(is_reference_mer());
        return ReferenceMer{m_offset, m_count & 0x7fff'ffff};
    }

    void set_count(unsigned int count) {
        m_count = count;
    }

    void set_offset(unsigned int offset) {
        assert(!is_reference_mer());
        m_offset = offset;
    }

private:
    unsigned int m_offset;
    unsigned int m_count;
};

using RandstrobeMap = robin_hood::unordered_map<uint64_t, KmerLookupEntry>;

typedef std::vector<uint64_t> hash_vector; //only used during index generation

struct IndexCreationStatistics {
    unsigned int flat_vector_size = 0;
    unsigned int tot_strobemer_count = 0;
    unsigned int tot_occur_once = 0;
    float frac_unique = 0;
    unsigned int tot_high_ab = 0;
    unsigned int tot_mid_ab = 0;
    unsigned int tot_distinct_strobemer_count = 0;
    unsigned int index_cutoff = 0;
    unsigned int filter_cutoff = 0;
    uint64_t unique_mers = 0;

    std::chrono::duration<double> elapsed_hash_index;
    std::chrono::duration<double> elapsed_generating_seeds;
    std::chrono::duration<double> elapsed_unique_hashes;
    std::chrono::duration<double> elapsed_sorting_seeds;
};

struct StrobemerIndex {
    StrobemerIndex(const References& references, const IndexParameters& parameters)
        : filter_cutoff(0)
        , parameters(parameters)
        , references(references) {}
    unsigned int filter_cutoff; //This also exists in mapping_params, but is calculated during index generation,
                                //therefore stored here since it needs to be saved with the index.
    mers_vector flat_vector;
    mutable IndexCreationStatistics stats;

    void write(const std::string& filename) const;
    void read(const std::string& filename);
    void populate(float f, size_t n_threads);
    void print_diagnostics(const std::string& logfile_name, int k) const;

    RandstrobeMap::const_iterator find(uint64_t key) const {
        return mers_index.find(key);
    }

    RandstrobeMap::const_iterator end() const {
        return mers_index.cend();
    }

    void add_entry(uint64_t key, unsigned int offset, unsigned int count) {
        KmerLookupEntry s{offset, count};
        mers_index[key] = s;
    }

private:
    ind_mers_vector add_randstrobes_to_hash_table();

    const IndexParameters& parameters;
    const References& references;
    RandstrobeMap mers_index; // k-mer -> (offset in flat_vector, occurence count )
};

#endif
