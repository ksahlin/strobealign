//  Created by Kristoffer Sahlin on 4/21/21.

#ifndef STROBEALIGN_INDEX_HPP
#define STROBEALIGN_INDEX_HPP

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

/*
 * This describes where a randstrobe occurs. Info stored:
 * - reference index
 * - position of the first strobe
 * - offset of the second strobe
*/

class RefRandstrobe {
public:
    RefRandstrobe() { }  // TODO should not be needed
    RefRandstrobe(uint32_t position, uint32_t packed) : position(position), m_packed(packed) {
    }
    uint32_t position;

    int reference_index() const {
        return m_packed >> bit_alloc;
    }

    int strobe2_offset() const {
        return m_packed & mask;
    }

    RefRandstrobeWithHash::packed_t packed() const {
        return m_packed;
    }

private:
    static const int bit_alloc = 8;
    static const int mask = (1 << bit_alloc) - 1;
    RefRandstrobeWithHash::packed_t m_packed;
};

using RefRandstrobeVector = std::vector<RefRandstrobe>;


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
    mutable IndexCreationStatistics stats;

    void write(const std::string& filename) const;
    void read(const std::string& filename);
    void populate(float f, size_t n_threads);
    void print_diagnostics(const std::string& logfile_name, int k) const;
    unsigned int find(uint64_t key) const;

    uint64_t get_hash(unsigned int position) const {
        if (position < randstrobes_vector.size()){
            return randstrobes_vector[position].hash;
        }else{
            return -1;
        }
    }

    unsigned int get_strobe1_position(unsigned int position) const {
        return randstrobes_vector[position].position;
    }

    int strobe2_offset(unsigned int position) const {
        return randstrobes_vector[position].packed & mask;
    }

    int reference_index(unsigned int position) const {
        return randstrobes_vector[position].packed >> bit_alloc;
    }

    unsigned int get_count(const unsigned int position) const {
        constexpr unsigned int MAX_LINEAR_SEARCH = 8;
        const auto hash = randstrobes_vector[position].hash;

        // For 95% of cases, the result will be small and a brute force search
        // is the best option. Once, we go over MAX_LINEAR_SEARCH, though, we
        // call get_count_line_search which performs a line search with
        // expanding step in O(log N) comparisons
        //
        unsigned int count = 1;
        while (get_hash(position + count) == hash
                && count < MAX_LINEAR_SEARCH) {
            ++count;
        }
        if (get_hash(position + count) != hash) {
            return count;
        }
        return count + get_count_line_search(position + count);
    }

    int k() const {
        return parameters.k;
    }

private:
    void add_randstrobes_to_vector(int randstrobe_hashes);
    unsigned int get_count_line_search(unsigned int) const;

    const IndexParameters& parameters;
    const References& references;
    std::vector<RefRandstrobeWithHash> randstrobes_vector;
    std::vector<unsigned int> hash_positions;
    static constexpr int bit_alloc = 8;
    static constexpr int mask = (1 << bit_alloc) - 1;
};

#endif
