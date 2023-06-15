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


struct IndexCreationStatistics {
    uint64_t tot_strobemer_count = 0;
    uint64_t tot_occur_once = 0;
    uint64_t tot_high_ab = 0;
    uint64_t tot_mid_ab = 0;
    uint64_t index_cutoff = 0;
    uint64_t filter_cutoff = 0;
    uint64_t unique_strobemers = 0;

    std::chrono::duration<double> elapsed_hash_index;
    std::chrono::duration<double> elapsed_generating_seeds;
    std::chrono::duration<double> elapsed_counting_hashes;
    std::chrono::duration<double> elapsed_sorting_seeds;
};

struct StrobemerIndex {
    using bucket_index_t = uint64_t;
    StrobemerIndex(const References& references, const IndexParameters& parameters, int bits=-1)
        : filter_cutoff(0)
        , parameters(parameters)
        , references(references)
        , bits(bits == -1 ? pick_bits(references.total_length()) : bits)
    {
        if (this->bits < 8 || this->bits > 31) {
            throw BadParameter("Bits must be between 8 and 31");
        }
    }
    unsigned int filter_cutoff; //This also exists in mapping_params, but is calculated during index generation,
                                //therefore stored here since it needs to be saved with the index.
    mutable IndexCreationStatistics stats;

    void write(const std::string& filename) const;
    void read(const std::string& filename);
    void populate(float f, size_t n_threads);
    void print_diagnostics(const std::string& logfile_name, int k) const;
    int pick_bits(size_t size) const;
    size_t find(randstrobe_hash_t key) const {
        constexpr int MAX_LINEAR_SEARCH = 4;
        const unsigned int top_N = key >> (64 - bits);
        bucket_index_t position_start = randstrobe_start_indices[top_N];
        bucket_index_t position_end = randstrobe_start_indices[top_N + 1];
        if (position_start == position_end) {
            return end();
        }

        if (position_end - position_start < MAX_LINEAR_SEARCH) {
            for ( ; position_start < position_end; ++position_start) {
                if (randstrobes[position_start].hash == key) return position_start;
                if (randstrobes[position_start].hash > key) return end();
            }
            return end();
        }
        auto cmp = [](const RefRandstrobe lhs, const RefRandstrobe rhs) {return lhs.hash < rhs.hash; };

        auto pos = std::lower_bound(randstrobes.begin() + position_start,
                                               randstrobes.begin() + position_end,
                                               RefRandstrobe{key, 0, 0},
                                               cmp);
        if (pos->hash == key) return pos - randstrobes.begin();
        return end();
    }

    randstrobe_hash_t get_hash(bucket_index_t position) const {
        if (position < randstrobes.size()) {
            return randstrobes[position].hash;
        } else {
            return end();
        }
    }
    
    bool is_filtered(bucket_index_t position) const {
        return get_hash(position) == get_hash(position + filter_cutoff);
    }

    unsigned int get_strobe1_position(bucket_index_t position) const {
        return randstrobes[position].position;
    }

    int strobe2_offset(bucket_index_t position) const {
        return randstrobes[position].strobe2_offset();
    }

    int reference_index(bucket_index_t position) const {
        return randstrobes[position].reference_index();
    }

    unsigned int get_count(bucket_index_t position) const {
        // For 95% of cases, the result will be small and a brute force search
        // is the best option. Once, we go over MAX_LINEAR_SEARCH, though, we
        // use a binary search to get the next position
        // In the human genome, if we assume that the frequency 
        // a hash will be queried is proportional to the frequency it appears in the table, 
        // with MAX_LINEAR_SEARCH=8, the actual value will be 96%.

        constexpr unsigned int MAX_LINEAR_SEARCH = 8;
        const auto key = randstrobes[position].hash;
        const unsigned int top_N = key >> (64 - bits);
        bucket_index_t position_end = randstrobe_start_indices[top_N + 1];
        uint64_t count = 1;

        if (position_end - position < MAX_LINEAR_SEARCH) {
            for (bucket_index_t position_start = position + 1; position_start < position_end; ++position_start) {
                if (randstrobes[position_start].hash == key){
                    count += 1;
                }
                else{
                    break;
                }
            }
            return count;
        }
        auto cmp = [](const RefRandstrobe lhs, const RefRandstrobe rhs) {return lhs.hash < rhs.hash; };

        auto pos = std::upper_bound(randstrobes.begin() + position,
                                               randstrobes.begin() + position_end,
                                               RefRandstrobe{key, 0, 0},
                                               cmp);
        return (pos - randstrobes.begin() - 1) - position + 1;
    }

    int end() const {
        return -1;
    }

    int k() const {
        return parameters.k;
    }

    int get_bits() const {
        return bits;
    }

private:
    void add_randstrobes_to_vector();

    const IndexParameters& parameters;
    const References& references;

    /*
     * The randstrobes vector contains all randstrobes sorted by hash.
     *
     * The randstrobe_start_indices vector points to entries in the
     * randstrobes vector. randstrobe_start_indices[x] is the index of the
     * first entry in randstrobes whose top *bits* bits of its hash value are
     * greater than or equal to x.
     *
     * randstrobe_start_indices has one extra guard entry at the end that
     * is always randstrobes.size().
     */

    std::vector<RefRandstrobe> randstrobes;
    std::vector<bucket_index_t> randstrobe_start_indices;
    int bits; // no. of bits of the hash to use when indexing a randstrobe bucket
};

#endif
