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
    uint64_t distinct_strobemers = 0;

    std::chrono::duration<double> elapsed_hash_index;
    std::chrono::duration<double> elapsed_generating_seeds;
    std::chrono::duration<double> elapsed_counting_hashes;
    std::chrono::duration<double> elapsed_sorting_seeds;
};

struct StrobemerIndex {
    using bucket_index_t = uint64_t;
    StrobemerIndex(const References& references, const IndexParameters& parameters, int bits=-1)
        : filter_cutoff(0)
        , partial_filter_cutoff(0)
        , parameters(parameters)
        , references(references)
        , bits(bits == -1 ? pick_bits(references.total_length()) : bits)
    {
        if (this->bits < 8 || this->bits > 31) {
            throw BadParameter("Bits must be between 8 and 31");
        }
    }
    unsigned int filter_cutoff;
    unsigned int partial_filter_cutoff;
    mutable IndexCreationStatistics stats;

    void write(const std::string& filename) const;
    void read(const std::string& filename);
    void populate(float f, unsigned n_threads);
    void print_diagnostics(const std::string& logfile_name, int k) const;
    int pick_bits(size_t size) const;

    // Find first entry that matches the given key
    size_t find_full(randstrobe_hash_t key) const {
        return find(key, 0);
    }

    /*
     * Find the first entry that matches the main hash (ignoring the aux_len
     * least significant bits)
     */
    size_t find_partial(randstrobe_hash_t key) const {
        return find(key, parameters.randstrobe.aux_len);
    }

    /*
     * Find first entry whose hash matches the given key, but ignore the
     * b least significant bits
     */
    size_t find(randstrobe_hash_t key, uint8_t b) const {
        const unsigned int aux_len = b;
        randstrobe_hash_t key_prefix = key >> aux_len;

        constexpr int MAX_LINEAR_SEARCH = 4;
        const unsigned int top_N = key >> (64 - bits);
        bucket_index_t position_start = randstrobe_start_indices[top_N];
        bucket_index_t position_end = randstrobe_start_indices[top_N + 1];
        if (position_start == position_end) {
            return end();
        }

        if (position_end - position_start < MAX_LINEAR_SEARCH) {
            for ( ; position_start < position_end; ++position_start) {
                if (randstrobes[position_start].hash >> aux_len == key_prefix) return position_start;
                if (randstrobes[position_start].hash >> aux_len > key_prefix) return end();
            }
            return end();
        }
        auto cmp = [&aux_len](const RefRandstrobe lhs, const RefRandstrobe rhs) {
            return (lhs.hash >> aux_len) < (rhs.hash >> aux_len); };

        auto pos = std::lower_bound(randstrobes.begin() + position_start,
                                    randstrobes.begin() + position_end,
                                    RefRandstrobe{key, 0, 0},
                                    cmp);
        if (pos->hash >> aux_len == key_prefix) return pos - randstrobes.begin();
        return end();
    }

    randstrobe_hash_t get_hash(bucket_index_t position) const {
        if (position < randstrobes.size()) {
            return randstrobes[position].hash;
        } else {
            return end();
        }
    }

    randstrobe_hash_t get_main_hash(bucket_index_t position) const {
        if (position < randstrobes.size()) {
            return randstrobes[position].hash >> parameters.randstrobe.aux_len;
        } else {
            return end();
        }
    }

    bool first_strobe_is_main(bucket_index_t position) const {
        return randstrobes[position].first_strobe_is_main();
    }
    
    bool is_filtered(bucket_index_t position) const {
        return get_hash(position) == get_hash(position + filter_cutoff);
    }

    bool is_partial_filtered(bucket_index_t position) const {
        const unsigned int shift = parameters.randstrobe.aux_len;
        return (get_hash(position) >> shift) == (get_hash(position + partial_filter_cutoff) >> shift);
    }

    unsigned int get_strobe1_position(bucket_index_t position) const {
        return randstrobes[position].position;
    }

    int strobe2_offset(bucket_index_t position) const {
        return randstrobes[position].strobe2_offset();
    }

    std::pair<int, int> strobe_extent_partial(bucket_index_t position) const {
        // Construct the match from the strobe that was selected as the main part of the hash
        int ref_start = get_strobe1_position(position);
        if (!first_strobe_is_main(position)) {
            ref_start += strobe2_offset(position);
        }
        return {ref_start, ref_start + k()};
    }

    int reference_index(bucket_index_t position) const {
        return randstrobes[position].reference_index();
    }

    RefRandstrobe get_randstrobe(bucket_index_t position) const {
        return randstrobes[position];
    }

    size_t size() const {
        return randstrobes.size();
    }

    unsigned int get_count_full(bucket_index_t position) const {
        return get_count(position, 0);
    }

    unsigned int get_count_partial(bucket_index_t position) const {
        return get_count(position, parameters.randstrobe.aux_len);
    }

    unsigned int get_count(bucket_index_t position, uint8_t b) const {
        // For 95% of cases, the result will be small and a brute force search
        // is the best option. Once, we go over MAX_LINEAR_SEARCH, though, we
        // use a binary search to get the next position
        // In the human genome, if we assume that the frequency
        // a hash will be queried is proportional to the frequency it appears in the table,
        // with MAX_LINEAR_SEARCH=8, the actual value will be 96%.

        // Since the result depends on position, this function must be used on the smallest position which points to the
        // seed with the given hash to yield the number of seeds with this hash.

        constexpr unsigned int MAX_LINEAR_SEARCH = 8;
        const unsigned int aux_len = b;

        const auto key = randstrobes[position].hash;
        randstrobe_hash_t key_prefix = key >> aux_len;

        const unsigned int top_N = key >> (64 - bits);
        bucket_index_t position_end = randstrobe_start_indices[top_N + 1];
        uint64_t count = 1;

        if (position_end - position < MAX_LINEAR_SEARCH) {
            for (bucket_index_t position_start = position + 1; position_start < position_end; ++position_start) {
                if (randstrobes[position_start].hash >> aux_len == key_prefix){
                    count += 1;
                }
                else{
                    break;
                }
            }
            return count;
        }
        auto cmp = [&aux_len](const RefRandstrobe lhs, const RefRandstrobe rhs) {return (lhs.hash >> aux_len) < (rhs.hash >> aux_len); };

        auto pos = std::upper_bound(randstrobes.begin() + position,
                                    randstrobes.begin() + position_end,
                                    RefRandstrobe{key, 0, 0},
                                    cmp);
        return (pos - randstrobes.begin() - 1) - position + 1;
    }

    size_t end() const {
        return -1;
    }

    int k() const {
        return parameters.syncmer.k;
    }

    int get_bits() const {
        return bits;
    }

    int get_aux_len() const {
        return parameters.randstrobe.aux_len;
    }

private:
    void assign_all_randstrobes(const std::vector<uint64_t>& randstrobe_counts, size_t n_threads);
    void assign_randstrobes(size_t ref_index, size_t offset);

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
