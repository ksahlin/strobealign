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
#include "robin_hood.h"
#include "exceptions.hpp"
#include "refs.hpp"
#include "randstrobes.hpp"

class ReferenceMer {
public:
    ReferenceMer() { }  // TODO should not be needed
    ReferenceMer(uint32_t position, int32_t packed) : position(position), packed(packed) {
    }
    uint32_t position;

    int reference_index() const {
        return packed >> bit_alloc;
    }

    int strobe2_offset() const {
        return packed & mask;
    }

private:
    const int bit_alloc = 8;
    const int mask = (1 << bit_alloc) - 1;
    int32_t packed;
};

typedef std::vector<ReferenceMer> mers_vector;

class KmerLookupEntry {
public:
    KmerLookupEntry() { }
    KmerLookupEntry(unsigned int offset, unsigned int count) : offset(offset), m_count(count) { }
    unsigned int offset;

    unsigned int count() const {
        return m_count;
    }
private:
    unsigned int m_count;
};

typedef robin_hood::unordered_map<uint64_t, KmerLookupEntry> kmer_lookup;

typedef std::vector<uint64_t> hash_vector; //only used during index generation

/* Settings that influence index creation */
class IndexParameters {
public:
    const int k;
    const int s;
    const int l;
    const int u;
    const uint64_t q;
    const int max_dist;
    const int t_syncmer;
    const int w_min;
    const int w_max;

    IndexParameters(int k, int s, int l, int u, int q, int max_dist)
        : k(k)
        , s(s)
        , l(l)
        , u(u)
        , q(q)
        , max_dist(max_dist)
        , t_syncmer((k - s) / 2 + 1)
        , w_min(std::max(1, k / (k - s + 1) + l))
        , w_max(k / (k - s + 1) + u) {
    }

    static IndexParameters from_read_length(int read_length, int c, int k = -1, int s = -1, int max_seed_len = -1);
    static IndexParameters read(std::istream& os);

    void write(std::ostream& os) const;
    bool operator==(const IndexParameters& other) const;
    bool operator!=(const IndexParameters& other) const { return !(*this == other); }

    void verify() const {
        if (k <= 7 || k > 32) {
            throw BadParameter("k not in [8,32]");
        }
        if (s > k) {
            throw BadParameter("s is larger than k");
        }
        if ((k - s) % 2 != 0) {
            throw BadParameter("(k - s) should be an even number to create canonical syncmers. Please set s to e.g. k-2, k-4, k-6, ...");
        }
        if (max_dist > 255) {
            throw BadParameter("maximum seed length (-m <max_dist>) is larger than 255");
        }
    }
};

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

    std::chrono::duration<double> elapsed_flat_vector;
    std::chrono::duration<double> elapsed_hash_index;
    std::chrono::duration<double> elapsed_generating_seeds;
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
    void populate(float f);
    void print_diagnostics(const std::string& logfile_name, int k) const;

    kmer_lookup::const_iterator find(uint64_t key) const {
        return mers_index.find(key);
    }

    kmer_lookup::const_iterator end() const {
        return mers_index.cend();
    }

    void add_entry(uint64_t key, unsigned int offset, unsigned int count) {
        KmerLookupEntry s{offset, count};
        mers_index[key] = s;
    }

private:
    const IndexParameters& parameters;
    const References& references;
    kmer_lookup mers_index; // k-mer -> (offset in flat_vector, occurence count )

    void index_vector(const hash_vector& h_vector, float f);
    ind_mers_vector generate_and_sort_seeds() const;
};


#endif
