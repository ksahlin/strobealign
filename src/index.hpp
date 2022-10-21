//
//  index.hpp
//  cpp_course
//
//  Created by Kristoffer Sahlin on 4/21/21.
//

#ifndef index_hpp
#define index_hpp

#include <chrono>  // for high_resolution_clock
#include <stdio.h>
#include <string>
#include <vector>
#include <deque>
#include <tuple>
#include "robin_hood.h"
#include "exceptions.hpp"
#include "refs.hpp"
#include "randstrobes.hpp"

struct ReferenceMer {
    uint32_t position;
    int32_t packed;
};

typedef std::vector<ReferenceMer> mers_vector;

struct KmerLookupEntry {
    unsigned int offset;
    unsigned int count;
};

typedef robin_hood::unordered_map<uint64_t, KmerLookupEntry> kmer_lookup;


struct hit {
    int query_s;
    int query_e;
    int ref_s;
    int ref_e;
    bool is_rc = false;
};

// Non-overlapping approximate match
struct nam {
    int nam_id;
    int query_s;
    int query_e;
    int query_prev_hit_startpos;
    int ref_s;
    int ref_e;
    int ref_prev_hit_startpos;
    int n_hits = 0;
    int ref_id;
    float score;
//    unsigned int previous_query_start;
//    unsigned int previous_ref_start;
    bool is_rc = false;
};

//struct aln_info {
//    std::string cigar;
//    unsigned int ed;
//    unsigned int ref_offset;
//    int sw_score;
//    int global_ed;
//    int length;
//};

struct alignment_params {
    int match;
    int mismatch;
    int gap_open;
    int gap_extend;
};


class i_dist_est {
public:
    float sample_size = 1;
    float mu = 300;
    float sigma = 100;
    float V = 10000;
    float SSE = 10000;

    // Add a new observation
    void update(int dist);
};

struct mapping_params {
    uint64_t q;
    int n;
    int k;
    int s;
    int t_syncmer;
    int w_min;
    int w_max;
    int max_secondary;
    float dropoff_threshold;
    int r;
    int m;
    int l;
    int u;
    int c;
    float f;
    int S;
    int M;
    int R;
    int max_dist;
    int maxTries;
    int max_seed_len;
    int rescue_cutoff;
    bool is_sam_out;

    void verify(){
        if (k <= 7 || k > 32) {
            throw BadMappingParameter("k not in [8,32]");
        }
        if (s > k) {
            throw BadMappingParameter("s is larger than k");
        }
        if ((k - s) % 2 != 0) {
            throw BadMappingParameter("(k - s) should be an even number to create canonical syncmers. Please set s to e.g. k-2, k-4, k-6, ...");
        }
        if (max_dist > 255) {
            throw BadMappingParameter("maximum seed length (-m <max_dist>) is larger than 255");
        }
    }
};

struct StrobemerIndex {
    StrobemerIndex() : filter_cutoff(0) {}
    unsigned int filter_cutoff; //This also exists in mapping_params, but is calculated during index generation,
                                //therefore stored here since it needs to be saved with the index.
    mers_vector flat_vector;
    kmer_lookup mers_index; // k-mer -> (offset in flat_vector, occurence count )

    void write(const References& references, const std::string& filename) const;
    void read(References& references, const std::string& filename);
    void populate(const References& references, const mapping_params& map_param);
private:
    ind_mers_vector generate_seeds(const References& references, const mapping_params& map_param) const;

};


#endif /* index_hpp */
