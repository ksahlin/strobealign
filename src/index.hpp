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
#include "xxhash.h"
#include "exceptions.hpp"
#include "refs.hpp"

uint64_t hash(const std::string& kmer);

typedef std::vector< uint64_t > hash_vector; //only used during index generation
typedef std::vector< std::tuple<uint32_t, int32_t >> mers_vector;
typedef std::vector< std::tuple<uint64_t, uint32_t, int32_t >> ind_mers_vector; //only used during index generation
//typedef std::vector< std::tuple<uint64_t, unsigned int, unsigned int>> mers_vector_reduced;
typedef robin_hood::unordered_map< uint64_t, std::tuple<unsigned int, unsigned int >> kmer_lookup;
typedef std::vector< std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, bool>> mers_vector_read;


struct st_index {
    st_index() : filter_cutoff(0) {}
    unsigned int filter_cutoff; //This also exists in mapping_params, but is calculated during index generation, 
                                //therefore stored here since it needs to be saved with the index.
    mers_vector flat_vector;
    kmer_lookup mers_index;
};


//mers_vector seq_to_kmers(int k, std::string &seq, unsigned int ref_index);
void seq_to_randstrobes2(ind_mers_vector& flat_vector, int n, int k, int w_min, int w_max, const std::string &seq, int ref_index, int s, int t, uint64_t q, int max_dist);
mers_vector_read seq_to_randstrobes2_read(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index, int s, int t, uint64_t q, int max_dist);
//mers_vector seq_to_randstrobes3(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index, int w);

void write_index(const st_index& index, const References& references, const std::string& filename);
void read_index(st_index& index, References& references, const std::string& filename);

uint64_t count_unique_elements(const hash_vector& h_vector);
unsigned int index_vector(const hash_vector& h_vector, kmer_lookup &mers_index, float f);

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


struct logging_variables {
    std::chrono::duration<double> tot_read_file;
    std::chrono::duration<double> tot_construct_strobemers;
    std::chrono::duration<double> tot_find_nams;
    std::chrono::duration<double> tot_time_rescue;
    std::chrono::duration<double> tot_find_nams_alt;
    std::chrono::duration<double> tot_sort_nams;
    std::chrono::duration<double> tot_extend;
    std::chrono::duration<double> tot_rc;
    std::chrono::duration<double> tot_write_file;

    unsigned int tot_ksw_aligned = 0;
    unsigned int tot_rescued = 0;
    unsigned int tot_all_tried = 0;
    unsigned int did_not_fit = 0;
    unsigned int tried_rescue = 0;

    logging_variables operator+=(const logging_variables& other) {
        this->tot_read_file += other.tot_read_file;
        this->tot_construct_strobemers += other.tot_construct_strobemers;
        this->tot_find_nams += other.tot_find_nams;
        this->tot_time_rescue += other.tot_time_rescue;
        this->tot_find_nams_alt += other.tot_find_nams_alt;
        this->tot_sort_nams += other.tot_sort_nams;
        this->tot_extend += other.tot_extend;
        this->tot_rc += other.tot_rc;
        this->tot_write_file += other.tot_write_file;
        this->tot_ksw_aligned += other.tot_ksw_aligned;
        this->tot_rescued += other.tot_rescued;
        this->tot_all_tried += other.tot_all_tried;
        this->did_not_fit += other.did_not_fit;
        this->tried_rescue += other.tried_rescue;
        return *this;
    }
};


struct i_dist_est {
    float sample_size = 1;
    float mu = 300;
    float sigma = 100;
    float V = 10000;
    float SSE = 10000;
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
    unsigned int filter_cutoff;
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


#endif /* index_hpp */



