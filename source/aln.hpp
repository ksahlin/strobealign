#ifndef aln_hpp
#define aln_hpp

// The declarations in the aln.hpp file

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <math.h>
#include <chrono>  // for high_resolution_clock
//#include <omp.h>
#include <zlib.h>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <inttypes.h>

#include "kseq++.hpp"
using namespace klibpp;
#include "robin_hood.h"
#include "index.hpp"
#include "ksw2.h"
#include "ssw_cpp.h"


//#include <chrono>
//#include <thread>


typedef robin_hood::unordered_map<int, std::string > idx_to_acc;
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

struct aln_info {
    std::string cigar;
    unsigned int ed;
    unsigned int ref_offset;
    int sw_score;
    int global_ed;
    int length;
};



void align_PE_read(std::thread::id thread_id, KSeq &record1, KSeq &record2, std::string &sam_out, logging_variables &log_vars, i_dist_est &isize_est, alignment_params &aln_params,
                          mapping_params &map_param, std::vector<unsigned int> &ref_lengths, std::vector<std::string> &ref_seqs, kmer_lookup &mers_index, mers_vector &flat_vector, idx_to_acc &acc_map );

inline aln_info ssw_align(std::string &ref, std::string &query, int read_len, int match_score, int mismatch_penalty, int gap_opening_penalty, int gap_extending_penalty);

//inline void align_PE(alignment_params &aln_params, std::string &sam_string, std::vector<nam> &all_nams1, std::vector<nam> &all_nams2, KSeq &record1, KSeq &record2,
//                            idx_to_acc &acc_map, int k, std::vector<unsigned int> &ref_len_map, std::vector<std::string> &ref_seqs, logging_variables  &log_vars, float dropoff,
//                            i_dist_est &isize_est,  int max_tries, int max_secondary);

//inline void align_SE(alignment_params &aln_params, std::string &sam_string, std::vector<nam> &all_nams, std::string &query_acc, idx_to_acc &acc_map, int k, int read_len,
//                     std::vector<unsigned int> &ref_len_map, std::vector<std::string> &ref_seqs, std::string &read, std::string &qual, logging_variables &log_vars, float dropoff,
//                     int max_tries );

//inline void align_SE_secondary_hits(alignment_params &aln_params, std::string &sam_string, std::vector<nam> &all_nams, std::string &query_acc, idx_to_acc &acc_map,
//                                           int k, int read_len, std::vector<unsigned int> &ref_len_map, std::vector<std::string> &ref_seqs, std::string &read, std::string &qual,
//                                           logging_variables &log_vars, float dropoff, int max_tries, int max_secondary );


//inline std::pair<float,int> find_nams_rescue(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> hits_fw,
//                                                    std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> hits_rc,
//                                                    std::vector<nam> &final_nams, robin_hood::unordered_map< unsigned int, std::vector<hit>> &hits_per_ref,
//                                                    mers_vector_read &query_mers, mers_vector &ref_mers, kmer_lookup &mers_index, int k, std::vector<std::string> &ref_seqs,
//                                                    std::string &read, unsigned int filter_cutoff );
//
//inline std::pair<float,int> find_nams(std::vector<nam> &final_nams, robin_hood::unordered_map< unsigned int, std::vector<hit>> &hits_per_ref,
//                                             mers_vector_read &query_mers, mers_vector &ref_mers, kmer_lookup &mers_index, int k, std::vector<std::string> &ref_seqs,
//                                             std::string &read, unsigned int filter_cutoff );

//inline void output_hits_paf_PE(std::string &paf_output, nam &n, std::string &query_acc, idx_to_acc &acc_map, int k, int read_len, std::vector<unsigned int> &ref_len_map);

inline void output_hits_paf(std::string &paf_output, std::vector<nam> &all_nams, std::string query_acc, idx_to_acc &acc_map, int k, int read_len, std::vector<unsigned int> &ref_len_map);

inline void get_best_map_location(std::vector<std::tuple<int,nam,nam>> joint_NAM_scores, std::vector<nam> &nams1, std::vector<nam> &nams2, i_dist_est &isize_est, nam &best_nam1,  nam &best_nam2 );


#endif // aln_hpp_