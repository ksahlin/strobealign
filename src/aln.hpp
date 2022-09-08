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
//#include "ksw2.h"
#include "ssw_cpp.h"


//#include <chrono>
//#include <thread>

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

void align_SE_read(std::thread::id thread_id, KSeq &record1, std::string &sam_out, logging_variables &log_vars, alignment_params &aln_params,
                   mapping_params &map_param, std::vector<unsigned int> &ref_lengths, std::vector<std::string> &ref_seqs, kmer_lookup &mers_index, mers_vector &flat_vector, idx_to_acc &acc_map );


#endif // aln_hpp_
