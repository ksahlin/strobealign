#ifndef aln_hpp
#define aln_hpp

#include <string>
#include <vector>
#include "kseq++.hpp"
#include "index.hpp"
#include "refs.hpp"

struct aln_info {
    std::string cigar;
    unsigned int ed;
    unsigned int ref_offset;
    int sw_score;
    int global_ed;
    int length;
};


void align_PE_read(klibpp::KSeq& record1, klibpp::KSeq& record2, std::string& outstring, logging_variables& log_vars, i_dist_est& isize_est, alignment_params& aln_params, mapping_params& map_param, const References& references, kmer_lookup& mers_index, mers_vector& flat_vector);

void align_SE_read(klibpp::KSeq& record, std::string& outstring, logging_variables& log_vars, alignment_params& aln_params, mapping_params& map_param, const References& references, kmer_lookup& mers_index, mers_vector& flat_vector);


#endif
