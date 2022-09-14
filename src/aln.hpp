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


void align_PE_read(klibpp::KSeq &record1, klibpp::KSeq &record2, std::string &sam_out, logging_variables &log_vars, i_dist_est &isize_est, alignment_params &aln_params,
                          mapping_params &map_param, std::vector<unsigned int> &ref_lengths, std::vector<std::string> &ref_seqs, kmer_lookup &mers_index, mers_vector &flat_vector, ref_names &acc_map );

void align_SE_read(klibpp::KSeq& record, std::string& outstring, logging_variables& log_vars, alignment_params& aln_params, mapping_params& map_param, std::vector< unsigned int >& ref_lengths, std::vector< std::string >& ref_seqs, kmer_lookup& mers_index, mers_vector& flat_vector, ref_names& acc_map );


#endif
