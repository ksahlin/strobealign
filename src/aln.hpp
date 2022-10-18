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

struct AlignmentStatistics {
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

    AlignmentStatistics operator+=(const AlignmentStatistics& other) {
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


void align_PE_read(klibpp::KSeq& record1, klibpp::KSeq& record2, std::string& outstring, AlignmentStatistics& statistics, i_dist_est& isize_est, alignment_params& aln_params, mapping_params& map_param, const References& references, kmer_lookup& mers_index, mers_vector& flat_vector);

void align_SE_read(klibpp::KSeq& record, std::string& outstring, AlignmentStatistics& statistics, alignment_params& aln_params, mapping_params& map_param, const References& references, kmer_lookup& mers_index, mers_vector& flat_vector);


#endif
