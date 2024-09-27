#ifndef STROBEALIGN_ALN_HPP
#define STROBEALIGN_ALN_HPP

#include <string>
#include <vector>
#include <random>
#include "kseq++/kseq++.hpp"
#include "index.hpp"
#include "refs.hpp"
#include "sam.hpp"
#include "aligner.hpp"
#include "insertsizedistribution.hpp"

struct AlignmentStatistics {
    std::chrono::duration<double> tot_read_file{0};
    std::chrono::duration<double> tot_construct_strobemers{0};
    std::chrono::duration<double> tot_find_nams{0};
    std::chrono::duration<double> tot_time_rescue{0};
    std::chrono::duration<double> tot_find_nams_alt{0};
    std::chrono::duration<double> tot_sort_nams{0};
    std::chrono::duration<double> tot_extend{0};
    std::chrono::duration<double> tot_write_file{0};

    uint64_t n_reads{0};
    uint64_t tot_aligner_calls{0};
    uint64_t tot_rescued{0};
    uint64_t tot_all_tried{0};
    uint64_t inconsistent_nams{0};
    uint64_t nam_rescue{0};

    AlignmentStatistics operator+=(const AlignmentStatistics& other) {
        this->tot_read_file += other.tot_read_file;
        this->tot_construct_strobemers += other.tot_construct_strobemers;
        this->tot_find_nams += other.tot_find_nams;
        this->tot_time_rescue += other.tot_time_rescue;
        this->tot_find_nams_alt += other.tot_find_nams_alt;
        this->tot_sort_nams += other.tot_sort_nams;
        this->tot_extend += other.tot_extend;
        this->tot_write_file += other.tot_write_file;
        this->n_reads += other.n_reads;
        this->tot_aligner_calls += other.tot_aligner_calls;
        this->tot_rescued += other.tot_rescued;
        this->tot_all_tried += other.tot_all_tried;
        this->inconsistent_nams += other.inconsistent_nams;
        this->nam_rescue += other.nam_rescue;
        return *this;
    }

    AlignmentStatistics operator+=(const Details& details) {
        this->nam_rescue += details.nam_rescue;
        this->tot_rescued += details.mate_rescue;
        this->tot_all_tried += details.tried_alignment;
        this->inconsistent_nams += details.nam_inconsistent;

        return *this;
    }
};

enum class OutputFormat {
    SAM,
    PAF,
    Abundance
};

struct MappingParameters {
    int r { 150 };
    int max_secondary { 0 };
    float dropoff_threshold { 0.5 };
    int rescue_level { 2 };
    int max_tries { 20 };
    int rescue_cutoff;
    bool use_mcs{false};  // multi-context seeds
    OutputFormat output_format {OutputFormat::SAM};
    CigarOps cigar_ops{CigarOps::M};
    bool output_unmapped { true };
    bool details{false};
    bool fastq_comments{false};

    void verify() const {
        if (max_tries < 1) {
            throw BadParameter("max_tries must be greater than zero");
        }
    }
};

void align_or_map_paired(
    const klibpp::KSeq& record1,
    const klibpp::KSeq& record2,
    Sam& sam,
    std::string& outstring,
    AlignmentStatistics& statistics,
    InsertSizeDistribution& isize_est,
    const Aligner& aligner,
    const MappingParameters& map_param,
    const IndexParameters& index_parameters,
    const References& references,
    const StrobemerIndex& index,
    std::minstd_rand& random_engine,
    std::vector<double> &abundances
);

void align_or_map_single(
    const klibpp::KSeq& record,
    Sam& sam,
    std::string& outstring,
    AlignmentStatistics& statistics,
    const Aligner& aligner,
    const MappingParameters& map_param,
    const IndexParameters& index_parameters,
    const References& references,
    const StrobemerIndex& index,
    std::minstd_rand& random_engine,
    std::vector<double> &abundances
);

// Private declarations, only needed for tests

bool has_shared_substring(const std::string& read_seq, const std::string& ref_seq, int k);

#endif
