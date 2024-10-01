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
#include "statistics.hpp"


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
