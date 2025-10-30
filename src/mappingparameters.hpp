#ifndef STROBEALIGN_MAPPINGPARAMETERS_HPP
#define STROBEALIGN_MAPPINGPARAMETERS_HPP

#include "mcsstrategy.hpp"
#include "sam.hpp"
#include "exceptions.hpp"

enum class OutputFormat {
    SAM,
    PAF,
    Abundance
};

struct ChainingParameters {
    int max_lookback;
    float diag_diff_penalty;
    float gap_length_penalty;
    float valid_score_threshold;
    uint max_ref_gap;
    float matches_weight;
};

struct MappingParameters {
    int r { 150 };
    int max_secondary { 0 };
    float dropoff_threshold { 0.5 };
    int rescue_level { 2 };
    int max_tries { 20 };
    int rescue_cutoff;
    McsStrategy mcs_strategy{McsStrategy::Rescue};
    OutputFormat output_format {OutputFormat::SAM};
    CigarOps cigar_ops{CigarOps::M};
    bool output_unmapped { true };
    bool details{false};
    bool fastq_comments{false};

    bool use_nams{false};
    ChainingParameters chaining_params;

    void verify() const {
        if (max_tries < 1) {
            throw BadParameter("max_tries must be greater than zero");
        }
    }
};

#endif
