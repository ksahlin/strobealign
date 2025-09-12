#ifndef BLOCK_ALIGNER_WRAPPER_H
#define BLOCK_ALIGNER_WRAPPER_H
#include <string>
#include <cstddef>
#include "block-aligner/c/block_aligner.h"
#include "aligner.hpp"
#include "cigar.hpp"

struct AlignmentResult {
    int score;
    size_t query_start;
    size_t query_end;
    size_t ref_start;
    size_t ref_end;
    Cigar cigar;
};

AlignmentResult global_alignment(const std::string_view& query, const std::string_view& ref, const AlignmentParameters& scoring_params);
AlignmentResult xdrop_alignment(const std::string_view& query, const std::string_view& ref, const AlignmentParameters& scoring_params, const bool reverse);

#endif
