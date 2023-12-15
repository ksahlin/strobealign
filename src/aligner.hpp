#ifndef STROBEALIGN_ALIGNER_HPP
#define STROBEALIGN_ALIGNER_HPP

#include <string>
#include <tuple>
#include "ssw/ssw_cpp.h"
#include "cigar.hpp"


struct AlignmentParameters {
    // match is a score, the others are penalties (all are nonnegative)
    int match;
    int mismatch;
    int gap_open;
    int gap_extend;
    int end_bonus;
};

struct AlignmentInfo {
    Cigar cigar;
    unsigned int edit_distance{0};
    unsigned int ref_start{0};
    unsigned int ref_end{0};
    unsigned int query_start{0};
    unsigned int query_end{0};
    int sw_score{0};

    int ref_span() const { return ref_end - ref_start; }
};

std::ostream& operator<<(std::ostream& os, const AlignmentInfo& info);

void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b);

struct Aligner {
public:
    Aligner(AlignmentParameters parameters)
        : parameters(parameters)
        , ssw_aligner(StripedSmithWaterman::Aligner(parameters.match, parameters.mismatch, parameters.gap_open, parameters.gap_extend))
    {
        ksw_gen_simple_mat(ksw_matrix_m, ksw_matrix, parameters.match, -parameters.mismatch);
    }

    AlignmentInfo align(const std::string& query, const std::string_view ref) const;
    AlignmentInfo ksw_extend(const std::string& query, const std::string& ref, bool right_align) const;

    AlignmentParameters parameters;

    unsigned calls_count() {
        return m_align_calls;
    }

private:
    const StripedSmithWaterman::Aligner ssw_aligner;
    const StripedSmithWaterman::Filter filter;
    mutable unsigned m_align_calls{0};  // no. of calls to the align() method
    const int8_t ksw_matrix_m{5};
    int8_t ksw_matrix[25];
};

inline int hamming_distance(const std::string& s, const std::string_view t) {
    if (s.length() != t.length()){
        return -1;
    }

    int mismatches = 0;
    for (size_t i = 0; i < s.length(); i++) {
        if (s[i] != t[i]) {
            mismatches++;
        }
    }

    return mismatches;
}

std::tuple<size_t, size_t, int> highest_scoring_segment(
    const std::string& query, const std::string_view ref, int match, int mismatch, int end_bonus
);

AlignmentInfo hamming_align(
    const std::string& query, const std::string_view ref, int match, int mismatch, int end_bonus
);

#endif
