#ifndef ALIGNMENT_HPP
#define ALIGNMENT_HPP

#include <string_view>
#include "cigar.hpp"

struct AlignmentParameters {
    // match is a score, the others are penalties (all are nonnegative)
    int match;
    int mismatch;
    int gap_open;
    int gap_extend;
    int end_bonus;

    // temporarly here for piecewise
    int x_drop_threshold;
    uint min_block;
    uint max_block;
};

std::ostream& operator<<(std::ostream& os, const AlignmentParameters& params);

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

inline int hamming_distance(const std::string &s, const std::string &t) {
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
    const std::string& query, const std::string& ref, int match, int mismatch, int end_bonus
);

AlignmentInfo hamming_align(
    const std::string &query, const std::string &ref, int match, int mismatch, int end_bonus
);

AlignmentInfo hamming_align_global(
    const std::string_view& query, const std::string_view& ref, uint match, uint mismatch
);

#endif
