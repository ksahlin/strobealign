#ifndef STROBEALIGN_ALIGNER_HPP
#define STROBEALIGN_ALIGNER_HPP

#include <string>
#include "ssw/ssw_cpp.h"


struct alignment_params {
    // match is a score, the others are penalties (all are nonnegative)
    int match;
    int mismatch;
    int gap_open;
    int gap_extend;
};

struct aln_info {
    std::string cigar;
    unsigned int ed{0};  // edit distance
    unsigned int ref_start{0};
    unsigned int ref_end{0};
    unsigned int query_start{0};
    unsigned int query_end{0};
    int sw_score{0};

    int ref_span() const { return ref_end - ref_start; }
};

struct Aligner {
public:
    Aligner(alignment_params parameters)
        : parameters(parameters)
        , ssw_aligner(StripedSmithWaterman::Aligner(parameters.match, parameters.mismatch, parameters.gap_open, parameters.gap_extend))
    { }

    aln_info align(const std::string &ref, const std::string &query) const;

    alignment_params parameters;

    unsigned calls_count() {
        return m_align_calls;
    }

private:
    const StripedSmithWaterman::Aligner ssw_aligner;
    const StripedSmithWaterman::Filter filter;
    mutable unsigned m_align_calls{0};  // no. of calls to the align() method
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

std::string compress_cigar(const std::string& ops);

std::pair<size_t, size_t> highest_scoring_segment(
    const std::string& query, const std::string& ref, int match, int mismatch
);

aln_info hamming_align(
    const std::string &query, const std::string &ref, int match, int mismatch
);

#endif
