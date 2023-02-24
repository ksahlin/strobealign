/*
 * Low-level alignment functions
 *
 * This is for anything that returns an aln_info object, currently
 * Aligner::align and hamming_align.
 */
#include <sstream>
#include <tuple>
#include <algorithm>
#include "aligner.hpp"

aln_info Aligner::align(const std::string &query, const std::string &ref) const {
    m_align_calls++;
    aln_info aln;
    int32_t maskLen = query.length() / 2;
    maskLen = std::max(maskLen, 15);
    if (ref.length() > 2000){
//        std::cerr << "ALIGNMENT TO REF LONGER THAN 2000bp - REPORT TO DEVELOPER. Happened for read: " <<  query << " ref len:" << ref.length() << std::endl;
        aln.edit_distance = 100000;
        aln.ref_start = 0;
        aln.sw_score = -1000000;
        return aln;
    }

    StripedSmithWaterman::Alignment alignment_ssw;

    // query must be NULL-terminated
    ssw_aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment_ssw, maskLen, 1);

    aln.edit_distance = alignment_ssw.mismatches;
    aln.cigar = Cigar(alignment_ssw.cigar);
    aln.sw_score = alignment_ssw.sw_score;
    aln.ref_start = alignment_ssw.ref_begin;
    // end positions are off by 1 in SSW
    aln.ref_end = alignment_ssw.ref_end + 1;
    aln.query_start = alignment_ssw.query_begin;
    aln.query_end = alignment_ssw.query_end + 1;
    return aln;
}

/*
 * Find highest-scoring segment between reference and query assuming only matches
 * and mismatches are allowed.
 */
std::pair<size_t, size_t> highest_scoring_segment(
    const std::string& query, const std::string& ref, int match, int mismatch
) {
    size_t n = query.length();

    size_t start = 0; // start of the current segment
    int score = 0; // accumulated score so far in the current segment

    size_t best_start = 0;
    size_t best_end = 0;
    int best_score = 0;
    for (size_t i = 0; i < n; ++i) {
        if (query[i] == ref[i]) {
            score += match;
        } else {
            score -= mismatch;
        }
        if (score < 0) {
            start = i + 1;
            score = 0;
        }
        if (score > best_score) {
            best_start = start;
            best_score = score;
            best_end = i + 1;
        }
    }
    return std::make_pair(best_start, best_end);
}

aln_info hamming_align(
    const std::string &query, const std::string &ref, int match, int mismatch
) {
    aln_info aln;
    if (query.length() != ref.length()) {
        return aln;
    }

    auto [segment_start, segment_end] = highest_scoring_segment(query, ref, match, mismatch);

    Cigar cigar;
    if (segment_start > 0) {
        cigar.push(CIGAR_SOFTCLIP, segment_start);
    }

    // Create CIGAR string and count mismatches
    int counter = 0;
    bool prev_is_match = false;
    int mismatches = 0;
    bool first = true;
    for (size_t i = segment_start; i < segment_end; i++) {
        bool is_match = query[i] == ref[i];
        mismatches += is_match ? 0 : 1;
        if (!first && is_match != prev_is_match) {
            cigar.push(prev_is_match ? CIGAR_EQ : CIGAR_X, counter);
            counter = 0;
        }
        counter++;
        prev_is_match = is_match;
        first = false;
    }
    if (!first) {
        cigar.push(prev_is_match ? CIGAR_EQ : CIGAR_X, counter);
    }
    int aln_score = (segment_end - segment_start - mismatches) * match - mismatches * mismatch;

    int soft_right = query.length() - segment_end;
    if (soft_right > 0) {
        cigar.push(CIGAR_SOFTCLIP, soft_right);
    }

    aln.cigar = std::move(cigar);
    aln.sw_score = aln_score;
    aln.edit_distance = mismatches;
    aln.ref_start = segment_start;
    aln.ref_end = segment_end;
    aln.query_start = segment_start;
    aln.query_end = segment_end;
    return aln;
}
