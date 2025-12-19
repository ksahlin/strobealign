#include <cstddef>
#include <tuple>
#include <cassert>
#include <ostream>
#include "alignment.hpp"


/*
 * Find highest-scoring segment between reference and query assuming only matches
 * and mismatches are allowed.
 *
 * The end_bonus is added to the score if the segment extends until the end
 * of the query, once for each end.
 */
std::tuple<size_t, size_t, int> highest_scoring_segment(
    const std::string& query, const std::string& ref, int match, int mismatch, int end_bonus
) {
    size_t n = query.length();

    size_t start = 0; // start of the current segment
    int score = end_bonus; // accumulated score so far in the current segment

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
    if (score + end_bonus > best_score) {
        best_score = score + end_bonus;
        best_end = query.length();
        best_start = start;
    }
    return std::make_tuple(best_start, best_end, best_score);
}

AlignmentInfo hamming_align(
    const std::string &query, const std::string &ref, int match, int mismatch, int end_bonus
) {
    AlignmentInfo aln;
    if (query.length() != ref.length()) {
        return aln;
    }

    auto [segment_start, segment_end, score] = highest_scoring_segment(query, ref, match, mismatch, end_bonus);

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

    int soft_right = query.length() - segment_end;
    if (soft_right > 0) {
        cigar.push(CIGAR_SOFTCLIP, soft_right);
    }

    aln.cigar = std::move(cigar);
    aln.sw_score = score;
    aln.edit_distance = mismatches;
    aln.ref_start = segment_start;
    aln.ref_end = segment_end;
    aln.query_start = segment_start;
    aln.query_end = segment_end;
    return aln;
}

AlignmentInfo hamming_align_global(
    const std::string_view& query, const std::string_view& ref, uint match, uint mismatch
) {
    AlignmentInfo result;
    
    if (query.size() != ref.size()) {
        return result;
    }
    
    const size_t len = query.size();
    result.edit_distance = 0;
    result.query_start = 0;
    result.query_end = len;
    result.ref_start = 0;
    result.ref_end = len;
    
    for (size_t i = 0; i < len; ++i) {
        if (query[i] == ref[i]) {
            result.sw_score += match;
            result.cigar.push(CIGAR_EQ, 1);
        } else {
            result.sw_score -= mismatch;
            result.edit_distance += 1;
            result.cigar.push(CIGAR_X, 1);
        }
    }
    
    return result;
}

std::ostream& operator<<(std::ostream& os, const AlignmentParameters& params) {
    os
        << "AlignmentParameters("
        << "match=" << params.match
        << ", mismatch=" << params.mismatch
        << ", gap_open=" << params.gap_open
        << ", gap_extend=" << params.gap_extend
        << ", end_bonus=" << params.end_bonus
        << ")";
    return os;
}

