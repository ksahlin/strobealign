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

aln_info Aligner::align(const std::string &ref, const std::string &query) const {
    m_align_calls++;
    aln_info aln;
    int32_t maskLen = query.length() / 2;
    maskLen = std::max(maskLen, 15);
    if (ref.length() > 2000){
//        std::cerr << "ALIGNMENT TO REF LONGER THAN 2000bp - REPORT TO DEVELOPER. Happened for read: " <<  query << " ref len:" << ref.length() << std::endl;
        aln.global_ed = 100000;
        aln.ed = 100000;
        aln.ref_offset = 0;
        aln.cigar = "*";
        aln.sw_score = -1000000;
        return aln;
    }

    StripedSmithWaterman::Alignment alignment_ssw;

    // query must be NULL-terminated
    ssw_aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment_ssw, maskLen, 1);

    aln.global_ed = alignment_ssw.global_ed;
    aln.ed = alignment_ssw.mismatches;
    aln.ref_offset = alignment_ssw.ref_begin;
    aln.cigar = alignment_ssw.cigar_string;
    aln.sw_score = alignment_ssw.sw_score;
    // ref_end is a 1-based position
    aln.length = alignment_ssw.ref_end - alignment_ssw.ref_begin + 1;
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
    const std::string &query, const std::string &ref, int match, int mismatch, int &soft_left, int &soft_right
) {
    aln_info aln;
    if (query.length() != ref.length()) {
        return aln;
    }

    auto [segment_start, segment_end] = highest_scoring_segment(query, ref, match, mismatch);

    std::stringstream cigar;
    if (segment_start > 0) {
        cigar << segment_start << 'S';
    }

    // Create CIGAR string and count mismatches
    int counter = 0;
    bool prev_is_match = false;
    int hamming_mod = 0;
    bool first = true;
    for (size_t i = segment_start; i < segment_end; i++) {
        bool is_match = query[i] == ref[i];
        hamming_mod += is_match ? 0 : 1;
        if (!first && is_match != prev_is_match) {
            cigar << counter << (prev_is_match ? '=' : 'X');
            counter = 0;
        }
        counter++;
        prev_is_match = is_match;
        first = false;
    }
    if (!first) {
        cigar << counter << (prev_is_match ? '=' : 'X');
    }
    int aln_score = (segment_end - segment_start - hamming_mod) * match - hamming_mod * mismatch;

    soft_left = segment_start;
    soft_right = query.length() - segment_end;
    if (soft_right > 0) {
        cigar << query.length() - segment_end << 'S';
    }

    aln.cigar = cigar.str();
    aln.sw_score = aln_score;
    aln.ed = hamming_mod;
    aln.global_ed = aln.ed + soft_left + soft_right;
    aln.length = query.length() - soft_left - soft_right;
    aln.ref_offset = soft_left;
    return aln;
}
