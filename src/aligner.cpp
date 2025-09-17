#include <algorithm>
#include <cassert>
#include "piecewisealigner.hpp"
#include "aligner.hpp"

std::optional<AlignmentInfo> Aligner::align(const std::string &query, const std::string &ref) const {
    m_align_calls++;
    AlignmentInfo aln;
    int32_t maskLen = query.length() / 2;
    maskLen = std::max(maskLen, 15);
//     if (ref.length() > 2000){
// //        std::cerr << "ALIGNMENT TO REF LONGER THAN 2000bp - REPORT TO DEVELOPER. Happened for read: " <<  query << " ref len:" << ref.length() << std::endl;
//         return {};
//     }

    StripedSmithWaterman::Alignment alignment_ssw;

    // query must be NULL-terminated
    auto flag = ssw_aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment_ssw, maskLen);
    if (flag != 0 || alignment_ssw.ref_begin == -1) {
        return {};
    }

    aln.edit_distance = alignment_ssw.mismatches;
    aln.cigar = Cigar(alignment_ssw.cigar);
    aln.sw_score = alignment_ssw.sw_score;
    aln.ref_start = alignment_ssw.ref_begin;
    // end positions are off by 1 in SSW
    aln.ref_end = alignment_ssw.ref_end + 1;
    aln.query_start = alignment_ssw.query_begin;
    aln.query_end = alignment_ssw.query_end + 1;

    // Try to extend to beginning of the query to get an end bonus
    auto qstart = aln.query_start;
    auto rstart = aln.ref_start;
    auto score = aln.sw_score;
    auto edits = aln.edit_distance;
    Cigar front_cigar;
    while (qstart > 0 && rstart > 0) {
        qstart--;
        rstart--;
        if (query[qstart] == ref[rstart]) {
            score += parameters.match;
            front_cigar.push(CIGAR_EQ, 1);
        } else {
            score -= parameters.mismatch;
            front_cigar.push(CIGAR_X, 1);
            edits++;
        }
    }
    if (qstart == 0 && score + parameters.end_bonus > aln.sw_score) {
        if (aln.query_start > 0) {
            assert((aln.cigar.m_ops[0] & 0xF) == CIGAR_SOFTCLIP);
            aln.cigar.m_ops.erase(aln.cigar.m_ops.begin());  // remove soft clipping
            front_cigar.reverse();
            front_cigar += aln.cigar;
            aln.cigar = std::move(front_cigar);
        }
        aln.query_start = 0;
        aln.ref_start = rstart;
        aln.sw_score = score + parameters.end_bonus;
        aln.edit_distance = edits;
    }

    // Try to extend to end of query to get an end bonus
    auto qend = aln.query_end;
    auto rend = aln.ref_end;
    score = aln.sw_score;
    edits = aln.edit_distance;
    Cigar back_cigar;
    while (qend < query.length() && rend < ref.length()) {
        if (query[qend] == ref[rend]) {
            score += parameters.match;
            back_cigar.push(CIGAR_EQ, 1);
        } else {
            score -= parameters.mismatch;
            back_cigar.push(CIGAR_X, 1);
            edits++;
        }
        qend++;
        rend++;
    }
    if (qend == query.length() && score + parameters.end_bonus > aln.sw_score) {
        if (aln.query_end < query.length()) {
            assert((aln.cigar.m_ops[aln.cigar.m_ops.size() - 1] & 0xf) == CIGAR_SOFTCLIP);
            aln.cigar.m_ops.pop_back();
            aln.cigar += back_cigar;
        }
        aln.query_end = query.length();
        aln.ref_end = rend;
        aln.sw_score = score + parameters.end_bonus;
        aln.edit_distance = edits;
    }

    return aln;
}


AlignmentInfo Aligner::align_piecewise(
    const std::string &query,
    const std::string &ref,
    const std::vector<Anchor>& anchors,
    int padding
) const {
    m_align_calls++;
    AlignmentInfo info = piecewise_aligner.piecewise_extension_alignment(ref, query, anchors, padding);
    return info;
}

