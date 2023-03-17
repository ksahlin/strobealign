/*
 * Low-level alignment functions
 *
 * This is for anything that returns an aln_info object, currently
 * Aligner::align and hamming_align.
 *
 * ksw_extend code is based on https://github.com/lh3/ksw2/blob/master/cli.c
 */
#include <sstream>
#include <tuple>
#include <algorithm>
#include <cassert>
#include <cstring>  // memset
#include <iostream>
#include "ksw2/ksw2.h"
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

aln_info hamming_align(
    const std::string &query, const std::string &ref, int match, int mismatch, int end_bonus
) {
    aln_info aln;
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

namespace {

unsigned char seq_nt4_table[256] = {
    0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

}  // namespace

void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
{
    int i, j;
    a = a < 0? -a : a;
    b = b > 0? -b : b;
    for (i = 0; i < m - 1; ++i) {
        for (j = 0; j < m - 1; ++j)
            mat[i * m + j] = i == j? a : b;
        mat[i * m + m - 1] = 0;
    }
    for (j = 0; j < m; ++j)
        mat[(m - 1) * m + j] = 0;
}

std::ostream& operator<<(std::ostream& os, const ksw_extz_t& ez) {
    os << "ksw_extz_t("
        //
        << "\n max:             " << ez.max   // max overall score
        << "\n coord max_q:     " << ez.max_q  // max extension coordinate
        << "\n coord max_t:     " << ez.max_t  // max extension coordinate

        << "\n score mqe:       " << ez.mqe // max score when reaching the end of query
        << "\n mqe_t:           " << ez.mqe_t // coordinate in target corresponding to mqe

        << "\n score mte:       " << ez.mte  // max score when reaching the end of target
        << "\n mte_q:           " << ez.mte_q // coordinate in query corresponding to mte

        << "\n score both ends: " << ez.score  // max score reaching both ends
        << "\n cigar:           " << Cigar(ez.cigar, ez.n_cigar)
        << "\n zdropped:        " << ez.zdropped
        << "\n reach_end:       " << ez.reach_end
        << "\n)";
    return os;
}

aln_info Aligner::ksw_extend(const std::string& query, const std::string& ref, bool right_align) const {
    int w = -1; // band width; -1 is inf
    int zdrop = -1; // -1 to disable
    int flag = KSW_EZ_EXTZ_ONLY;
    if (right_align) {
        flag |= KSW_EZ_RIGHT;
    }
    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));

    ez.max_q = ez.max_t = ez.mqe_t = ez.mte_q = -1;
    ez.max = 0; ez.mqe = ez.mte = KSW_NEG_INF;
    ez.n_cigar = 0;
    int qlen = query.length();
    int tlen = ref.length();
    uint8_t *qseq = (uint8_t*)calloc(qlen + 33, 1);
    uint8_t *tseq = (uint8_t*)calloc(tlen + 33, 1);
    for (int i = 0; i < qlen; ++i)
        qseq[i] = seq_nt4_table[(uint8_t)query[i]];
    for (int i = 0; i < tlen; ++i)
        tseq[i] = seq_nt4_table[(uint8_t)ref[i]];

    ksw_extz2_sse(
        nullptr, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, ksw_matrix_m, ksw_matrix, parameters.gap_open, parameters.gap_extend, w, zdrop, parameters.end_bonus, flag, &ez
    );
    free(qseq);
    free(tseq);


    aln_info info;
    auto cigar = Cigar(ez.cigar, ez.n_cigar).to_eqx(query, ref);
    info.edit_distance = cigar.edit_distance();
    info.cigar = std::move(cigar);
    info.ref_start = 0;
    info.query_start = 0;
    if (ez.reach_end) {
        info.ref_end = ez.mqe_t + 1;
        info.query_end = query.size();
        info.sw_score = ez.mqe + parameters.end_bonus;
    } else {
        info.ref_end = ez.max_t + 1;
        info.query_end = ez.max_q + 1;
        info.sw_score = ez.max;
    }

    kfree(km, ez.cigar);
    return info;
}
