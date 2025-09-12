#include <cstdint>
#include <string>
#include "block-aligner/c/block_aligner.h"
#include "baligner.hpp"

Cigar build_cigar(const Cigar* cigar, size_t cigar_len) {
    Cigar result;
    for (size_t i = 0; i < cigar_len; i++) {
        OpLen oplen = block_get_cigar(cigar, i);
        uint8_t op_code;
        switch (oplen.op) {
            case Operation::M:
                op_code = CIGAR_MATCH;
                break;
            case Operation::Eq:
                op_code = CIGAR_EQ;
                break;
            case Operation::X:
                op_code = CIGAR_X;
                break;
            case Operation::I:
                op_code = CIGAR_INS;
                break;
            case Operation::D:
                op_code = CIGAR_DEL;
                break;
            default:
                continue;
        }
        result.push(op_code, oplen.len);
    }
    return result;
}

Cigar build_cigar_sawp_ID(const Cigar* cigar, size_t cigar_len) {
    Cigar result;
    for (size_t i = 0; i < cigar_len; i++) {
        OpLen oplen = block_get_cigar(cigar, i);
        uint8_t op_code;
        switch (oplen.op) {
            case Operation::M:
                op_code = CIGAR_MATCH;
                break;
            case Operation::Eq:
                op_code = CIGAR_EQ;
                break;
            case Operation::X:
                op_code = CIGAR_X;
                break;
            case Operation::I:
                op_code = CIGAR_DEL;
                break;
            case Operation::D:
                op_code = CIGAR_INS;
                break;
            default:
                continue;
        }
        result.push(op_code, oplen.len);
    }
    return result;
}

Cigar build_cigar_reverse_sawp_ID(const Cigar* cigar, size_t cigar_len) {
    Cigar result;
    for (size_t i = cigar_len; i > 0; --i) {
        OpLen oplen = block_get_cigar(cigar, i - 1);
        uint8_t op_code;
        switch (oplen.op) {
            case Operation::M:
                op_code = CIGAR_MATCH;
                break;
            case Operation::Eq:
                op_code = CIGAR_EQ;
                break;
            case Operation::X:
                op_code = CIGAR_X;
                break;
            case Operation::I:
                op_code = CIGAR_DEL;
                break;
            case Operation::D:
                op_code = CIGAR_INS;
                break;
            default:
                continue;
        }
        result.push(op_code, oplen.len);
    }
    return result;
}

AlignmentResult global_alignment(const std::string_view& query, const std::string_view& ref, const AlignmentParameters& params) {
    AlignmentResult result;

    if (query.empty() || ref.empty()) {
        return result;
    }

    SizeRange range = {.min = uintptr_t(params.min_block), .max = uintptr_t(params.max_block)};
    Gaps gaps = {.open = int8_t(-params.gap_open), .extend = int8_t(-params.gap_extend)};

    PaddedBytes* q_padded = block_new_padded_aa(query.length(), range.max);
    PaddedBytes* r_padded = block_new_padded_aa(ref.length(), range.max);
    block_set_bytes_padded_aa(q_padded, reinterpret_cast<const uint8_t*>(query.data()), query.length(), range.max);
    block_set_bytes_padded_aa(r_padded, reinterpret_cast<const uint8_t*>(ref.data()), ref.length(), range.max);

    BlockHandle block = nullptr;
    AlignResult res;
    Cigar* cigar_ptr = nullptr;

    AAMatrix* dna_matrix = block_new_simple_aamatrix(int8_t(params.match), int8_t(-params.mismatch));
    block = block_new_aa_trace(query.length(), ref.length(), range.max);
    block_align_aa_trace(block, q_padded, r_padded, dna_matrix, gaps, range, 0);
    res = block_res_aa_trace(block);
    cigar_ptr = block_new_cigar(res.query_idx, res.reference_idx);
    block_cigar_eq_aa_trace(block, q_padded, r_padded, res.query_idx, res.reference_idx, cigar_ptr);
    size_t cigar_len = block_len_cigar(cigar_ptr);
    
    block_free_aamatrix(dna_matrix);
    block_free_aa_trace(block);

    result.score = res.score;
    result.query_start = 0;
    result.query_end = res.query_idx;
    result.ref_start = 0;
    result.ref_end = res.reference_idx;
    result.cigar = build_cigar(cigar_ptr, cigar_len);
    
    block_free_cigar(cigar_ptr);
    block_free_padded_aa(q_padded);
    block_free_padded_aa(r_padded);

    return result;
}

AlignmentResult xdrop_alignment(const std::string_view& query, const std::string_view& ref, const AlignmentParameters& params, bool reverse) {
    AlignmentResult result;

    if (query.empty() || ref.empty()) {
        return result;
    }

    std::string_view query_used = query;
    std::string_view ref_used = ref;
    std::string query_rev;
    std::string ref_rev;

    if (reverse) {
        query_rev.assign(query.rbegin(), query.rend());
        ref_rev.assign(ref.rbegin(), ref.rend());
        query_used = query_rev;
        ref_used = ref_rev;
    }

    SizeRange range = {.min = uintptr_t(params.min_block), .max = uintptr_t(params.max_block)};
    const int32_t x_drop_threshold = int32_t(params.x_drop_threshold);

    PaddedBytes* q_padded = block_new_padded_aa(query_used.length(), range.max);
    PaddedBytes* r_padded = block_new_padded_aa(ref_used.length(), range.max);
    block_set_bytes_padded_aa(q_padded, reinterpret_cast<const uint8_t*>(query_used.data()), query_used.length(), range.max);
    block_set_bytes_padded_aa(r_padded, reinterpret_cast<const uint8_t*>(ref_used.data()), ref_used.length(), range.max);

    BlockHandle block = nullptr;
    AlignResult res;
    Cigar* cigar_ptr = nullptr;

    AAProfile* query_profile = block_new_aaprofile(query_used.length(), range.max, int8_t(-params.gap_extend));
    
    for (size_t i = 1; i <= query_used.length(); i++) {
        for (int c = 'A'; c <= 'Z'; c++) {
            if (c == query_used[i - 1]) {
                block_set_aaprofile(query_profile, i, c, int8_t(params.match));
            } else {
                block_set_aaprofile(query_profile, i, c, int8_t(-params.mismatch));
            }
        }
    }

    block_set_all_gap_open_C_aaprofile(query_profile, int8_t(-params.gap_open) - int8_t(-params.gap_extend));
    block_set_all_gap_close_C_aaprofile(query_profile, 0);
    block_set_all_gap_open_R_aaprofile(query_profile, int8_t(-params.gap_open) - int8_t(-params.gap_extend));

    const size_t bonus_pos = query_used.length();
    for (int c = 'A'; c <= 'Z'; c++) {
        int8_t current_score = block_get_aaprofile(query_profile, bonus_pos, c);
        block_set_aaprofile(query_profile, bonus_pos, c, current_score + params.end_bonus);
    }

    block = block_new_aa_trace_xdrop(ref_used.length(), query_used.length(), range.max);
    block_align_profile_aa_trace_xdrop(block, r_padded, query_profile, range, x_drop_threshold);
    res = block_res_aa_trace_xdrop(block);
    cigar_ptr = block_new_cigar(res.query_idx, res.reference_idx);
    block_cigar_eq_aa_trace_xdrop(block, r_padded, q_padded, res.query_idx, res.reference_idx, cigar_ptr);
    size_t cigar_len = block_len_cigar(cigar_ptr);
   
    block_free_aaprofile(query_profile);
    block_free_aa_trace_xdrop(block);

    result.score = res.score;

    if (reverse) {
        result.cigar = build_cigar_reverse_sawp_ID(cigar_ptr, cigar_len);
        result.query_start = query.length() - res.reference_idx;
        result.query_end = query.length();
        result.ref_start = ref.length() - res.query_idx;
        result.ref_end = ref.length();
    } else {
        result.cigar = build_cigar_sawp_ID(cigar_ptr, cigar_len);
        result.query_start = 0;
        result.query_end = res.reference_idx;
        result.ref_start = 0;
        result.ref_end = res.query_idx;
    }

    block_free_cigar(cigar_ptr);
    block_free_padded_aa(q_padded);
    block_free_padded_aa(r_padded);

    return result;
}
