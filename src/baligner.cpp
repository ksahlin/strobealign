#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include "baligner.hpp"

std::vector<OpLen> build_cigar_vector(const Cigar* cigar, size_t cigar_len) {
    std::vector<OpLen> cigar_vec;
    cigar_vec.reserve(cigar_len);
    for (size_t i = 0; i < cigar_len; i++) {
        cigar_vec.push_back(block_get_cigar(cigar, i));
    }
    return cigar_vec;
}

std::vector<OpLen> reverse_cigar_vector(const Cigar* cigar, size_t cigar_len) {
    std::vector<OpLen> reversed_cigar_vec;
    reversed_cigar_vec.reserve(cigar_len);
    for (int i = cigar_len - 1; i >= 0; i--) {
        reversed_cigar_vec.push_back(block_get_cigar(cigar, i));
    }
    return reversed_cigar_vec;
}

std::string reverse_string(const std::string& s) {
    std::string temp = s;
    std::reverse(temp.begin(), temp.end());
    return temp;
}

enum class AlignmentMode {
    Global,
    FreeQueryEnd,
    FreeQueryStart
};

AlignmentResult run_block_alignment(const std::string& query, const std::string& ref, AlignmentMode mode, const AlignmentScoring& scoring_params) {
    AlignmentResult result;

    if (query.length() == 0 || ref.length() == 0) {
        result.score = std::numeric_limits<int>::min();
        result.cigar.clear();
        result.query_start = 0; result.query_end = 0;
        result.ref_start = 0; result.ref_end = 0;
        return result;
    }

    std::string processed_query = query;
    std::string processed_ref = ref;
    size_t original_query_len = query.length();
    size_t original_ref_len = ref.length();

    if (mode == AlignmentMode::FreeQueryStart) {
        processed_query = reverse_string(query);
        processed_ref = reverse_string(ref);
    }

    SizeRange range = {.min = 32, .max = 256};
    Gaps gaps = {.open = scoring_params.gap_open, .extend = scoring_params.gap_extend};
    AAMatrix* dna_matrix = block_new_simple_aamatrix(scoring_params.match, scoring_params.mismatch);

    PaddedBytes* q_padded = block_new_padded_aa(processed_query.length(), range.max);
    PaddedBytes* r_padded = block_new_padded_aa(processed_ref.length(), range.max);

    block_set_bytes_padded_aa(q_padded, (const uint8_t*)processed_query.c_str(), processed_query.length(), range.max);
    block_set_bytes_padded_aa(r_padded, (const uint8_t*)processed_ref.c_str(), processed_ref.length(), range.max);

    BlockHandle block = nullptr;
    AlignResult res;
    Cigar* cigar_ptr = nullptr;

    const int32_t x_drop_threshold = 800; // make it a param at some point

    if (mode == AlignmentMode::Global) {
        block = block_new_aa_trace(original_query_len, original_ref_len, range.max);
        block_align_aa_trace(block, q_padded, r_padded, dna_matrix, gaps, range, x_drop_threshold);
        res = block_res_aa_trace(block);
        cigar_ptr = block_new_cigar(res.query_idx, res.reference_idx);
        block_cigar_eq_aa_trace(block, q_padded, r_padded, res.query_idx, res.reference_idx, cigar_ptr);
        block_free_aa_trace(block);
    } else {
        block = block_new_aa_trace_xdrop(processed_query.length(), processed_ref.length(), range.max);
        block_align_aa_trace_xdrop(block, q_padded, r_padded, dna_matrix, gaps, range, x_drop_threshold);
        res = block_res_aa_trace_xdrop(block);
        cigar_ptr = block_new_cigar(res.query_idx, res.reference_idx);
        block_cigar_eq_aa_trace_xdrop(block, q_padded, r_padded, res.query_idx, res.reference_idx, cigar_ptr);
        block_free_aa_trace_xdrop(block);
    }

    result.score = res.score;
    size_t cigar_len = block_len_cigar(cigar_ptr);

    if (mode == AlignmentMode::FreeQueryStart) {
        result.query_start = original_query_len - res.query_idx;
        result.query_end = original_query_len;
        result.ref_start = original_ref_len - res.reference_idx;
        result.ref_end = original_ref_len;
        result.cigar = reverse_cigar_vector(cigar_ptr, cigar_len);
    } else {
        result.query_start = 0;
        result.query_end = res.query_idx;
        result.ref_start = 0;
        result.ref_end = res.reference_idx;
        result.cigar = build_cigar_vector(cigar_ptr, cigar_len);
    }

    block_free_cigar(cigar_ptr);
    block_free_padded_aa(q_padded);
    block_free_padded_aa(r_padded);
    block_free_aamatrix(dna_matrix);

    return result;
}

AlignmentResult global_alignment(const std::string& query, const std::string& ref, const AlignmentScoring& scoring_params) {
    return run_block_alignment(query, ref, AlignmentMode::Global, scoring_params);
}

AlignmentResult free_query_end_alignment(const std::string& query, const std::string& ref, const AlignmentScoring& scoring_params) {
    return run_block_alignment(query, ref, AlignmentMode::FreeQueryEnd, scoring_params);
}

AlignmentResult free_query_start_alignment(const std::string& query, const std::string& ref, const AlignmentScoring& scoring_params) {
    return run_block_alignment(query, ref, AlignmentMode::FreeQueryStart, scoring_params);
}

