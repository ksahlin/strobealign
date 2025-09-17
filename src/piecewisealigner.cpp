#include "piecewisealigner.hpp"

namespace Piecewise {

Aligner::Aligner(const AlignmentParameters& params, int k) 
    : params(params),
      k(k),
      range{.min = static_cast<uintptr_t>(params.min_block), .max = static_cast<uintptr_t>(params.max_block)},
      gaps{.open = static_cast<int8_t>(-params.gap_open), .extend = static_cast<int8_t>(-params.gap_extend)},
      matrix(nullptr),
      x_drop_threshold(static_cast<int32_t>(params.x_drop_threshold))
{
    matrix = block_new_simple_aamatrix(static_cast<int8_t>(params.match), static_cast<int8_t>(-params.mismatch));
}

Aligner::~Aligner() {
    if (matrix) {
        block_free_aamatrix(matrix);
        matrix = nullptr;
    }
}

Cigar Aligner::build_cigar(const Cigar* cigar, size_t cigar_len) const {
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

Cigar Aligner::build_cigar_swap_ID(const Cigar* cigar, size_t cigar_len) const {
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

Cigar Aligner::build_cigar_reverse_swap_ID(const Cigar* cigar, size_t cigar_len) const {
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

AlignmentResult Aligner::global_alignment(const std::string_view& query, const std::string_view& ref) const {
    AlignmentResult result;

    if (query.empty() || ref.empty()) {
        return result;
    }

    PaddedBytes* q_padded = block_new_padded_aa(query.length(), range.max);
    PaddedBytes* r_padded = block_new_padded_aa(ref.length(), range.max);
    block_set_bytes_padded_aa(q_padded, reinterpret_cast<const uint8_t*>(query.data()), query.length(), range.max);
    block_set_bytes_padded_aa(r_padded, reinterpret_cast<const uint8_t*>(ref.data()), ref.length(), range.max);

    BlockHandle block = block_new_aa_trace(query.length(), ref.length(), range.max);
    block_align_aa_trace(block, q_padded, r_padded, matrix, gaps, range, 0);
    AlignResult res = block_res_aa_trace(block);
    
    Cigar* cigar_ptr = block_new_cigar(res.query_idx, res.reference_idx);
    block_cigar_eq_aa_trace(block, q_padded, r_padded, res.query_idx, res.reference_idx, cigar_ptr);
    size_t cigar_len = block_len_cigar(cigar_ptr);
    
    result.score = res.score;
    result.query_start = 0;
    result.query_end = res.query_idx;
    result.ref_start = 0;
    result.ref_end = res.reference_idx;
    result.cigar = build_cigar(cigar_ptr, cigar_len);
    
    block_free_cigar(cigar_ptr);
    block_free_aa_trace(block);
    block_free_padded_aa(q_padded);
    block_free_padded_aa(r_padded);

    return result;
}

AlignmentResult Aligner::xdrop_alignment(const std::string_view& query, const std::string_view& ref, bool reverse) const {
    AlignmentResult result;

    if (query.empty() || ref.empty()) {
        return result;
    }

    std::string_view query_used = query;
    std::string_view ref_used = ref;
    std::string query_rev, ref_rev;

    if (reverse) {
        query_rev.assign(query.rbegin(), query.rend());
        ref_rev.assign(ref.rbegin(), ref.rend());
        query_used = query_rev;
        ref_used = ref_rev;
    }

    PaddedBytes* q_padded = block_new_padded_aa(query_used.length(), range.max);
    PaddedBytes* r_padded = block_new_padded_aa(ref_used.length(), range.max);
    block_set_bytes_padded_aa(q_padded, reinterpret_cast<const uint8_t*>(query_used.data()), query_used.length(), range.max);
    block_set_bytes_padded_aa(r_padded, reinterpret_cast<const uint8_t*>(ref_used.data()), ref_used.length(), range.max);

    AAProfile* query_profile = block_new_aaprofile(query_used.length(), range.max, static_cast<int8_t>(-params.gap_extend));
    if (!query_profile) {
        block_free_padded_aa(q_padded);
        block_free_padded_aa(r_padded);
        return result;
    }
    
    for (size_t i = 1; i <= query_used.length(); i++) {
        for (int c = 'A'; c <= 'Z'; c++) {
            if (c == query_used[i - 1]) {
                block_set_aaprofile(query_profile, i, c, static_cast<int8_t>(params.match));
            } else {
                block_set_aaprofile(query_profile, i, c, static_cast<int8_t>(-params.mismatch));
            }
        }
    }

    block_set_all_gap_open_C_aaprofile(query_profile, static_cast<int8_t>(-params.gap_open) - static_cast<int8_t>(-params.gap_extend));
    block_set_all_gap_close_C_aaprofile(query_profile, 0);
    block_set_all_gap_open_R_aaprofile(query_profile, static_cast<int8_t>(-params.gap_open) - static_cast<int8_t>(-params.gap_extend));

    const size_t bonus_pos = query_used.length();
    for (int c = 'A'; c <= 'Z'; c++) {
        int8_t current_score = block_get_aaprofile(query_profile, bonus_pos, c);
        block_set_aaprofile(query_profile, bonus_pos, c, current_score + params.end_bonus);
    }

    BlockHandle block = block_new_aa_trace_xdrop(ref_used.length(), query_used.length(), range.max);
    if (!block) {
        block_free_aaprofile(query_profile);
        block_free_padded_aa(q_padded);
        block_free_padded_aa(r_padded);
        return result;
    }

    block_align_profile_aa_trace_xdrop(block, r_padded, query_profile, range, x_drop_threshold);
    AlignResult res = block_res_aa_trace_xdrop(block);
    
    Cigar* cigar_ptr = block_new_cigar(res.query_idx, res.reference_idx);
    block_cigar_eq_aa_trace_xdrop(block, r_padded, q_padded, res.query_idx, res.reference_idx, cigar_ptr);
    size_t cigar_len = block_len_cigar(cigar_ptr);
    
    result.score = res.score;

    if (reverse) {
        result.cigar = build_cigar_reverse_swap_ID(cigar_ptr, cigar_len);
        result.query_start = query.length() - res.reference_idx;
        result.query_end = query.length();
        result.ref_start = ref.length() - res.query_idx;
        result.ref_end = ref.length();
    } else {
        result.cigar = build_cigar_swap_ID(cigar_ptr, cigar_len);
        result.query_start = 0;
        result.query_end = res.reference_idx;
        result.ref_start = 0;
        result.ref_end = res.query_idx;
    }
    
    block_free_cigar(cigar_ptr);
    block_free_aaprofile(query_profile);
    block_free_aa_trace_xdrop(block);
    block_free_padded_aa(q_padded);
    block_free_padded_aa(r_padded);

    return result;
}


void Aligner::align_before_first_anchor(
    const std::string& reference,
    const std::string& query,
    const Anchor& first_anchor,
    const int padding,
    AlignmentInfo* result
) const {
    if (first_anchor.query_start > 0 && first_anchor.ref_start > 0) {
        const std::string_view query_part(query.data(), first_anchor.query_start);
        const size_t ref_start = std::max(0, static_cast<int>(first_anchor.ref_start) - (static_cast<int>(query_part.length()) + padding));
        const std::string_view ref_part(reference.data() + ref_start, first_anchor.ref_start - ref_start);

        const AlignmentResult pre_align = xdrop_alignment(query_part, ref_part, true);

        if (pre_align.score == 0) {
            result->query_start = first_anchor.query_start;
            result->ref_start = first_anchor.ref_start;
            result->cigar.push(CIGAR_SOFTCLIP, result->query_start);
        } else {
            result->sw_score += pre_align.score;
            result->query_start = pre_align.query_start;
            result->ref_start = ref_start + pre_align.ref_start;
            if (result->query_start > 0) {
                result->cigar.push(CIGAR_SOFTCLIP, result->query_start);
            }
            result->cigar += pre_align.cigar;
        }
    } else {
        result->query_start = first_anchor.query_start;
        result->ref_start = first_anchor.ref_start;
        if (result->query_start == 0) {
            result->sw_score += params.end_bonus;
        } else {
            result->cigar.push(CIGAR_SOFTCLIP, result->query_start);
        }
    }
}

void Aligner::align_after_last_anchor(
    const std::string& reference,
    const std::string& query,
    const Anchor& last_anchor,
    const int padding,
    AlignmentInfo* result
) const {
    const size_t last_anchor_end_query = last_anchor.query_start + k;
    const size_t last_anchor_end_ref = last_anchor.ref_start + k;
    
    if (last_anchor_end_query < query.length() && last_anchor_end_ref < reference.length()) {
        const std::string_view query_part(query.data() + last_anchor_end_query, query.length() - last_anchor_end_query);
        const size_t ref_part_end = std::min(reference.length(), last_anchor_end_ref + query_part.length() + padding);
        const std::string_view ref_part(reference.data() + last_anchor_end_ref, ref_part_end - last_anchor_end_ref);

        const AlignmentResult post_align = xdrop_alignment(query_part, ref_part, false);

        if (post_align.score == 0) {
            result->query_end = last_anchor_end_query;
            result->ref_end = last_anchor_end_ref;
        } else {
            result->sw_score += post_align.score;
            result->query_end = last_anchor_end_query + post_align.query_end;
            result->ref_end = last_anchor_end_ref + post_align.ref_end;
            result->cigar += post_align.cigar;
        }
    } else {
        result->query_end = last_anchor_end_query;
        result->ref_end = last_anchor_end_ref;
        if (result->query_end == query.length()) {
            result->sw_score += params.end_bonus;
        }
    }

    if (result->query_end < query.length()) {
        result->cigar.push(CIGAR_SOFTCLIP, query.length() - result->query_end);
    }
}

AlignmentInfo Aligner::piecewise_extension_alignment(
    const std::string& reference,
    const std::string& query,
    const std::vector<Anchor>& anchors,
    const int padding
) const {
    AlignmentInfo result;

    align_before_first_anchor(reference, query, anchors[0], padding, &result);

    result.sw_score += k * params.match;
    result.cigar.push(CIGAR_EQ, k);

    for (size_t i = 1; i < anchors.size(); ++i) {
        const Anchor& anchor = anchors[i];
        const Anchor& prev_anchor = anchors[i - 1];

        const int curr_start_query = anchor.query_start;
        const int curr_start_ref = anchor.ref_start;
        const int prev_end_query = prev_anchor.query_start + k;
        const int prev_end_ref = prev_anchor.ref_start + k;

        const int ref_diff = curr_start_ref - prev_end_ref;
        const int query_diff = curr_start_query - prev_end_query;

        // magic heuristic to prune off annoying anchors on the query end
        if (int(query.length()) - curr_start_query <= 200 && ref_diff - query_diff >= (int(query.length()) - prev_end_query)/2) {
            align_after_last_anchor(reference, query, prev_anchor, padding + ref_diff - query_diff, &result);
            result.edit_distance = result.cigar.edit_distance();
            return result;
        }

        if (ref_diff > 0 && query_diff > 0){
            const std::string_view query_part(query.data() + prev_end_query, query_diff);
            const std::string_view ref_part(reference.data() + prev_end_ref, ref_diff);

            if (ref_diff == query_diff) {
                const AlignmentInfo hamming_aligned = hamming_align_global(query_part, ref_part, params.match, params.mismatch);

                if (hamming_aligned.sw_score >= params.match * static_cast<float>(query_part.size()) * 0.85f) {
                    result.sw_score += hamming_aligned.sw_score;
                    result.cigar += hamming_aligned.cigar;

                    result.sw_score += k * params.match;
                    result.cigar.push(CIGAR_EQ, k);
                    continue;
                }
            }

            const AlignmentResult aligned = global_alignment(query_part, ref_part);

            result.sw_score += aligned.score;
            result.cigar += aligned.cigar;

            result.sw_score += k * params.match;
            result.cigar.push(CIGAR_EQ, k);
        } else {
             if (ref_diff < query_diff) {
                const size_t inserted_part = -ref_diff + query_diff;
                result.sw_score += -params.gap_open + (inserted_part - 1) * -params.gap_extend;
                result.cigar.push(CIGAR_INS, inserted_part);

                const size_t matching_part = k + ref_diff;
                result.sw_score += matching_part * params.match;
                result.cigar.push(CIGAR_EQ, matching_part);
            } else if (ref_diff > query_diff) {
                const size_t deleted_part = -query_diff + ref_diff;
                result.sw_score += -params.gap_open + (deleted_part - 1) * -params.gap_extend;
                result.cigar.push(CIGAR_DEL, deleted_part);

                const size_t matching_part = k + query_diff;
                result.sw_score += matching_part * params.match;
                result.cigar.push(CIGAR_EQ, matching_part);
            } else {
                const size_t matching_part = k + ref_diff;
                result.sw_score += matching_part * params.match;
                result.cigar.push(CIGAR_EQ, matching_part);
            }
        }
    }

    align_after_last_anchor(reference, query, anchors.back(), padding, &result);
    
    result.edit_distance = result.cigar.edit_distance();
    return result;
}

}
