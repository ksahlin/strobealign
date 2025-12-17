#include "piecewisealigner.hpp"
#include <vector>

namespace Piecewise {

Aligner::Aligner(const AlignmentParameters& params, int k) 
    : params(params),
      k(k),
      gaps{.open = static_cast<int8_t>(-params.gap_open), .extend = static_cast<int8_t>(-params.gap_extend)},
      matrix(nullptr),
      x_drop_threshold(static_cast<int32_t>(params.x_drop_threshold))
{
    // create a matrix, that will be used for all blockaligner calls
    matrix = block_new_simple_aamatrix(static_cast<int8_t>(params.match), static_cast<int8_t>(-params.mismatch));
}

Aligner::~Aligner() {
    // freeing the matrix pointer
    if (matrix) {
        block_free_aamatrix(matrix);
        matrix = nullptr;
    }
}

// Converts a blockaligner Cigar to a strobealign Cigar
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

// Converts a blockaligner Cigar to a strobealign Cigar
// Also exchanged Insertions and Deletions
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

// Converts a blockaligner Cigar to a strobealign Cigar in the reversed order
// Also exchanged Insertions and Deletions
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

// returns the next power of 2
// more recent versions of c++ have built-in functions.
constexpr size_t next_power_of_2(size_t n) {
    if (n == 0) return 1;
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    if constexpr (sizeof(size_t) == 8) // if 64bits arch
        n |= n >> 32;
    return n + 1;
}

// returns a full global alignment using blockaligner
AlignmentResult Aligner::global_alignment(const std::string_view& query, const std::string_view& ref) const {
    AlignmentResult result;
    if (query.empty() || ref.empty()) return result;

    // blockaligner block ranges, set to the maximum values of the sequence for full global alignment
    // must be powers of 2
    SizeRange range;
    range.min = next_power_of_2(std::max<uintptr_t>(32, std::max(query.length(), ref.length())));
    range.max = next_power_of_2(std::max<uintptr_t>(128, std::max(query.length(), ref.length())));

    // creates the sequences for blockaligner
    PaddedBytes* q_padded = block_new_padded_aa(query.length(), range.max);
    PaddedBytes* r_padded = block_new_padded_aa(ref.length(), range.max);
    block_set_bytes_padded_aa(q_padded, reinterpret_cast<const uint8_t*>(query.data()), query.length(), range.max);
    block_set_bytes_padded_aa(r_padded, reinterpret_cast<const uint8_t*>(ref.data()), ref.length(), range.max);

    // makes a global alignment call with backtracking trace result
    BlockHandle block = block_new_aa_trace(query.length(), ref.length(), range.max);
    block_align_aa_trace(block, q_padded, r_padded, matrix, gaps, range, 0); // 0 is the x-drop thershold
    AlignResult res = block_res_aa_trace(block);

    // retrieve CIGAR
    Cigar* cigar_ptr = block_new_cigar(res.query_idx, res.reference_idx);
    block_cigar_eq_aa_trace(block, q_padded, r_padded, res.query_idx, res.reference_idx, cigar_ptr);
    size_t cigar_len = block_len_cigar(cigar_ptr);

    result.score = res.score;
    result.query_start = 0;
    result.query_end = res.query_idx;
    result.ref_start = 0;
    result.ref_end = res.reference_idx;
    result.cigar = build_cigar(cigar_ptr, cigar_len);

    // free memory pointers
    block_free_cigar(cigar_ptr);
    block_free_aa_trace(block);
    block_free_padded_aa(q_padded);
    block_free_padded_aa(r_padded);

    return result;
}

// returns a x-drop alignment using blockaligner
// we align the reference to the query profile to allow end_bonus scoring
AlignmentResult Aligner::xdrop_alignment(const std::string_view& query, const std::string_view& ref, bool reverse) const {
    AlignmentResult result;
    if (query.empty() || ref.empty()) return result;

    // blockaligner block ranges, set to 20% and 50% of the longest sequence
    // must be powers of 2
    SizeRange range;
    range.min = next_power_of_2(std::max<uintptr_t>(32, std::max(query.length(), ref.length()) / 5 ));
    range.max = next_power_of_2(std::max<uintptr_t>(128, std::max(query.length(), ref.length()) / 2 ));;

    std::string_view query_used = query;
    std::string_view ref_used = ref;
    std::string query_rev, ref_rev;

    // in case we x-drop in reverse (aligning before the first anchor),
    // we exchange the query for the reversed refernece and vice-versa
    if (reverse) {
        query_rev.assign(query.rbegin(), query.rend());
        ref_rev.assign(ref.rbegin(), ref.rend());
        query_used = query_rev;
        ref_used = ref_rev;
    }

    // creates the sequences for blockaligner
    PaddedBytes* q_padded = block_new_padded_aa(query_used.length(), range.max);
    PaddedBytes* r_padded = block_new_padded_aa(ref_used.length(), range.max);
    block_set_bytes_padded_aa(q_padded, reinterpret_cast<const uint8_t*>(query_used.data()), query_used.length(), range.max);
    block_set_bytes_padded_aa(r_padded, reinterpret_cast<const uint8_t*>(ref_used.data()), ref_used.length(), range.max);

    // creates a profile to align the reference on the query, so we can give the query ends a bonus score
    AAProfile* query_profile = block_new_aaprofile(query_used.length(), range.max, static_cast<int8_t>(-params.gap_extend));
    if (!query_profile) {
        block_free_padded_aa(q_padded);
        block_free_padded_aa(r_padded);
        return result;
    }

    // gives a score for each nucleotide combinaison
    for (size_t i = 1; i <= query_used.length(); i++) {
        for (int c : {'A', 'C', 'G', 'T', 'N'}) {
            if (c == query_used[i - 1]) {
                block_set_aaprofile(query_profile, i, c, static_cast<int8_t>(params.match));
            } else {
                block_set_aaprofile(query_profile, i, c, static_cast<int8_t>(-params.mismatch));
            }
        }
    }

    // set gap costs
    block_set_all_gap_open_C_aaprofile(query_profile, static_cast<int8_t>(-params.gap_open) - static_cast<int8_t>(-params.gap_extend));
    block_set_all_gap_close_C_aaprofile(query_profile, 0);
    block_set_all_gap_open_R_aaprofile(query_profile, static_cast<int8_t>(-params.gap_open) - static_cast<int8_t>(-params.gap_extend));

    // gives a bonus score to the end
    const size_t bonus_pos = query_used.length();
    for (int c : {'A', 'C', 'G', 'T', 'N'}) {
        int8_t current_score = block_get_aaprofile(query_profile, bonus_pos, c);
        block_set_aaprofile(query_profile, bonus_pos, c, current_score + params.end_bonus);
    }

    // makes a x-drop alignment call with backtracking trace result
    BlockHandle block = block_new_aa_trace_xdrop(ref_used.length(), query_used.length(), range.max);
    block_align_profile_aa_trace_xdrop(block, r_padded, query_profile, range, x_drop_threshold);
    AlignResult res = block_res_aa_trace_xdrop(block);

    // retrieve CIGAR
    Cigar* cigar_ptr = block_new_cigar(res.query_idx, res.reference_idx);
    block_cigar_eq_aa_trace_xdrop(block, r_padded, q_padded, res.query_idx, res.reference_idx, cigar_ptr);
    size_t cigar_len = block_len_cigar(cigar_ptr);

    result.score = res.score;

    if (reverse) {
        // returns the cigar reversed with insertions and deletions swapped
        result.cigar = build_cigar_reverse_swap_ID(cigar_ptr, cigar_len);
        result.query_start = query.length() - res.reference_idx;
        result.query_end = query.length();
        result.ref_start = ref.length() - res.query_idx;
        result.ref_end = ref.length();
    } else {
        // returns the cigar with insertions and deletions swapped
        result.cigar = build_cigar_swap_ID(cigar_ptr, cigar_len);
        result.query_start = 0;
        result.query_end = res.reference_idx;
        result.ref_start = 0;
        result.ref_end = res.query_idx;
    }

    // free memory pointers
    block_free_cigar(cigar_ptr);
    block_free_aaprofile(query_profile);
    block_free_aa_trace_xdrop(block);
    block_free_padded_aa(q_padded);
    block_free_padded_aa(r_padded);

    return result;
}

// Heuristic to remove the most amount of anchors flagged as spurious,
// a spurious anchor being likely not part of the best scoring alignment
void Aligner::remove_spurious_anchors(
    std::vector<Anchor>& anchors
) const {
    if (anchors.size() < 2) {
        return;
    }

    // First pruning: 
    // We remove all cluster of anchors creating 2 canceling indels within the tolerance range.
    // This will only happen inside the chain, so even if they were part of the best scoring path
    // they will be retrieved in the global alignment.
    
    const int diagonal_tolerance = 5;
    
    int tracked_indel = 0;
    size_t deviation_start = -1;
    
    for (size_t i = 1; i < anchors.size(); ++i) {
        const int query_diff = static_cast<int>(anchors[i].query_start) - static_cast<int>(anchors[i-1].query_start);
        const int ref_diff = static_cast<int>(anchors[i].ref_start) - static_cast<int>(anchors[i-1].ref_start);

        const int indel = query_diff - ref_diff;
       
        // 2 anchors are consider to deviate from the diagonal if
        // they create an indel bigger than the diagonal tolerance
        if (std::abs(indel) > diagonal_tolerance) {
            if (deviation_start == -1) {
                deviation_start = i;
                tracked_indel = indel;
            } else {
                if (std::abs(tracked_indel + indel) <= diagonal_tolerance) {
                    // if we find an indel canceling the tracked indel, we
                    // remove all anchors covered between the 2 indel positions
                    anchors.erase(anchors.begin() + deviation_start, anchors.begin() + i);
                    i -= i - deviation_start;
                    deviation_start = -1;
                    tracked_indel = 0;
                } else {
                    deviation_start = i;
                    tracked_indel = indel;
                }
            }
        }
    } 

    // Second pruning: 
    // We remove anchors of the ends of the chain if they create any indel of min_indel_size.
    // If they were part of the best scoring path, they should be retireved by the x-drop alignment.
    // the idea is that spurious anchors are much more difficult to detect on the ends of the chain,
    // so we look at a small ratio of anchors and remove any anchors creating deviations.

    const double edge_prune_ratio = 0.11; 
    const int min_indel_size = 0;

    const size_t max_prune_count = static_cast<size_t>(std::ceil(anchors.size() * edge_prune_ratio));

    // anchors at the beginning of the chain
    for (size_t i = 1; i < anchors.size() && i <= max_prune_count; ++i) {
        const int query_diff = static_cast<int>(anchors[i].query_start) - static_cast<int>(anchors[i - 1].query_start);
        const int ref_diff = static_cast<int>(anchors[i].ref_start) - static_cast<int>(anchors[i - 1].ref_start);
       
        const int indel = query_diff - ref_diff;

        if (std::abs(indel) > min_indel_size) {
            anchors.erase(anchors.begin(), anchors.begin() + i);
            break;
        }
    }

    const size_t max_prune_count_end = static_cast<size_t>(std::floor(anchors.size() * (1.0 - edge_prune_ratio)));
    const size_t last_idx = anchors.size() - 1;

    // anchors at the end of the chain
    for (size_t i = last_idx; i >= max_prune_count_end; --i) {
        const int query_diff = static_cast<int>(anchors[i].query_start) - static_cast<int>(anchors[i - 1].query_start);
        const int ref_diff = static_cast<int>(anchors[i].ref_start) - static_cast<int>(anchors[i - 1].ref_start);
        
        const int indel = query_diff - ref_diff;

        if (std::abs(indel) > min_indel_size) {
            anchors.erase(anchors.begin() + i, anchors.end());
            break;
        }
    }
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
    std::vector<Anchor>& anchors,
    const int padding
) const {
    AlignmentInfo result;

    remove_spurious_anchors(anchors);

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

        // Early return for piecewise global alignment heuristic:
        // Because of spurious strobemer matches on the end of the query, chains tend to have a cluster of spurious
        // anchors on the ends of the query, this heuristic stop trusting the anchors if they create a suspiciously
        // large deletion on the end of the query
        if (int(query.length()) - curr_start_query <= 200 && ref_diff - query_diff >= (int(query.length()) - prev_end_query)/2 - k) {
            align_after_last_anchor(reference, query, prev_anchor, padding + ref_diff - query_diff, &result);
            result.edit_distance = result.cigar.edit_distance();
            return result;
        }

        if (ref_diff > 0 && query_diff > 0){
            // No overlap between anchors, we do hamming and/or global alignment
            const std::string_view query_part(query.data() + prev_end_query, query_diff);
            const std::string_view ref_part(reference.data() + prev_end_ref, ref_diff);

            if (ref_diff == query_diff) {
                const AlignmentInfo hamming_aligned = hamming_align_global(query_part, ref_part, params.match, params.mismatch);

                // if we have 85% matches with hamming, we consider it the valid alignment
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
            // Overlap between anchors, no need to align
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
