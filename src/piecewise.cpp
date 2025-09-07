#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <vector>
#include "baligner.hpp"
#include "nam.hpp"
#include "revcomp.hpp"
#include "piecewise.hpp"
#include "revcomp.hpp"
#include "aligner.hpp"
#include "sam.hpp"

std::pair<Cigar, int> merge_cigar_elements_with_edit_distance(
    const std::vector<OpLen>& elements,
    uint query_start,
    uint query_end,
    int query_length
) {
    Cigar cigar;
    int edit_dist = 0;

    if (query_start > 0) {
        cigar.push(CIGAR_SOFTCLIP, query_start);
    }

    if (!elements.empty()) {
        auto operation_to_cigar = [](uint8_t op) -> uint8_t {
            switch (op) {
                case Operation::M:  return CIGAR_MATCH;
                case Operation::Eq: return CIGAR_EQ;
                case Operation::X:  return CIGAR_X;
                case Operation::I:  return CIGAR_INS;
                case Operation::D:  return CIGAR_DEL;
                default:
                    assert(false && "Unknown operation");
                    return CIGAR_MATCH;
            }
        };

        uint8_t last_op = elements[0].op;
        uintptr_t last_len = elements[0].len;

        for (size_t i = 1; i < elements.size(); ++i) {
            if (elements[i].op == last_op) {
                last_len += elements[i].len;
                continue;
            } 
            cigar.push(operation_to_cigar(last_op), static_cast<int>(last_len));

            if (last_op == Operation::X || last_op == Operation::D || last_op == Operation::I) {
                edit_dist += last_len;
            }

            last_op = elements[i].op;
            last_len = elements[i].len;
        }

        cigar.push(operation_to_cigar(last_op), static_cast<int>(last_len));
        if (last_op == Operation::X || last_op == Operation::D || last_op == Operation::I) {
            edit_dist += last_len;
        }
    }

    if (query_end < query_length) {
        cigar.push(CIGAR_SOFTCLIP, query_length - query_end);
    }

    return {cigar, edit_dist};
}

AlignmentInfo piecewise_extension_alignment(
    const std::string& reference,
    const std::string& query,
    const std::vector<Anchor>& anchors,
    const int k,
    const int padding,
    const AlignmentParameters& params
) {
    AlignmentInfo result;
    result.sw_score = 0;
    std::vector<OpLen> temp_cigar_elements;
    const Anchor& first_anchor = anchors[0];
    if (first_anchor.query_start > 0 && first_anchor.ref_start > 0) {
        const std::string query_part = query.substr(0, first_anchor.query_start);
        const size_t ref_start = std::max(0, static_cast<int>(first_anchor.ref_start) - (static_cast<int>(query_part.length()) + padding));
        const std::string ref_part = reference.substr(ref_start, first_anchor.ref_start - ref_start);
        AlignmentResult pre_align = xdrop_query_start_alignment(query_part, ref_part, params);

        if (pre_align.score == 0) {
            result.query_start = first_anchor.query_start;
            result.ref_start = first_anchor.ref_start;
        } else {
            result.sw_score += pre_align.score;
            result.query_start = pre_align.query_start;
            result.ref_start = ref_start + pre_align.ref_start;
            temp_cigar_elements.insert(temp_cigar_elements.end(), pre_align.cigar.begin(), pre_align.cigar.end());
        }
    } else {
        result.query_start = first_anchor.query_start;
        result.ref_start = first_anchor.ref_start;
        if (result.query_start == 0) {
            result.sw_score += params.end_bonus;
        }
    }

    result.sw_score += k * params.match;
    temp_cigar_elements.push_back({Operation::Eq, (uintptr_t)k});

    for (size_t i = 1; i < anchors.size(); ++i) {
        const Anchor& anchor = anchors[i];
        const Anchor& prev_anchor = anchors[i - 1];

        const int curr_start_query = anchor.query_start;
        const int curr_start_ref = anchor.ref_start;
        const int prev_end_query = prev_anchor.query_start + k;
        const int prev_end_ref = prev_anchor.ref_start + k;

        const int ref_diff = curr_start_ref - prev_end_ref;
        const int query_diff = curr_start_query - prev_end_query;


        if (ref_diff > 0 && query_diff > 0){
            const std::string query_part = query.substr(prev_end_query, query_diff);
            const std::string ref_part = reference.substr(prev_end_ref, ref_diff);

            AlignmentResult aligned = global_alignment(query_part, ref_part, params);
            result.sw_score += aligned.score;
            temp_cigar_elements.insert(temp_cigar_elements.end(), aligned.cigar.begin(), aligned.cigar.end());

            result.sw_score += k * params.match;
            temp_cigar_elements.push_back({Operation::Eq, (uintptr_t)k});
        } else {
             if (ref_diff < query_diff) {
                const size_t inserted_part = -ref_diff + query_diff;
                result.sw_score += -params.gap_open + (inserted_part - 1) * -params.gap_extend;
                temp_cigar_elements.push_back({Operation::I, (uintptr_t)inserted_part});

                const size_t matching_part = k + ref_diff;
                result.sw_score += matching_part * params.match;
                temp_cigar_elements.push_back({Operation::Eq, (uintptr_t)matching_part});
            } else if (ref_diff > query_diff) {
                const size_t deleted_part = -query_diff + ref_diff;
                result.sw_score += -params.gap_open + (deleted_part - 1) * -params.gap_extend;
                temp_cigar_elements.push_back({Operation::D, (uintptr_t)deleted_part});

                const size_t matching_part = k + query_diff;
                result.sw_score += matching_part * params.match;
                temp_cigar_elements.push_back({Operation::Eq, (uintptr_t)matching_part});
            } else {
                const size_t matching_part = k + ref_diff;
                result.sw_score += matching_part * params.match;
                temp_cigar_elements.push_back({Operation::Eq, (uintptr_t)matching_part});
            }
        }
    }

    const Anchor& last_anchor = anchors.back();
    const size_t last_anchor_end_query = last_anchor.query_start + k;
    const size_t last_anchor_end_ref = last_anchor.ref_start + k;
    if (last_anchor_end_query < query.length() && last_anchor_end_ref < reference.length()) {
        const std::string query_part = query.substr(last_anchor_end_query);
        const size_t ref_part_end = std::min(reference.length(), last_anchor_end_ref + query_part.length() + padding);
        const std::string ref_part = reference.substr(last_anchor_end_ref, ref_part_end - last_anchor_end_ref);

        AlignmentResult post_align = xdrop_query_end_alignment(query_part, ref_part, params);

        if (post_align.score == 0) {
            result.query_end = last_anchor_end_query;
            result.ref_end = last_anchor_end_ref;
        } else {
            result.sw_score += post_align.score;
            result.query_end = last_anchor_end_query + post_align.query_end;
            result.ref_end = last_anchor_end_ref + post_align.ref_end;
            temp_cigar_elements.insert(temp_cigar_elements.end(), post_align.cigar.begin(), post_align.cigar.end());
        }
    } else {
        result.query_end = last_anchor_end_query;
        result.ref_end = last_anchor_end_ref;
        if (result.query_end == uint(query.length())) {
            result.sw_score += params.end_bonus;
        }
    }

    const auto cigar_ed = merge_cigar_elements_with_edit_distance(
        temp_cigar_elements,
        result.query_start,
        result.query_end,
        static_cast<int>(query.length())
    );

    result.cigar = cigar_ed.first;
    result.edit_distance = cigar_ed.second;

    return result;
}

inline Alignment extend_seed_piecewise(
    const AlignmentParameters& scoring_params,
    const Chain& chain,
    const int k,
    const References& references,
    const Read& read
) {
    const std::string query = chain.is_revcomp ? read.rc : read.seq;
    const std::string& ref = references.sequences[chain.ref_id];

    const auto projected_ref_start = std::max(0, int(chain.ref_start) - int(chain.query_start));
    const auto projected_ref_end = std::min(chain.ref_end + query.size() - chain.query_end, ref.size());

    AlignmentInfo info;
    int result_ref_start;
    bool gapped = true;
    if (projected_ref_end - projected_ref_start == query.size()) {
        std::string ref_segm_ham = ref.substr(projected_ref_start, query.size());
        auto hamming_dist = hamming_distance(query, ref_segm_ham);

        if (hamming_dist >= 0 && (((float) hamming_dist / query.size()) < 0.05) ) { //Hamming distance worked fine, no need to ksw align
            info = hamming_align(query, ref_segm_ham, scoring_params.match, scoring_params.mismatch, scoring_params.end_bonus);
            result_ref_start = projected_ref_start + info.ref_start;
            gapped = false;
        }
    }
    if (gapped) {
        info = piecewise_extension_alignment(ref, query, chain.anchors, k, read.size()/10, scoring_params);
        result_ref_start = info.ref_start;
    }
    int softclipped = info.query_start + (query.size() - info.query_end);
    Alignment alignment;
    alignment.cigar = std::move(info.cigar);
    alignment.edit_distance = info.edit_distance;
    alignment.global_ed = info.edit_distance + softclipped;
    alignment.score = info.sw_score;
    alignment.ref_start = result_ref_start;
    alignment.length = info.ref_span();
    alignment.is_revcomp = chain.is_revcomp;
    alignment.is_unaligned = false;
    alignment.ref_id = chain.ref_id;
    alignment.gapped = gapped;

    return alignment;
}

using namespace klibpp;


#ifdef TRACE
bool reverse_nam_if_needed(Nam& nam, const Read& read, const References& references, int k) {
    auto read_len = read.size();
    std::string ref_start_kmer = references.sequences[nam.ref_id].substr(nam.ref_start, k);
    std::string ref_end_kmer = references.sequences[nam.ref_id].substr(nam.ref_end-k, k);

    std::string seq, seq_rc;
    if (nam.is_revcomp) {
        seq = read.rc;
        seq_rc = read.seq;
    } else {
        seq = read.seq;
        seq_rc = read.rc;
    }
    std::string read_start_kmer = seq.substr(nam.query_start, k);
    std::string read_end_kmer = seq.substr(nam.query_end-k, k);
    if (ref_start_kmer == read_start_kmer && ref_end_kmer == read_end_kmer) {
        return true;
    }

    // False forward or false reverse (possible due to symmetrical hash values)
    //    we need two extra checks for this - hopefully this will remove all the false hits we see (true hash collisions should be very few)
    int q_start_tmp = read_len - nam.query_end;
    int q_end_tmp = read_len - nam.query_start;
    // false reverse hit, change coordinates in nam to forward
    read_start_kmer = seq_rc.substr(q_start_tmp, k);
    read_end_kmer = seq_rc.substr(q_end_tmp - k, k);
    if (ref_start_kmer == read_start_kmer && ref_end_kmer == read_end_kmer) {
        nam.is_revcomp = !nam.is_revcomp;
        nam.query_start = q_start_tmp;
        nam.query_end = q_end_tmp;
        return true;
    }
    return false;
}

/*
 Extend a NAM so that it covers the entire read and return the resulting
 alignment.
*/
inline Alignment extend_seed(
    const Aligner& aligner,
    const Nam &nam,
    const References& references,
    const Read& read,
    bool consistent_nam
) {
    const std::string query = nam.is_revcomp ? read.rc : read.seq;
    const std::string& ref = references.sequences[nam.ref_id];

    const auto projected_ref_start = nam.projected_ref_start();
    const auto projected_ref_end = std::min(nam.ref_end + query.size() - nam.query_end, ref.size());

    AlignmentInfo info;
    int result_ref_start;
    bool gapped = true;
    if (projected_ref_end - projected_ref_start == query.size() && consistent_nam) {
        std::string ref_segm_ham = ref.substr(projected_ref_start, query.size());
        auto hamming_dist = hamming_distance(query, ref_segm_ham);

        if (hamming_dist >= 0 && (((float) hamming_dist / query.size()) < 0.05) ) { //Hamming distance worked fine, no need to ksw align
            info = hamming_align(query, ref_segm_ham, aligner.parameters.match, aligner.parameters.mismatch, aligner.parameters.end_bonus);
            result_ref_start = projected_ref_start + info.ref_start;
            gapped = false;
        }
    }
    if (gapped) {
        const int diff = std::abs(nam.ref_span() - nam.query_span());
        const int ext_left = std::min(50, projected_ref_start);
        const int ref_start = projected_ref_start - ext_left;
        const int ext_right = std::min(std::size_t(50), ref.size() - nam.ref_end);
        const auto ref_segm_size = read.size() + diff + ext_left + ext_right;
        const auto ref_segm = ref.substr(ref_start, ref_segm_size);
        auto opt_info = aligner.align(query, ref_segm);
        if (opt_info) {
            info = opt_info.value();
            result_ref_start = ref_start + info.ref_start;
        } else {
            // TODO This function should instead return an std::optional<Alignment>
            Alignment alignment;
            alignment.is_unaligned = true;
            alignment.edit_distance = 100000;
            alignment.ref_start = 0;
            alignment.score = -100000;

            return alignment;
        }
    }
    int softclipped = info.query_start + (query.size() - info.query_end);
    Alignment alignment;
    alignment.cigar = std::move(info.cigar);
    alignment.edit_distance = info.edit_distance;
    alignment.global_ed = info.edit_distance + softclipped;
    alignment.score = info.sw_score;
    alignment.ref_start = result_ref_start;
    alignment.length = info.ref_span();
    alignment.is_revcomp = nam.is_revcomp;
    alignment.is_unaligned = false;
    alignment.ref_id = nam.ref_id;
    alignment.gapped = gapped;

    return alignment;
}
#endif


void align_single_piecewse(
    const Aligner &aligner, //just to compare with SSW
    const AlignmentParameters& scoring_params,
    Sam& sam,
    const std::vector<Chain>& chains,
    const KSeq& record,
    const int k,
    const References& references,
    Details& details,
    float dropoff_threshold,
    int max_tries,
    unsigned max_secondary,
    std::minstd_rand& random_engine
) {
    if (chains.empty()) {
        sam.add_unmapped(record);
        return;
    }

    Read read(record.seq);
    std::vector<Alignment> alignments;
    int tries = 0;
    Chain n_max = chains[0];

    int best_edit_distance = std::numeric_limits<int>::max();
    int best_score = 0;
    int second_best_score = 0;
    int alignments_with_best_score = 0;
    size_t best_index = 0;

    Alignment best_alignment;
    best_alignment.is_unaligned = true;

#ifdef TRACE
    int considered = 0; 
#endif
    for (auto &chain : chains) {
        float score_dropoff = (float) chain.score / n_max.score;
        if (tries >= max_tries || (tries > 1 && best_edit_distance == 0) || score_dropoff < dropoff_threshold) {
            break;
        }
#ifdef TRACE
    considered++; 
#endif
        // bool consistent_nam = reverse_nam_if_needed(nam, read, references, k);
        // details.inconsistent_nams += !consistent_nam;
        auto alignment = extend_seed_piecewise(scoring_params, chain, k, references, read);
        details.tried_alignment++;
        if (alignment.is_unaligned) {
            tries++;
            continue;
        }
        details.gapped += alignment.gapped;

        if (max_secondary > 0) {
            alignments.emplace_back(alignment);
        }

        if (alignment.score >= best_score) {
            second_best_score = best_score;
            bool update_best = false;
            if (alignment.score > best_score) {
                alignments_with_best_score = 1;
                update_best = true;
            } else {
                assert(alignment.score == best_score);
                // Two or more alignments have the same best score - count them
                alignments_with_best_score++;

                // Pick one randomly using reservoir sampling
                std::uniform_int_distribution<> distrib(1, alignments_with_best_score);
                if (distrib(random_engine) == 1) {
                    update_best = true;
                }
            }
            if (update_best) {
                best_score = alignment.score;
                best_alignment = std::move(alignment);
                best_index = tries;
                if (max_secondary == 0) {
                    best_edit_distance = best_alignment.global_ed;
                }
            }
        } else if (alignment.score > second_best_score) {
            second_best_score = alignment.score;
        }
        tries++;
    }
    if (best_alignment.is_unaligned) {
        sam.add_unmapped(record);
        return;
    }
    details.best_alignments = alignments_with_best_score;
    uint8_t mapq = (60.0 * (best_score - second_best_score) + best_score - 1) / best_score;
    bool is_primary = true;
    sam.add(best_alignment, record, read.rc, mapq, is_primary, details);
#ifdef TRACE
    std::cerr << "Cigars:[";
    int curr = 0;
    for (const auto& chain : chains) {
        auto alignment = extend_seed_piecewise(scoring_params, chain, k, references, read);
        
        Nam nam = Nam {
            .nam_id = -1,
            .query_start = (int)chain.query_start,
            .query_end = (int)chain.query_end,
            .query_prev_match_startpos = -1,
            .ref_start = (int)chain.ref_start,
            .ref_end = (int)chain.ref_end,
            .ref_prev_match_startpos = -1,
            .n_matches = (int)chain.anchors.size(),
            .ref_id = (int)chain.ref_id,
            .score = chain.score,
            .is_revcomp = chain.is_revcomp
        };
        bool consistent_nam = reverse_nam_if_needed(nam, read, references, k);
        auto alignment_ssw = extend_seed(aligner, nam, references, read, consistent_nam);


        std::cerr << "(" <<alignment.cigar << ",was_considered:"<< (curr < considered) <<  ",rstart=" << alignment.ref_start << ",SSW:" << alignment_ssw.cigar << ",SSW_rstart=" << alignment_ssw.ref_start << ")";
        curr++;
    }
    std::cerr << "]\n";
#endif
    if (max_secondary == 0) {
        return;
    }

    // Secondary alignments

    // Remove the alignment that was already output
    if (alignments.size() > 1) {
        std::swap(alignments[best_index], alignments[alignments.size() - 1]);
    }
    alignments.resize(alignments.size() - 1);

    // Sort remaining alignments by score, highest first
    std::sort(alignments.begin(), alignments.end(),
        [](const Alignment& a, const Alignment& b) -> bool {
            return a.score > b.score;
        }
    );

    // Output secondary alignments
    size_t n = 0;
    for (const auto& alignment : alignments) {
        if (
            n >= max_secondary
            || alignment.score - best_score > 2*scoring_params.mismatch + scoring_params.gap_open
        ) {
            break;
        }
        bool is_primary = false;
        sam.add(alignment, record, read.rc, mapq, is_primary, details);
        n++;
    }
}
