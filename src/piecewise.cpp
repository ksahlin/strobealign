#include <algorithm>
#include <random>
#include <string>
#include <vector>
#include "nam.hpp"
#include "piecewisealigner.hpp"
#include "revcomp.hpp"
#include "piecewise.hpp"
#include "revcomp.hpp"
#include "aligner.hpp"
#include "sam.hpp"
#include "logger.hpp"

static Logger& logger = Logger::get();

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
        Piecewise::Aligner aln(scoring_params, k);
        info = aln.piecewise_extension_alignment(ref, query, chain.anchors, read.size()/10);
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

    int considered = 0; 
    for (auto &chain : chains) {
        float score_dropoff = (float) chain.score / n_max.score;
        if (tries >= max_tries || (tries > 1 && best_edit_distance == 0) || score_dropoff < dropoff_threshold) {
            break;
        }
        considered++;

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

    if (logger.level() <= LOG_TRACE){
        logger.trace() << "Cigars:[";
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


            logger.trace() << "(" <<alignment.cigar << ",was_considered:"<< (curr < considered) <<  ",rstart=" << alignment.ref_start << ",SSW:" << alignment_ssw.cigar << ",SSW_rstart=" << alignment_ssw.ref_start << ")";
            curr++;
        }
        logger.trace() << "]\n";
    }
    
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
