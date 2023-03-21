#include "aln.hpp"

#include <iostream>
#include <math.h>
#include <sstream>
#include "revcomp.hpp"
#include "timer.hpp"
#include "nam.hpp"
#include "paf.hpp"
#include "aligner.hpp"

using namespace klibpp;

static inline alignment get_alignment(
    const Aligner& aligner,
    const Nam &n,
    const References& references,
    const Read& read,
    bool fits
);

static inline bool score(const Nam &a, const Nam &b) {
    return a.score > b.score;
}

/*
 * Determine whether the NAM represents a match to the forward or
 * reverse-complemented sequence by checking in which orientation the
 * first and last strobe in the NAM match
 *
 * - If first and last strobe match in forward orientation, return true.
 * - If first and last strobe match in reverse orientation, update the NAM
 *   in place and return true.
 * - If first and last strobe do not match consistently, return false.
 */
bool reverse_nam_if_needed(Nam& n, const Read& read, const References& references, int k) {
    auto read_len = read.size();
    std::string ref_start_kmer = references.sequences[n.ref_id].substr(n.ref_s, k);
    std::string ref_end_kmer = references.sequences[n.ref_id].substr(n.ref_e-k, k);

    std::string seq, seq_rc;
    if (n.is_rc) {
        seq = read.rc;
        seq_rc = read.seq;
    } else {
        seq = read.seq;
        seq_rc = read.rc;
    }
    std::string read_start_kmer = seq.substr(n.query_s, k);
    std::string read_end_kmer = seq.substr(n.query_e-k, k);
    if (ref_start_kmer == read_start_kmer && ref_end_kmer == read_end_kmer) {
        return true;
    }

    // False forward or false reverse (possible due to symmetrical hash values)
    //    we need two extra checks for this - hopefully this will remove all the false hits we see (true hash collisions should be very few)
    int q_start_tmp = read_len - n.query_e;
    int q_end_tmp = read_len - n.query_s;
    // false reverse hit, change coordinates in nam to forward
    read_start_kmer = seq_rc.substr(q_start_tmp, k);
    read_end_kmer = seq_rc.substr(q_end_tmp - k, k);
    if (ref_start_kmer == read_start_kmer && ref_end_kmer == read_end_kmer) {
        n.is_rc = !n.is_rc;
        n.query_s = q_start_tmp;
        n.query_e = q_end_tmp;
        return true;
    }
    return false;
}

static inline void align_SE_secondary_hits(
    const Aligner& aligner,
    Sam& sam,
    std::vector<Nam>& all_nams,
    const KSeq& record,
    int k,
    const References& references,
    AlignmentStatistics &statistics,
    float dropoff,
    int max_tries,
    unsigned max_secondary
) {
    Read read(record.seq);

    if (all_nams.empty()) {
        sam.add_unmapped(record);
        return;
    }

    std::vector<alignment> alignments;
    int cnt = 0;
    float score_dropoff;
    Nam n_max = all_nams[0];

    int best_align_dist = ~0U >> 1;
    int best_align_sw_score = -1000;

    alignment best_sam_aln;
    best_sam_aln.sw_score = -100000;
    best_sam_aln.is_unaligned = true;
    int min_mapq_diff = best_align_dist;
    for (auto &n : all_nams) {
        score_dropoff = (float) n.n_hits / n_max.n_hits;

        if ( (cnt >= max_tries) || best_align_dist == 0 || score_dropoff < dropoff){ // only consider top 20 hits as minimap2 and break if alignment is exact match to reference or the match below droppoff cutoff.
            break;
        }
        statistics.tot_all_tried ++;
        bool fits = reverse_nam_if_needed(n, read, references, k);
        if (!fits){
            statistics.did_not_fit++;
        }
        auto sam_aln = get_alignment(aligner, n, references, read, fits);

        int diff_to_best = std::abs(best_align_sw_score - sam_aln.sw_score);
        min_mapq_diff = std::min(min_mapq_diff, diff_to_best);

        if (sam_aln.sw_score > best_align_sw_score) {
            min_mapq_diff = std::max(0, sam_aln.sw_score - best_align_sw_score); // new distance to next best match
            best_align_sw_score = sam_aln.sw_score;
            best_sam_aln = sam_aln;
            if (max_secondary == 1) {
                best_align_dist = sam_aln.global_ed;
            }
        }
//        sam_aln.ed = 10000; // init
        if (max_secondary > 1) {
            alignments.push_back(sam_aln);
        }
        cnt++;
    }

    if (max_secondary == 1) {
        best_sam_aln.mapq = std::min(min_mapq_diff, 60);
        sam.add(best_sam_aln, record, read.rc);
        return;
    }
    // Sort alignments by score, highest first
    std::sort(alignments.begin(), alignments.end(),
        [](const alignment& a, const alignment& b) -> bool {
            return a.sw_score > b.sw_score;
        }
    );

    auto max_out = std::min(alignments.size(), static_cast<size_t>(max_secondary));
    bool is_secondary = false;
    for (size_t i = 0; i < max_out; ++i) {
        auto sam_aln = alignments[i];
        if ((sam_aln.sw_score - best_align_sw_score) > (2*aligner.parameters.mismatch + aligner.parameters.gap_open) ){
            break;
        }
        if (is_secondary) {
            sam_aln.mapq = 255;
        } else {
            sam_aln.mapq = std::min(min_mapq_diff, 60);
        }
        sam.add(sam_aln, record, read.rc, is_secondary);
        is_secondary = true;
    }
}

static inline alignment align_segment(
    const Aligner& aligner,
    const std::string &read_segm,
    const std::string &ref_segm,
    int ref_start,
    int ext_left,
    int ext_right,
    bool aln_did_not_fit,
    bool is_rc
) {
    alignment sam_aln_segm;
    auto ref_segm_len = ref_segm.size();
    auto read_segm_len = read_segm.size();
    // The ref_segm includes an extension of ext_left bases upstream and ext_right bases downstream
    int ref_segm_len_ham = ref_segm_len - ext_left - ext_right; // we send in the already extended ref segment to save time. This is not true in center alignment if merged match have diff length
    if (ref_segm_len_ham == read_segm_len && !aln_did_not_fit) {
        std::string ref_segm_ham = ref_segm.substr(ext_left, read_segm_len);

        auto hamming_dist = hamming_distance(read_segm, ref_segm_ham);

        if (hamming_dist >= 0 && (((float) hamming_dist / read_segm_len) < 0.05) ) { //Hamming distance worked fine, no need to ksw align
            auto info = hamming_align(read_segm, ref_segm_ham, aligner.parameters.match, aligner.parameters.mismatch, aligner.parameters.end_bonus);
            sam_aln_segm.cigar = info.cigar.to_string();
            sam_aln_segm.ed = info.edit_distance;
            sam_aln_segm.sw_score = info.sw_score;
            sam_aln_segm.ref_start = ref_start + ext_left + info.query_start;
            sam_aln_segm.is_rc = is_rc;
            sam_aln_segm.is_unaligned = false;
            sam_aln_segm.aln_score = info.sw_score;
            sam_aln_segm.aln_length = read_segm_len;
            return sam_aln_segm;
        }
    }
    auto info = aligner.align(read_segm, ref_segm);
    sam_aln_segm.cigar = info.cigar.to_string();
    sam_aln_segm.ed = info.edit_distance;
    sam_aln_segm.sw_score = info.sw_score;
    sam_aln_segm.ref_start = ref_start + info.ref_start;
    sam_aln_segm.is_rc = is_rc;
    sam_aln_segm.is_unaligned = false;
    sam_aln_segm.aln_score = info.sw_score;
    sam_aln_segm.aln_length = info.ref_span();
    return sam_aln_segm;
}


/*
 Extend a NAM so that it covers the entire read and return the resulting
 alignment.
*/
/*
 Only the following fields of the 'nam' struct are used:
 - is_rc (r/w)
 - ref_id (read)
 - ref_s (read)
 - ref_e (read)
 - query_s (r/w)
 - query_e (r/w)
 This is almost like a 'hit', except for ref_id.

 the nam is sent afterwards into
 - get_MAPQ, which only uses .score and .n_hits
 - rescue_mate, which ...?
*/

static inline alignment get_alignment(
    const Aligner& aligner,
    const Nam &n,
    const References& references,
    const Read& read,
    bool fits
) {
    const std::string query = n.is_rc ? read.rc : read.seq;
    const std::string& ref = references.sequences[n.ref_id];

    const int projected_ref_start = std::max(0, n.ref_s - n.query_s);
    const int projected_ref_end = std::min(n.ref_e + query.size() - n.query_e, ref.size());

    aln_info info;
    int result_ref_start;
    bool has_result = false;
    if (projected_ref_end - projected_ref_start == query.size() && fits) {
        std::string ref_segm_ham = ref.substr(projected_ref_start, query.size());
        auto hamming_dist = hamming_distance(query, ref_segm_ham);

        if (hamming_dist >= 0 && (((float) hamming_dist / query.size()) < 0.05) ) { //Hamming distance worked fine, no need to ksw align
            info = hamming_align(query, ref_segm_ham, aligner.parameters.match, aligner.parameters.mismatch, aligner.parameters.end_bonus);
            result_ref_start = projected_ref_start + info.ref_start;
            has_result = true;
        }
    }
    if (!has_result) {
        const int diff = std::abs(n.ref_span() - n.query_span());
        const int ext_left = std::min(50, projected_ref_start);
        const int ref_start = projected_ref_start - ext_left;
        const int ext_right = std::min(std::size_t(50), ref.size() - n.ref_e);
        const auto ref_segm_size = read.size() + diff + ext_left + ext_right;
        const auto ref_segm = ref.substr(ref_start, ref_segm_size);
        info = aligner.align(query, ref_segm);
        result_ref_start = ref_start + info.ref_start;
    }
    int softclipped = info.query_start + (query.size() - info.query_end);
    alignment sam_aln;
    sam_aln.cigar = info.cigar.to_string();
    sam_aln.ed = info.edit_distance;
    sam_aln.global_ed = info.edit_distance + softclipped;
    sam_aln.sw_score = info.sw_score;
    sam_aln.aln_score = info.sw_score;
    sam_aln.ref_start = result_ref_start;
    sam_aln.aln_length = info.ref_span();
    sam_aln.is_rc = n.is_rc;
    sam_aln.is_unaligned = false;
    sam_aln.ref_id = n.ref_id;
    return sam_aln;
}

static inline alignment get_alignment_unused(
    const Aligner& aligner,
    const Nam &n,
    const References& references,
    const Read& read,
    int k,
    bool fits
) {
    alignment sam_aln;

    const auto read_len = read.size();
    const bool aln_did_not_fit = !fits;
    const int diff = std::abs(n.ref_span() - n.query_span());

    const std::string r_tmp = n.is_rc ? read.rc : read.seq;
    const bool is_rc = n.is_rc;

    int ext_left;
    int ext_right;
    int ref_tmp_segm_size;
    const int ref_len = references.lengths[n.ref_id];
    int ref_tmp_start;
    std::string read_segm;

    // test full hamming based alignment first
    ref_tmp_start = std::max(0, n.ref_s - n.query_s);
    int ref_start = std::max(0, ref_tmp_start);
    ref_tmp_segm_size = read_len + diff;
    auto ref_segm_size = std::min(ref_tmp_segm_size, ref_len - ref_start + 1);

    std::string ref_segm = references.sequences[n.ref_id].substr(ref_start, ref_segm_size);
    if (ref_segm_size == read_len && fits) {
        int hamming_dist = hamming_distance(r_tmp, ref_segm);

        if (hamming_dist >= 0 && (((float) hamming_dist / ref_segm_size) < 0.05) ) { //Hamming distance worked fine, no need to ksw align
            const auto info = hamming_align(r_tmp, ref_segm, aligner.parameters.match, aligner.parameters.mismatch, aligner.parameters.end_bonus);
            sam_aln.cigar = info.cigar.to_string();
            sam_aln.ed = info.edit_distance;
            sam_aln.sw_score = info.sw_score; // aln_params.match*(read_len-hamming_dist) - aln_params.mismatch*hamming_dist;
            sam_aln.ref_start = ref_start + ext_left + info.query_start;
            sam_aln.is_rc = is_rc;
            sam_aln.is_unaligned = false;
            sam_aln.aln_score = info.sw_score;
            return sam_aln;
        }
    }

    //// Didn't work with global Hamming - split into parts

    // identify one or two split points within the read if the segment is  are larger than T
    int T = 20;
    // find the most central split point Use convex function result sum of squares
    int left_outer =  pow (n.query_s, 2) + pow(read_len - n.query_s, 2);
    int right_inner = pow (n.query_e - k, 2) + pow (read_len - (n.query_e - k), 2);


    int global_max_bp = left_outer < right_inner ? n.query_s : n.query_e - k;
    int break_point = (global_max_bp >= T) && (global_max_bp <= (read_len - T)) ? global_max_bp : -1;
    if (break_point > 0 ){
//            std::cerr << "MAX BREAKPOINT " << break_point << " candidates: " <<  n.query_s  << " " << n.query_e - k << std::endl;
        int left_region_bp = break_point + k;
        int right_region_bp = break_point;
//            std::cerr << "left_region_bp " << left_region_bp << " right_region_bp: " << right_region_bp << std::endl;
        int right_ref_start_bp = -1;
        if (break_point == n.query_s){
            right_ref_start_bp = n.ref_s;
        } else if (break_point == (n.query_e - k)) {
            right_ref_start_bp = n.ref_e-k;
        } else  {
            std::cerr << "BUUUUUUG " << std::endl;
        }

        // Left region align
        read_segm = r_tmp.substr(0, left_region_bp);
        ref_tmp_start = std::max(0, n.ref_s - n.query_s);
        ext_left = std::min(50, ref_tmp_start);
        ext_right = 0;
        ref_start = ref_tmp_start - ext_left;
        ref_segm_size = left_region_bp + ext_left + diff;
        ref_segm = references.sequences[n.ref_id].substr(ref_start, ref_segm_size);
//            std::cerr << " "  << std::endl;
//            std::cerr << "GOING IN LEFT: " << " read segm len " << read_segm.length() << " ref segm len " << ref_segm_size  << " ext_left: " << ext_left << std::endl;
//            std::cerr << diff << std::endl;
//            std::cerr << read_segm << std::endl;
//            std::cerr << ref_segm << std::endl;

        auto sam_aln_segm_left = align_segment(aligner, read_segm, ref_segm, ref_start, ext_left, ext_right, aln_did_not_fit, is_rc);
//            std::cerr << "LEFT CIGAR: " << sam_aln_segm_left.cigar << std::endl;

        //Right region align
        read_segm = r_tmp.substr(right_region_bp, read_len - right_region_bp );
        ref_tmp_segm_size = right_ref_start_bp + (read_len + diff - right_region_bp) < ref_len ? (read_len + diff - right_region_bp) : ref_len - right_ref_start_bp;
        ext_left = 0;
        ext_right = std::min(50, ref_len - (right_ref_start_bp + ref_tmp_segm_size));
        ref_segm_size = ref_tmp_segm_size + ext_right;
        ref_segm = references.sequences[n.ref_id].substr(right_ref_start_bp, ref_segm_size);
//            std::cerr << " "  << std::endl;
//            std::cerr << "GOING IN RIGHT: " << " read segm len " << read_segm.length() << " ref segm len " << ref_segm_size  << " ext_right: " << ext_right << std::endl;
//            std::cerr << diff << std::endl;
//            std::cerr << read_segm << std::endl;
//            std::cerr << ref_segm << std::endl;
//            std::cerr << "read_segm.length(): " << read_segm.length() << " ref_segm_size " << ref_segm_size << std::endl;
        auto sam_aln_segm_right = align_segment(aligner, read_segm, ref_segm, ref_start, ext_left, ext_right, aln_did_not_fit, is_rc);
//            std::cerr << "RIGHT CIGAR: " << sam_aln_segm_right.cigar << std::endl;


        // Stitch together
        sam_aln.ref_id = n.ref_id;
        sam_aln.cigar = sam_aln_segm_left.cigar + sam_aln_segm_right.cigar;
        sam_aln.ed = sam_aln_segm_left.ed + sam_aln_segm_right.ed;
        sam_aln.sw_score = sam_aln_segm_left.sw_score + sam_aln_segm_right.sw_score;
        sam_aln.ref_start =  sam_aln_segm_left.ref_start;
        sam_aln.is_rc = n.is_rc;
        sam_aln.is_unaligned = false;
        sam_aln.aln_score = sam_aln.sw_score;
    } else {
//            std::cerr << "NOOOO MAX BREAKPOINT " << break_point << " candidates: "  <<  n.query_s  << " " << n.query_e - k << std::endl;
        // full align
        ref_tmp_start = std::max(0, n.ref_s - n.query_s);
        ext_left = std::min(50, ref_tmp_start);
        ref_start = ref_tmp_start - ext_left;

        ref_tmp_segm_size = read_len + diff;
        ext_right = std::min(50, ref_len - (n.ref_e +1));

        ref_segm_size = ref_tmp_segm_size + ext_left + ext_right;
        ref_segm = references.sequences[n.ref_id].substr(ref_start, ref_segm_size);
//        std::cerr << " ref_tmp_start " << ref_tmp_start << " ext left " << ext_left << " ext right " << ext_right << " ref_tmp_segm_size " << ref_tmp_segm_size << " ref_segm_size " << ref_segm_size << " ref_segm " << ref_segm << std::endl;
        sam_aln = align_segment(aligner, r_tmp, ref_segm, ref_start, ext_left, ext_right, aln_did_not_fit, is_rc);
        sam_aln.ref_id = n.ref_id;
    }


        // TODO: Several breakpoints. To complicated for probably little gain if short reads, maybe implement later..
//        // Left and right breakpoint
//        int left_bp_read = n.query_s + k;
//        int right_bp_read = n.query_e - k;
//        left_bp_read = left_bp_read > T ? left_bp_read : -1;
//        right_bp_read = read_len - right_bp_read > T ? read_len - right_bp_read : -1;
//
//        // Decide where to break
//        bool break1, break2;
//        if ( ((read_len - right_bp_read - left_bp_read) > T) && (left_bp_read >= T) && (right_bp_read >= T) ){ // Use two breakpoints ---------------X|------------------X|---------------
//            break1 = true;
//            break2 = true;
//            std::cerr << "CASE1 " << std::endl;
//        } else if (  (left_bp_read >= T) && (right_bp_read >= T) ) { // Use one breakpoint  ---------------------------X--------X--------------------------
//            break1 = left_bp_read >= right_bp_read ? true : false;  // ------------------------------------X|--------X--------------------
//            break2 = left_bp_read < right_bp_read ? true : false;   // ---------------------X--------|X-----------------------------------
//            std::cerr << "CASE2 " << std::endl;
//        } else if ( (left_bp_read >= T) && ( (read_len - left_bp_read) >= T) ) { // Use one breakpoint  -----------------------X|-------------------X-------
//            break1 = true;
//            break2 = false;
//            std::cerr << "CASE3 " << std::endl;
//        } else if ( ((right_bp_read) > T) && (read_len - right_bp_read >= T)  ) { // Use one breakpoint  ----------X-------------------X|---------------------------
//            break1 = false;
//            break2 = true;
//            std::cerr << "CASE4 " << std::endl;
//        } else { // No breakpoints ------X--------------------------------------X--------
//            break1 = false;
//            break2 = false;
//            std::cerr << "CASE5 " << std::endl;
//        }
//
//        std::cerr << "BREAKPOINTS: LEFT:  " << left_bp_read << " RIGHT:  " << right_bp_read << " BREAKING LEFT:  " << break1 << " BREAKING RIGHT: " << break2 << std::endl;
//        std::cerr << "MAPPING POS n.ref_s:  " << n.ref_s << " n.ref_e:  " << n.ref_e << " n.query_s  " << n.query_s << " n.query_e: " << n.query_e << std::endl;
//        std::cerr << " " << std::endl;


//        // left alignment
//        ref_end = n.ref_s + k;
//        ref_start = ref_end - ext_left;
//        ext_left = ref_end < 50 ? ref_end : 50;
//
//        ref_tmp_segm_size = n.query_s + k;
//        ext_right = 0;
//
//        ref_segm_size = ref_tmp_segm_size + ext_left + ext_right;
//        ref_segm = ref_seqs[n.ref_id].substr(ref_end - ext_left, ref_segm_size);
//        alignment sam_aln_segm_left;
//        sam_aln_segm_left.ref_id = n.ref_id;
//        read_segm = r_tmp.substr(0, n.query_s+k);
//
//        std::cerr << "LEFT: ref_start:  " << ref_start << " LEFT: ref_end:  " << ref_end << " ref_segm_size  " << ref_segm_size << " read segm size: " << read_segm.length() << std::endl;
//
//        align_segment(aln_params, read_segm, ref_segm, read_segm.length(), ref_segm_size, ref_start, ext_left, ext_right, aln_did_not_fit, is_rc, sam_aln_segm_left, tot_ksw_aligned);
//
//
//        // center alignment
//        ref_tmp_start = n.ref_s - n.query_s;
//        ref_start = ref_tmp_start > 0 ? ref_tmp_start : 0;
//        ext_left = 0;
//
//        ref_tmp_segm_size =  n.ref_e - n.ref_s;
//        ext_right = 0;
//
//        ref_segm_size = ref_tmp_segm_size + ext_left + ext_right;
//        ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_segm_size);
//        alignment sam_aln_segm_center;
//        sam_aln_segm_center.ref_id = n.ref_id;
//        read_segm = r_tmp.substr(n.query_s, n.query_e - n.query_s);
//        std::cerr << "CENTER: ref_start:  " << ref_start <<  " CENTER: ref_end:  " << ref_end << " ref_segm_size  " << ref_segm_size << " read segm size: " << read_segm.length() << std::endl;
//
//        align_segment(aln_params, read_segm, ref_segm, read_segm.length(), ref_segm_size, ref_start, ext_left, ext_right, aln_did_not_fit, is_rc, sam_aln_segm_center, tot_ksw_aligned);
//
//
//        // right alignment
//        ext_left = 0;
//        ref_start = n.ref_e - k;
//
//        ext_right = ref_len - (n.ref_e +1) < 50 ? ref_len - (n.ref_e +1) : 50;
//        ref_tmp_segm_size = n.ref_e + ext_right - ref_start;
//
//        ref_segm_size = ref_tmp_segm_size + ext_left + ext_right;
//        ref_segm = ref_seqs[n.ref_id].substr(n.ref_e  - ext_left, ref_segm_size);
//        alignment sam_aln_segm_right;
//        sam_aln_segm_left.ref_id = n.ref_id;
//        read_segm = r_tmp.substr(0, n.query_s+k);
//        std::cerr << "RIGHT: ref_start:  " << ref_start <<  " RIGHT: ref_end:  " << ref_end << " ref_segm_size  " << ref_segm_size << " read segm size: " << read_segm.length() << std::endl;
//
//        align_segment(aln_params, read_segm, ref_segm, read_segm.length(), ref_segm_size, ref_start, ext_left, ext_right, aln_did_not_fit, is_rc, sam_aln_segm_right, tot_ksw_aligned);
//
//        std::cout << sam_aln_segm_left.cigar << " " << sam_aln_segm_center.cigar << " " << sam_aln_segm_right.cigar << std::endl;
//
//        sam_aln.ref_id = n.ref_id;
//        sam_aln.cigar = sam_aln_segm_left.cigar + sam_aln_segm_center.cigar + sam_aln_segm_right.cigar;
//        sam_aln.ed = sam_aln_segm_left.ed + sam_aln_segm_center.ed + sam_aln_segm_right.ed;
//        sam_aln.sw_score = sam_aln_segm_left.sw_score + sam_aln_segm_center.sw_score + sam_aln_segm_right.sw_score;
//        sam_aln.ref_start =   sam_aln_segm_left.ref_start;
//        sam_aln.is_rc = sam_aln_segm_left.is_rc;
//        sam_aln.is_unaligned = false;
//        sam_aln.aln_score = sam_aln.sw_score;
//        std::cout << "Joint: " << sam_aln.cigar << std::endl;

    return sam_aln;
}




static inline int get_MAPQ(const std::vector<Nam> &all_nams, const Nam &n_max) {
    const float s1 = n_max.score;
    if (all_nams.size() <= 1) {
        return 60;
    }
    const Nam n_second = all_nams[1];
    const float s2 = n_second.score;
    // from minimap2: MAPQ = 40(1−s2/s1) ·min{1,|M|/10} · log s1
    const float min_matches = std::min(n_max.n_hits / 10.0, 1.0);
    const int uncapped_mapq = 40 * (1 - s2 / s1) * min_matches * log(s1);
    return std::min(uncapped_mapq, 60);
}


static inline std::pair<int, int> joint_mapq_from_alignment_scores(float S1, float S2) {
    int mapq;
    if (S1 == S2) { // At least two identical placements
        mapq = 0;
    } else {
        const int diff = S1 - S2; // (1.0 - (S1 - S2) / S1);
//        float log10_p = diff > 6 ? -6.0 : -diff; // Corresponds to: p_error= 0.1^diff // change in sw score times rough illumina error rate. This is highly heauristic, but so seem most computations of mapq scores
        if (S1 > 0 && S2 > 0) {
            mapq = std::min(60, diff);
//            mapq1 = -10 * log10_p < 60 ? -10 * log10_p : 60;
        } else if (S1 > 0 && S2 <= 0) {
            mapq = 60;
        } else { // both negative SW one is better
            mapq = 1;
        }
    }
    return std::make_pair(mapq, mapq);
}


static inline float normal_pdf(float x, float m, float s)
{
    static const float inv_sqrt_2pi = 0.3989422804014327;
    const float a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

static inline bool score_sw(const alignment &a, const alignment &b)
{
    return a.sw_score > b.sw_score;
}

static inline bool sort_scores(const std::tuple<double, alignment, alignment> &a,
                               const std::tuple<double, alignment, alignment> &b)
{
    return std::get<0>(a) > std::get<0>(b);
}




static inline void get_best_scoring_pair(
    const std::vector<alignment> &aln_scores1, const std::vector<alignment> &aln_scores2, std::vector<std::tuple<double,alignment,alignment>> &high_scores, float mu, float sigma
) {
    for (auto &a1 : aln_scores1) {
        for (auto &a2 : aln_scores2) {
            float dist = std::abs(a1.ref_start - a2.ref_start);
            double score = a1.sw_score + a2.sw_score;
            if ((a1.is_rc ^ a2.is_rc) && (dist < mu + 4 * sigma)) {
                score += log(normal_pdf(dist, mu, sigma));
            }
            else { // individual score
                // 10 corresponds to a value of log(normal_pdf(dist, mu, sigma)) of more than 4 stddevs away
                score -= 10;
            }
            high_scores.push_back(std::make_tuple(score, a1, a2));
        }
    }
    std::sort(high_scores.begin(), high_scores.end(), sort_scores); // Sorting by highest score first
}

static inline std::vector<std::tuple<int,Nam,Nam>> get_best_scoring_NAM_locations(
    const std::vector<Nam> &all_nams1,
    const std::vector<Nam> &all_nams2,
    float mu,
    float sigma
) {
    std::vector<std::tuple<int,Nam,Nam>> joint_NAM_scores;
    if (all_nams1.empty() && all_nams2.empty()) {
        return joint_NAM_scores;
    }

    robin_hood::unordered_set<int> added_n1;
    robin_hood::unordered_set<int> added_n2;
    int joint_hits;
    int hjss = 0; // highest joint score seen
//            std::cerr << "Scoring" << std::endl;
    int a,b;
    for (auto &n1 : all_nams1) {
        for (auto &n2 : all_nams2) {
            if ((n1.n_hits + n2.n_hits) < hjss/2){
                break;
            }

            if ((n1.is_rc ^ n2.is_rc) && (n1.ref_id == n2.ref_id)) {
                a = std::max(0, n1.ref_s - n1.query_s);
                b = std::max(0, n2.ref_s - n2.query_s);
                bool r1_r2 = n2.is_rc && (a < b) && ((b-a) < mu+10*sigma); // r1 ---> <---- r2
                bool r2_r1 = n1.is_rc && (b < a) && ((a-b) < mu+10*sigma); // r2 ---> <---- r1

                if ( r1_r2 || r2_r1 ){
//                    int diff1 = (n1.query_e - n1.query_s) - (n1.ref_e - n1.ref_s);
//                    int  n1_penalty = diff1 > 0 ? diff1 : - diff1;
//                    int diff2 = (n2.query_e - n2.query_s) - (n2.ref_e - n2.ref_s);
//                    int  n2_penalty = diff2 > 0 ? diff2 : - diff2;
                    joint_hits = n1.n_hits + n2.n_hits; // - n1_penalty - n2_penalty; // trying out idea about penalty but it needs to be on the individual seed level - to late on merged match level.
                    std::tuple<int, Nam, Nam> t (joint_hits, n1, n2);
                    joint_NAM_scores.push_back(t);
                    added_n1.insert(n1.nam_id);
                    added_n2.insert(n2.nam_id);
                    if (joint_hits > hjss) {
                        hjss = joint_hits;
                    }
                }
            }
        }
    }

    Nam dummy_nan;
    dummy_nan.ref_s = -1;
    if (!all_nams1.empty()) {
        int hjss1 = hjss > 0 ? hjss : all_nams1[0].n_hits;
        for (auto &n1 : all_nams1) {
            if (n1.n_hits  < hjss1/2){
                break;
            }
            if (added_n1.find(n1.nam_id) != added_n1.end()){
                continue;
            }
//            int diff1 = (n1.query_e - n1.query_s) - (n1.ref_e - n1.ref_s);
//            int  n1_penalty = diff1 > 0 ? diff1 : - diff1;
            joint_hits = n1.n_hits;
            std::tuple<int, Nam, Nam> t{joint_hits, n1, dummy_nan};
            joint_NAM_scores.push_back(t);
        }
    }

    if ( !all_nams2.empty() ){
        int hjss2 = hjss  > 0 ? hjss : all_nams2[0].n_hits;
        //    int hjss2 = all_nams2[0].n_hits;
        for (auto &n2 : all_nams2) {
            if (n2.n_hits  < hjss2/2){
                break;
            }
            if (added_n2.find(n2.nam_id) != added_n2.end()){
                continue;
            }
//            int diff2 = (n2.query_e - n2.query_s) - (n2.ref_e - n2.ref_s);
//            int  n2_penalty = diff2 > 0 ? diff2 : - diff2;
            joint_hits = n2.n_hits;
            //                        std::cerr << S << " individual score " << x << " " << std::endl;
            std::tuple<int, Nam, Nam> t{joint_hits, dummy_nan, n2};
            joint_NAM_scores.push_back(t);
        }
    }

    added_n1.clear();
    added_n2.clear();

    std::sort(
        joint_NAM_scores.begin(),
        joint_NAM_scores.end(),
        [](const std::tuple<int, Nam, Nam>& a, const std::tuple<int, Nam, Nam>& b) -> bool {
            return std::get<0>(a) > std::get<0>(b);
        }
    ); // Sort by highest score first

    return joint_NAM_scores;
}

/*
 * Determine (roughly) whether the read sequence has some l-mer (with l = k*2/3)
 * in common with the reference sequence
 */
bool has_shared_substring(const std::string& read_seq, const std::string& ref_seq, int k) {
    int sub_size = 2 * k / 3;
    int step_size = k / 3;
    std::string submer;
    for (size_t i = 0; i + k <= read_seq.size(); i += step_size) {
        submer = read_seq.substr(i, sub_size);
        if (ref_seq.find(submer) != std::string::npos) {
            return true;
        }
    }
    return false;
}

static inline void rescue_mate(
    const Aligner& aligner,
    Nam &n,
    const References& references,
    const Read& guide,
    const Read& read,
    alignment &sam_aln,
    float mu,
    float sigma,
    unsigned int &tot_rescued,
    int k
) {
    int a, b, ref_start,ref_len,ref_end;
    std::string r_tmp;
    bool a_is_rc;
    auto read_len = read.size();

    reverse_nam_if_needed(n, guide, references, k);
    if (n.is_rc){
        r_tmp = read.seq;
        a = n.ref_s - n.query_s - (mu+5*sigma);
        b = n.ref_s - n.query_s + read_len/2; // at most half read overlap
        a_is_rc = false;
    } else {
        r_tmp = read.rc; // mate is rc since fr orientation
        a = n.ref_e + (read_len - n.query_e) - read_len/2; // at most half read overlap
        b = n.ref_e + (read_len - n.query_e) + (mu+5*sigma);
        a_is_rc = true;
    }

    ref_len = references.lengths[n.ref_id];
    ref_start = std::max(0, std::min(a,ref_len));
    ref_end = std::min(ref_len, std::max(0, b));

    if (ref_end < ref_start + k){
        sam_aln.cigar = "*";
        sam_aln.ed = read_len;
        sam_aln.sw_score = 0;
        sam_aln.aln_score = 0;
        sam_aln.ref_start =  0;
        sam_aln.is_rc = n.is_rc;
        sam_aln.ref_id = n.ref_id;
        sam_aln.is_unaligned = true;
//        std::cerr << "RESCUE: Caught Bug3! ref start: " << ref_start << " ref end: " << ref_end << " ref len:  " << ref_len << std::endl;
        return;
    }
    std::string ref_segm = references.sequences[n.ref_id].substr(ref_start, ref_end - ref_start);

    if (!has_shared_substring(r_tmp, ref_segm, k)){
        sam_aln.cigar = "*";
        sam_aln.ed = read_len;
        sam_aln.sw_score = 0;
        sam_aln.aln_score = 0;
        sam_aln.ref_start =  0;
        sam_aln.is_rc = n.is_rc;
        sam_aln.ref_id = n.ref_id;
        sam_aln.is_unaligned = true;
//        std::cerr << "Avoided!" << std::endl;
        return;
//        std::cerr << "Aligning anyway at: " << ref_start << " to " << ref_end << "ref len:" << ref_len << " ref_id:" << n.ref_id << std::endl;
    }
    auto info = aligner.align(r_tmp, ref_segm);

//    if (info.ed == 100000){
//        std::cerr<< "________________________________________" << std::endl;
//        std::cerr<< "RESCUE MODE: " << mu << "  " << sigma << std::endl;
//        std::cerr<< read << "   " << read_rc << std::endl;
//        std::cerr << r_tmp << " " << n.n_hits << " " << n.score << " " << " " << sam_aln.ed << " "  <<  n.query_s << " "  << n.query_e << " "<<  n.ref_s << " "  << n.ref_e << " " << n.is_rc << " " << " " << sam_aln.cigar << " " << info.sw_score << std::endl;
//        std::cerr << "a " << a << " b " << b << " ref_start " <<  ref_start << " ref_end " << ref_end << "  ref_end - ref_start "  <<  ref_end - ref_start << "  n.is_flipped " <<  n.is_flipped << std::endl;
//        std::cerr<< "________________________________________" << std::endl;
//    }
//    info = parasail_align(ref_segm, ref_segm.size(), r_tmp, read_len, 1, 4, 6, 1);

//    ksw_extz_t ez;
//    const char *ref_ptr = ref_segm.c_str();
//    const char *read_ptr = r_tmp.c_str();
//    info = ksw_align(ref_ptr, ref_segm.size(), read_ptr, r_tmp.size(), 1, 4, 6, 1, ez);
//    std::cerr << "Cigar: " << info.cigar << std::endl;

    sam_aln.cigar = info.cigar.to_string();
    sam_aln.ed = info.edit_distance;
    sam_aln.sw_score = info.sw_score;
    sam_aln.aln_score = sam_aln.sw_score;
    sam_aln.ref_start =  ref_start + info.ref_start;
    sam_aln.is_rc = a_is_rc;
    sam_aln.ref_id = n.ref_id;
    sam_aln.is_unaligned = info.cigar.empty();
    sam_aln.aln_length = info.ref_span();
    tot_rescued ++;
}


void rescue_read(
    const Read& read2,  // read to be rescued
    const Read& read1,  // read that has NAMs
    const Aligner& aligner,
    const References& references,
    std::vector<Nam> &all_nams1,
    int max_tries,
    float dropoff,
    AlignmentStatistics &statistics,
    int k,
    float mu,
    float sigma,
    size_t max_secondary,
    double secondary_dropoff,
    Sam& sam,
    const klibpp::KSeq& record1,
    const klibpp::KSeq& record2,
    bool swap_r1r2  // TODO get rid of this
) {
    float score_dropoff1;
    Nam n_max1 = all_nams1[0];
    int cnt1 = 0;

    std::vector<alignment> aln_scores1;
    std::vector<alignment> aln_scores2;
    for (auto& n : all_nams1) {
        score_dropoff1 = (float) n.n_hits / n_max1.n_hits;
        if ( (cnt1 >= max_tries) || score_dropoff1 < dropoff){ // only consider top 20 hits as minimap2 and break if alignment is exact match to reference or the match below droppoff cutoff.
            break;
        }
        //////// the actual testing of base pair alignment part start /////////
//            std::cerr << query_acc1 << " force rescue"  << std::endl;

        const bool fits = reverse_nam_if_needed(n, read1, references, k);
        if (!fits) {
            statistics.did_not_fit++;
        }
        alignment a1 = get_alignment(aligner, n, references, read1, fits);
        aln_scores1.push_back(a1);

        //////// Force SW alignment to rescue mate /////////
        alignment a2;
//            std::cerr << query_acc2 << " force rescue" << std::endl;
        rescue_mate(aligner, n, references, read1, read2, a2, mu, sigma, statistics.tot_rescued, k);
        aln_scores2.push_back(a2);

        cnt1 ++;
        statistics.tot_all_tried ++;
    }
    std::sort(aln_scores1.begin(), aln_scores1.end(), score_sw);
    std::sort(aln_scores2.begin(), aln_scores2.end(), score_sw);

    // Calculate best combined score here
    std::vector<std::tuple<double,alignment,alignment>> high_scores; // (score, aln1, aln2)
    get_best_scoring_pair(aln_scores1, aln_scores2, high_scores, mu, sigma );

    // Calculate joint MAPQ score
    int mapq1, mapq2;
    if (high_scores.size() > 1) {
        auto best_aln_pair = high_scores[0];
        auto S1 = std::get<0>(best_aln_pair);
        auto second_aln_pair = high_scores[1];
        auto S2 = std::get<0>(second_aln_pair);
        std::tie(mapq1, mapq2) = joint_mapq_from_alignment_scores(S1, S2);
    } else {
        mapq1 = 60;
        mapq2 = 60;
    }

    // append both alignments to string here
    if (max_secondary == 0){
        auto best_aln_pair = high_scores[0];
        alignment sam_aln1 = std::get<1>(best_aln_pair);
        alignment sam_aln2 = std::get<2>(best_aln_pair);
//            get_MAPQ(all_nams1, n_max1, mapq1);
//            mapq2 = 0;
        if (swap_r1r2) {
            sam.add_pair(sam_aln2, sam_aln1, record2, record1, read2.rc, read1.rc, mapq2, mapq1, is_proper_pair(sam_aln2, sam_aln1, mu, sigma), true);
        } else {
            sam.add_pair(sam_aln1, sam_aln2, record1, record2, read1.rc, read2.rc, mapq1, mapq2, is_proper_pair(sam_aln1, sam_aln2, mu, sigma), true);
        }
    } else {
        int max_out = std::min(high_scores.size(), max_secondary);
        bool is_primary = true;
        auto best_aln_pair = high_scores[0];
        auto s_max = std::get<0>(best_aln_pair);
//            get_MAPQ(all_nams1, n_max1, mapq1);
        for (int i = 0; i < max_out; ++i) {
            if (i > 0) {
                is_primary = false;
                mapq1 = 0;
                mapq2 = 0;
            }
            auto aln_pair = high_scores[i];
            auto s_score = std::get<0>(aln_pair);
            alignment sam_aln1 = std::get<1>(aln_pair);
            alignment sam_aln2 = std::get<2>(aln_pair);
            if (s_max - s_score < secondary_dropoff) {
                if (swap_r1r2) {
                    bool is_proper = is_proper_pair(sam_aln2, sam_aln1, mu, sigma);
                    sam.add_pair(sam_aln2, sam_aln1, record2, record1, read2.rc, read1.rc, mapq2, mapq1, is_proper, is_primary);
                } else {
                    bool is_proper = is_proper_pair(sam_aln1, sam_aln2, mu, sigma);
                    sam.add_pair(sam_aln1, sam_aln2, record1, record2, read1.rc, read2.rc, mapq1, mapq2, is_proper, is_primary);
                }
            } else {
                break;
            }
        }
    }
}

/* Compute paired-end mapping score given top alignments */
static std::pair<int, int> joint_mapq_from_high_scores(const std::vector<std::tuple<double,alignment,alignment>>& high_scores) {
    if (high_scores.size() <= 1) {
        return std::make_pair(60, 60);
    }
    // Calculate joint MAPQ score
    int n_mappings = high_scores.size();
    auto best_aln_pair = high_scores[0];
    auto S1 = std::get<0>(best_aln_pair);
    auto a1_m1 = std::get<1>(best_aln_pair);
    auto a1_m2 = std::get<2>(best_aln_pair);
    int a1_start_m1 = a1_m1.ref_start;
    int a1_start_m2 = a1_m2.ref_start;
    int a1_ref_id_m1 = a1_m1.ref_id;
    int a1_ref_id_m2 = a1_m2.ref_id;

    auto second_aln_pair = high_scores[1];
    auto S2 = std::get<0>(second_aln_pair);
    auto a2_m1 = std::get<1>(second_aln_pair);
    auto a2_m2 = std::get<2>(second_aln_pair);
    int a2_start_m1 = a2_m1.ref_start;
    int a2_start_m2 = a2_m2.ref_start;
    int a2_ref_id_m1 = a2_m1.ref_id;
    int a2_ref_id_m2 = a2_m2.ref_id;
    bool same_pos = (a1_start_m1 == a2_start_m1) && (a1_start_m2 == a2_start_m2);
    bool same_ref = (a1_ref_id_m1 == a2_ref_id_m1) && (a1_ref_id_m2 == a2_ref_id_m2);
    if (!same_pos || !same_ref) {
        return joint_mapq_from_alignment_scores(S1, S2);
    } else if (n_mappings > 2) {
        // individually highest alignment score was the same alignment as the joint highest score - calculate mapq relative to third best
        auto third_aln_pair = high_scores[2];
        auto S2 = std::get<0>(third_aln_pair);
        return joint_mapq_from_alignment_scores(S1, S2);
    } else {
        // there was no other alignment
        return std::make_pair(60, 60);
    }
}

inline void align_PE(
    const Aligner& aligner,
    Sam &sam,
    std::vector<Nam> &all_nams1,
    std::vector<Nam> &all_nams2,
    const KSeq &record1,
    const KSeq &record2,
    int k,
    const References& references,
    AlignmentStatistics &statistics,
    float dropoff,
    i_dist_est &isize_est,
    int max_tries,
    size_t max_secondary
) {
    const auto mu = isize_est.mu;
    const auto sigma = isize_est.sigma;
    Read read1(record1.seq);
    Read read2(record2.seq);
    double secondary_dropoff = 2 * aligner.parameters.mismatch + aligner.parameters.gap_open;

    if (all_nams1.empty() && all_nams2.empty()) {
         // None of the reads have any NAMs
        sam.add_unmapped_pair(record1, record2);
        return;
    }

    if (!all_nams1.empty() && all_nams2.empty()) {
        // Only read 1 has NAMS: attempt to rescue read 2
        rescue_read(
            read2,
            read1,
            aligner,
            references,
            all_nams1,
            max_tries,
            dropoff,
            statistics,
            k,
            mu,
            sigma,
            max_secondary,
            secondary_dropoff,
            sam,
            record1,
            record2,
            false
        );
        return;
    }

    if (all_nams1.empty() && !all_nams2.empty()) {
        // Only read 2 has NAMS: attempt to rescue read 1
        rescue_read(
            read1,
            read2,
            aligner,
            references,
            all_nams2,
            max_tries,
            dropoff,
            statistics,
            k,
            mu,
            sigma,
            max_secondary,
            secondary_dropoff,
            sam,
            record2,
            record1,
            true
        );
        return;
    }

    // If we get here, both reads have NAMs
    assert(!all_nams1.empty() && !all_nams2.empty());

    int cnt = 0;
    double S = 0.0;
    Nam n_max1 = all_nams1[0];
    Nam n_max2 = all_nams2[0];

    float score_dropoff1 = all_nams1.size() > 1 ? (float) all_nams1[1].n_hits / n_max1.n_hits : 0.0;
    float score_dropoff2 = all_nams2.size() > 1 ? (float) all_nams2[1].n_hits / n_max2.n_hits : 0.0;
    score_dropoff1 = n_max1.n_hits > 2 ? score_dropoff1 : 1.0;
    score_dropoff2 = n_max2.n_hits > 2 ? score_dropoff2 : 1.0;

    int a = std::max(0, n_max1.ref_s - n_max1.query_s);
    int b = std::max(0, n_max2.ref_s - n_max2.query_s);
    bool r1_r2 = n_max2.is_rc && (a < b) && ((b-a) < 2000); // r1 ---> <---- r2
    bool r2_r1 = n_max1.is_rc && (b < a) && ((a-b) < 2000); // r2 ---> <---- r1

    if (score_dropoff1 < dropoff && score_dropoff2 < dropoff && (n_max1.is_rc ^ n_max2.is_rc) && (r1_r2 || r2_r1)) { //( ((n_max1.ref_s - n_max2.ref_s) < mu + 4*sigma ) || ((n_max2.ref_s - n_max1.ref_s ) < mu + 4*sigma ) ) &&

        bool fits1 = reverse_nam_if_needed(n_max1, read1, references, k);
        if (!fits1) {
            statistics.did_not_fit++;
        }
        bool fits2 = reverse_nam_if_needed(n_max2, read2, references, k);
        if (!fits2) {
            statistics.did_not_fit++;
        }

        auto sam_aln1 = get_alignment(aligner, n_max1, references, read1, fits1);
        statistics.tot_all_tried ++;
        auto sam_aln2 = get_alignment(aligner, n_max2, references, read2, fits2);
        statistics.tot_all_tried ++;
        int mapq1 = get_MAPQ(all_nams1, n_max1);
        int mapq2 = get_MAPQ(all_nams2, n_max2);
        bool is_proper = is_proper_pair(sam_aln1, sam_aln2, mu, sigma);
        sam.add_pair(sam_aln1, sam_aln2, record1, record2, read1.rc, read2.rc, mapq1, mapq2, is_proper, true);

        if ((isize_est.sample_size < 400) && ((sam_aln1.ed + sam_aln2.ed) < 3) && is_proper) {
            isize_est.update(std::abs(sam_aln1.ref_start - sam_aln2.ref_start));
        }
        return;
    }

    // do full search of highest scoring pair
    // std::cerr << "Joint search" << std::endl;

    //////////////////////////// NEW ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    // Get top hit counts for all locations. The joint hit count is the sum of hits of the two mates. Then align as long as score dropoff or cnt < 20

    // (score, aln1, aln2)
    std::vector<std::tuple<int,Nam,Nam>> joint_NAM_scores = get_best_scoring_NAM_locations(all_nams1, all_nams2, mu, sigma);
    auto nam_max = joint_NAM_scores[0];
    auto max_score = std::get<0>(nam_max);

    robin_hood::unordered_map<int,alignment> is_aligned1;
    robin_hood::unordered_map<int,alignment> is_aligned2;
    auto n1_max = all_nams1[0];

    bool fits1 = reverse_nam_if_needed(n1_max, read1, references, k);
    if (!fits1) {
        statistics.did_not_fit++;
    }
    auto a1_indv_max = get_alignment(aligner, n1_max, references, read1,
                    fits1);
//            a1_indv_max.sw_score = -10000;
    is_aligned1[n1_max.nam_id] = a1_indv_max;
    statistics.tot_all_tried ++;
    auto n2_max = all_nams2[0];
    bool fits2 = reverse_nam_if_needed(n2_max, read2, references, k);
    if (!fits2) {
        statistics.did_not_fit++;
    }
    auto a2_indv_max = get_alignment(aligner, n2_max, references, read2,
                    fits2);
//            a2_indv_max.sw_score = -10000;
    is_aligned2[n2_max.nam_id] = a2_indv_max;
    statistics.tot_all_tried ++;

//            int a, b;
    std::string r_tmp;
//            int min_ed1, min_ed2 = 1000;
//            bool new_opt1, new_opt2 = false;
//            bool a1_is_rc, a2_is_rc;
//            int ref_start, ref_len, ref_end;
//            std::cerr << "LOOOOOOOOOOOOOOOOOOOL " << min_ed << std::endl;
    std::vector<std::tuple<double,alignment,alignment>> high_scores; // (score, aln1, aln2)
    for (auto &[score_, n1, n2] : joint_NAM_scores) {
        score_dropoff1 = (float) score_ / max_score;
        if ( (cnt >= max_tries) || (score_dropoff1 < dropoff) ){ // only consider top 20 if there are more.
            break;
        }

        //////// the actual testing of base pair alignment part start ////////
        //////////////////////////////////////////////////////////////////////
        alignment a1;
        if (n1.ref_s >= 0) {
            if (is_aligned1.find(n1.nam_id) != is_aligned1.end() ){
                a1 = is_aligned1[n1.nam_id];
            } else {
                bool fits = reverse_nam_if_needed(n1, read1, references, k);
                if (!fits) {
                    statistics.did_not_fit++;
                }
                a1 = get_alignment(aligner, n1, references, read1, fits);
                is_aligned1[n1.nam_id] = a1;
                statistics.tot_all_tried++;
            }
        } else {
            //////// Force SW alignment to rescue mate /////////
//                    std::cerr << query_acc2 << " RESCUE MATE 1" << a1.is_rc << " " n1.is_rc << std::endl;
            rescue_mate(aligner, n2, references, read2, read1, a1, mu, sigma, statistics.tot_rescued, k);
//                    is_aligned1[n1.nam_id] = a1;
            statistics.tot_all_tried ++;
        }


//                a1_indv_max = a1.sw_score >  a1_indv_max.sw_score ? a1 : a1_indv_max;
//                min_ed = a1.ed < min_ed ? a1.ed : min_ed;

        if (a1.sw_score >  a1_indv_max.sw_score){
            a1_indv_max = a1;
//                    cnt = 0;
        }

        alignment a2;
        if(n2.ref_s >= 0) {
            if (is_aligned2.find(n2.nam_id) != is_aligned2.end() ){
//                    std::cerr << "Already aligned a2! " << std::endl;
                a2 = is_aligned2[n2.nam_id];
            } else {
//                    std::cerr << query_acc2 << std::endl;
                bool fits = reverse_nam_if_needed(n2, read2, references, k);
                if (!fits) {
                    statistics.did_not_fit++;
                }
                a2 = get_alignment(aligner, n2, references, read2,
                                fits);
                is_aligned2[n2.nam_id] = a2;
                statistics.tot_all_tried++;
            }
        } else{
//                    std::cerr << "RESCUE HERE2" << std::endl;
            //////// Force SW alignment to rescue mate /////////
//                    std::cerr << query_acc1 << " RESCUE MATE 2" << a1.is_rc << " " n1.is_rc << std::endl;
            rescue_mate(aligner, n1, references, read1, read2, a2, mu, sigma, statistics.tot_rescued, k);
//                    is_aligned2[n2.nam_id] = a2;
            statistics.tot_all_tried ++;
        }
//                a2_indv_max = a2.sw_score >  a2_indv_max.sw_score ? a2 : a2_indv_max;
//                min_ed = a2.ed < min_ed ? a2.ed : min_ed;

        if (a2.sw_score >  a2_indv_max.sw_score){
            a2_indv_max = a2;
//                    cnt = 0;
        }

        bool r1_r2 = a2.is_rc && (a1.ref_start < a2.ref_start) && ((a2.ref_start - a1.ref_start) < mu+5*sigma); // r1 ---> <---- r2
        bool r2_r1 = a1.is_rc && (a2.ref_start < a1.ref_start) && ((a1.ref_start - a2.ref_start) < mu+5*sigma); // r2 ---> <---- r1

        if (r1_r2 || r2_r1) {
            float x = std::abs(a1.ref_start - a2.ref_start);
            S = (double)a1.sw_score + (double)a2.sw_score + log(normal_pdf(x, mu, sigma));  //* (1 - s2 / s1) * min_matches * log(s1);
//                    std::cerr << " CASE1: " << S << " " <<  log( normal_pdf(x, mu, sigma ) ) << " " << (double)a1.sw_score << " " << (double)a2.sw_score << std::endl;
        } else{ // individual score
            S = (double)a1.sw_score + (double)a2.sw_score - 20; // 20 corresponds to a value of log( normal_pdf(x, mu, sigma ) ) of more than 5 stddevs away (for most reasonable values of stddev)
//                    std::cerr << " CASE2: " << S << " " << (double)a1.sw_score << " " << (double)a2.sw_score << std::endl;
        }

        std::tuple<double, alignment, alignment> aln_tuple (S, a1, a2);
        high_scores.push_back(aln_tuple);

        cnt ++;
    }

    // Finally, add highest scores of both mates as individually mapped
    S = (double)a1_indv_max.sw_score + (double)a2_indv_max.sw_score - 20; // 20 corresponds to  a value of log( normal_pdf(x, mu, sigma ) ) of more than 5 stddevs away (for most reasonable values of stddev)
    std::tuple<double, alignment, alignment> aln_tuple (S, a1_indv_max, a2_indv_max);
    high_scores.push_back(aln_tuple);
    std::sort(high_scores.begin(), high_scores.end(), sort_scores); // Sorting by highest score first

//            if (mapq1 != 60){
//                std::cerr << query_acc1 << " " << mapq1 << std::endl;
//            }

//            std::cerr << x << " " << mu << " " << sigma << " " << log( normal_pdf(x, mu, sigma ) ) << std::endl;
//            std::cerr << 200 << " " << 200 << " " << 30 << " " << log( normal_pdf(200, 200, 30 ) ) << std::endl;
//            std::cerr << 200 << " " << 200 << " " << 200 << " " << log( normal_pdf(200, 200, 200 ) ) << std::endl;
//            std::cerr << 350 << " " << 200 << " " << 30 << " " << log( normal_pdf(350, 200, 30 ) ) << std::endl;
//            std::cerr << 1000 << " " << 200 << " " << 200 << " " << log( normal_pdf(400, 200, 200 ) ) << std::endl;

//            for (auto hsp: high_scores){
//                auto score_ = std::get<0>(hsp);
//                auto s1_tmp = std::get<1>(hsp);
//                auto s2_tmp = std::get<2>(hsp);
//                std::cerr << "HSP SCORE: " << score_ << " " << s1_tmp.ref_start << " " << s2_tmp.ref_start << " " << s1_tmp.sw_score <<  " " << s2_tmp.sw_score << std::endl;
//            }
    int mapq1, mapq2;
    std::tie(mapq1, mapq2) = joint_mapq_from_high_scores(high_scores);

    auto best_aln_pair = high_scores[0];
    auto sam_aln1 = std::get<1>(best_aln_pair);
    auto sam_aln2 = std::get<2>(best_aln_pair);
    if (max_secondary == 0) {
//            get_MAPQ_aln(sam_aln1, sam_aln2);
        bool is_proper = is_proper_pair(sam_aln1, sam_aln2, mu, sigma);
        sam.add_pair(sam_aln1, sam_aln2, record1, record2, read1.rc, read2.rc,
                        mapq1, mapq2, is_proper, true);
    } else {
        int max_out = std::min(high_scores.size(), max_secondary);
        // remove eventual duplicates - comes from, e.g., adding individual best alignments above (if identical to joint best alignment)
        float s_max = std::get<0>(best_aln_pair);
        int prev_start_m1 = sam_aln1.ref_start;
        int prev_start_m2 = sam_aln2.ref_start;
        int prev_ref_id_m1 = sam_aln1.ref_id;
        int prev_ref_id_m2 = sam_aln2.ref_id;
        bool is_primary = true;
        for (int i = 0; i < max_out; ++i) {
            auto aln_pair = high_scores[i];
            sam_aln1 = std::get<1>(aln_pair);
            sam_aln2 = std::get<2>(aln_pair);
            float s_score = std::get<0>(aln_pair);
            if (i > 0) {
                is_primary = false;
                mapq1 = 255;
                mapq2 = 255;
                bool same_pos = (prev_start_m1 == sam_aln1.ref_start) && (prev_start_m2 == sam_aln2.ref_start);
                bool same_ref = (prev_ref_id_m1 == sam_aln1.ref_id) && (prev_ref_id_m2 == sam_aln2.ref_id);
                if ( same_pos && same_ref ){
                    continue;
                }
            }

            if (s_max - s_score < secondary_dropoff) {
                bool is_proper = is_proper_pair(sam_aln1, sam_aln2, mu, sigma);
                sam.add_pair(sam_aln1, sam_aln2, record1, record2, read1.rc, read2.rc,
                                mapq1, mapq2, is_proper, is_primary);
            } else {
                break;
            }

            prev_start_m1 = sam_aln1.ref_start;
            prev_start_m2 = sam_aln2.ref_start;
            prev_ref_id_m1 = sam_aln1.ref_id;
            prev_ref_id_m2 = sam_aln2.ref_id;
        }
    }
}

inline void get_best_map_location(std::vector<Nam> &nams1, std::vector<Nam> &nams2, i_dist_est &isize_est, Nam &best_nam1,  Nam &best_nam2 ) {
    std::vector<std::tuple<int,Nam,Nam>> joint_NAM_scores = get_best_scoring_NAM_locations(nams1, nams2, isize_est.mu, isize_est.sigma);
    Nam n1_joint_max, n2_joint_max, n1_indiv_max, n2_indiv_max;
    float score_joint = 0;
    float score_indiv = 0;
    best_nam1.ref_s = -1; //Unmapped until proven mapped
    best_nam2.ref_s = -1; //Unmapped until proven mapped

    if (joint_NAM_scores.empty()) {
        return;
    }
    // get best joint score
    for (auto &t : joint_NAM_scores) { // already sorted by descending score
        auto n1 = std::get<1>(t);
        auto n2 = std::get<2>(t);
        if ((n1.ref_s >=0) && (n2.ref_s >=0) ){ // Valid pair
            score_joint =  n1.score + n2.score;
            n1_joint_max = n1;
            n2_joint_max = n2;
            break;
        }
    }

    // get individual best scores
    if (!nams1.empty()) {
        auto n1_indiv_max = nams1[0];
        score_indiv += n1_indiv_max.score - (n1_indiv_max.score/2.0); //Penalty for being mapped individually
        best_nam1 = n1_indiv_max;
    }
    if (!nams2.empty()) {
        auto n2_indiv_max = nams2[0];
        score_indiv += n2_indiv_max.score - (n2_indiv_max.score/2.0); //Penalty for being mapped individually
        best_nam2 = n2_indiv_max;
    }
    if ( score_joint > score_indiv ){ // joint score is better than individual
        best_nam1 = n1_joint_max;
        best_nam2 = n2_joint_max;
    }

    if (isize_est.sample_size < 400 && score_joint > score_indiv) {
        isize_est.update(std::abs(n1_joint_max.ref_s - n2_joint_max.ref_s));
    }
}

/* Add a new observation */
void i_dist_est::update(int dist)
{
    if (dist >= 2000) {
        return;
    }
    const float e = dist - mu;
    mu += e / sample_size; // (1.0/(sample_size +1.0)) * (sample_size*mu + d);
    SSE += e * (dist - mu);
    if (sample_size > 1) {
        //d < 1000 ? ((sample_size +1.0)/sample_size) * ( (V*sample_size/(sample_size +1)) + ((mu-d)*(mu-d))/sample_size ) : V;
        V = SSE / (sample_size - 1.0);
    } else {
        V = SSE;
    }
    sigma = std::sqrt(V);
    sample_size = sample_size + 1.0;
    if (mu < 0) {
        std::cerr << "mu negative, mu: " << mu << " sigma: " << sigma << " SSE: " << SSE << " sample size: " << sample_size << std::endl;
    }
    if (SSE < 0) {
        std::cerr << "SSE negative, mu: " << mu << " sigma: " << sigma << " SSE: " << SSE << " sample size: " << sample_size << std::endl;
    }
}


void align_PE_read(
    const KSeq &record1,
    const KSeq &record2,
    Sam& sam,
    std::string& outstring,
    AlignmentStatistics &statistics,
    i_dist_est &isize_est,
    const Aligner &aligner,
    const mapping_params &map_param,
    const IndexParameters& index_parameters,
    const References& references,
    const StrobemerIndex& index
) {
    Timer strobe_timer;
    auto query_randstrobes1 = randstrobes_query(
        index_parameters.k, index_parameters.w_min, index_parameters.w_max, record1.seq, index_parameters.s, index_parameters.t_syncmer,
        index_parameters.q, index_parameters.max_dist);
    auto query_randstrobes2 = randstrobes_query(
        index_parameters.k, index_parameters.w_min, index_parameters.w_max, record2.seq, index_parameters.s, index_parameters.t_syncmer,
        index_parameters.q, index_parameters.max_dist);
    statistics.tot_construct_strobemers += strobe_timer.duration();

    // Find NAMs
    Timer nam_timer;
    auto [nonrepetitive_fraction1, nams1] = find_nams(query_randstrobes1, index);
    auto [nonrepetitive_fraction2, nams2] = find_nams(query_randstrobes2, index);
    statistics.tot_find_nams += nam_timer.duration();

    if (map_param.R > 1) {
        Timer rescue_timer;
        if (nams1.empty() || nonrepetitive_fraction1 < 0.7) {
            statistics.tried_rescue += 1;
            nams1 = find_nams_rescue(query_randstrobes1, index, map_param.rescue_cutoff);
        }

        if (nams2.empty() || nonrepetitive_fraction2 < 0.7) {
            statistics.tried_rescue += 1;
            nams2 = find_nams_rescue(query_randstrobes2, index, map_param.rescue_cutoff);
        }
        statistics.tot_time_rescue += rescue_timer.duration();
    }

    Timer nam_sort_timer;
    std::sort(nams1.begin(), nams1.end(), score);
    std::sort(nams2.begin(), nams2.end(), score);
    statistics.tot_sort_nams += nam_timer.duration();

    Timer extend_timer;
    if (!map_param.is_sam_out) {
        Nam nam_read1;
        Nam nam_read2;
        get_best_map_location(nams1, nams2, isize_est,
                              nam_read1,
                              nam_read2);
        output_hits_paf_PE(outstring, nam_read1, record1.name,
                           references,
                           index_parameters.k,
                           record1.seq.length());
        output_hits_paf_PE(outstring, nam_read2, record2.name,
                           references,
                           index_parameters.k,
                           record2.seq.length());
    } else {
        align_PE(aligner, sam, nams1, nams2, record1,
                 record2,
                 index_parameters.k,
                 references, statistics,
                 map_param.dropoff_threshold, isize_est, map_param.maxTries, map_param.max_secondary);
    }
    statistics.tot_extend += extend_timer.duration();
}


void align_SE_read(
    const KSeq &record,
    Sam& sam,
    std::string &outstring,
    AlignmentStatistics &statistics,
    const Aligner &aligner,
    const mapping_params &map_param,
    const IndexParameters& index_parameters,
    const References& references,
    const StrobemerIndex& index
) {
    Timer strobe_timer;
    auto query_randstrobes = randstrobes_query(index_parameters.k, index_parameters.w_min, index_parameters.w_max, record.seq, index_parameters.s, index_parameters.t_syncmer, index_parameters.q, index_parameters.max_dist);
    statistics.tot_construct_strobemers += strobe_timer.duration();

    // Find NAMs
    Timer nam_timer;
    auto [nonrepetitive_fraction, nams] = find_nams(query_randstrobes, index);
    statistics.tot_find_nams += nam_timer.duration();

    if (map_param.R > 1) {
        Timer rescue_timer;
        if (nams.empty() || nonrepetitive_fraction < 0.7) {
            statistics.tried_rescue += 1;
            nams = find_nams_rescue(query_randstrobes, index, map_param.rescue_cutoff);
        }
        statistics.tot_time_rescue += rescue_timer.duration();
    }

    Timer nam_sort_timer;
    std::sort(nams.begin(), nams.end(), score);
    statistics.tot_sort_nams += nam_sort_timer.duration();

    Timer extend_timer;
    if (!map_param.is_sam_out) {
        output_hits_paf(outstring, nams, record.name, references, index_parameters.k,
                        record.seq.length());
    } else {
        align_SE_secondary_hits(
            aligner, sam, nams, record, index_parameters.k,
            references, statistics, map_param.dropoff_threshold, map_param.maxTries,
            map_param.max_secondary + 1
        );
    }
    statistics.tot_extend += extend_timer.duration();
}
