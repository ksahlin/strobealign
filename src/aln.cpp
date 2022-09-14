#include "aln.hpp"

#include <iostream>
#include <math.h>
#include <sstream>
#include "ssw_cpp.h"
#include "sam.hpp"

using namespace klibpp;

static inline bool score(const nam &a, const nam &b)
{
    return ( a.score > b.score );
}

inline aln_info ssw_align(const std::string &ref, const std::string &query, int match_score, int mismatch_penalty, int gap_opening_penalty, int gap_extending_penalty) {

    aln_info aln;
    int32_t maskLen = strlen(query.c_str())/2;
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

    StripedSmithWaterman::Aligner aligner(match_score, mismatch_penalty, gap_opening_penalty, gap_extending_penalty);
//    StripedSmithWaterman::Aligner aligner;
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment alignment_ssw;
    aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment_ssw, maskLen, 1);
    // Have to give up this optimization untill the 'Command terminated abnormally' bug is fixed in ssw library
//     if (read_len*match_score < 255){
//         std::cerr << "Here: "  << read_len*match_score << " " << ref.length() << std::endl;
//         try
//         {
//             aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment_ssw, maskLen, 0);
//         }
//         catch (...)
//         {
//             aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment_ssw, maskLen, 1);
//         }
//
//     } else {
//            aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment_ssw, maskLen, 1);
//     }
//    std::cerr << passed << std::endl;
//    if(!passed){
//        std::cerr << "Failed" << std::endl;
//        std::cerr << "read: " << query << std::endl;
//        std::cerr << "ref: "  << ref << std::endl;
//    }


//    std::cerr << "===== SSW result =====" << std::endl;
//    std::cerr << "Best Smith-Waterman score:\t" << alignment_ssw.sw_score << std::endl
//         << "Next-best Smith-Waterman score:\t" << alignment_ssw.sw_score_next_best << std::endl
//         << "Reference start:\t" << alignment_ssw.ref_begin << std::endl
//         << "Reference end:\t" << alignment_ssw.ref_end << std::endl
//         << "Query start:\t" << alignment_ssw.query_begin << std::endl
//         << "Query end:\t" << alignment_ssw.query_end << std::endl
//         << "Next-best reference end:\t" << alignment_ssw.ref_end_next_best << std::endl
//         << "Number of mismatches:\t" << alignment_ssw.mismatches << std::endl
//         << "Cigar: " << alignment_ssw.cigar_string << std::endl;

    aln.global_ed = alignment_ssw.global_ed;
    aln.ed = alignment_ssw.mismatches;
    aln.ref_offset = alignment_ssw.ref_begin;
    aln.cigar = alignment_ssw.cigar_string;
    aln.sw_score = alignment_ssw.sw_score;
    aln.length = alignment_ssw.ref_end - alignment_ssw.ref_begin;
    return aln;
}



//static inline bool sort_hits(const hit &a, const hit &b)
//{
//    // first sort on ref ID, then on reference starts, then on query starts
//    return (a.ref_id < b.ref_id) ||
//           ( (a.ref_id == b.ref_id) && (a.ref_s < b.ref_s) ) ||
//           ((a.ref_id == b.ref_id) && (a.ref_s == b.ref_s) && (a.query_s < b.query_s )) ;
//}
//static inline std::vector<nam> find_nams_alt(mers_vector_read &query_mers, mers_vector_reduced &ref_mers, kmer_lookup &mers_index, int k, std::vector<std::string> &ref_seqs, std::string &read, unsigned int hit_upper_window_lim, unsigned int filter_cutoff ) {
////    std::cerr << "ENTER FIND NAMS " <<  std::endl;
//    std::vector<hit> all_hits;
//    for (auto &q : query_mers)
////    for (size_t i = 0; i < query_mers.size(); ++i)
//    {
////        std::cerr << "Q " << h.query_s << " " << h.query_e << " read length:" << read_length << std::endl;
//        uint64_t mer_hashv = std::get<0>(q);
//        if (mers_index.find(mer_hashv) != mers_index.end()) { //  In  index
//            hit h;
//            h.query_s = std::get<2>(q);
//            h.query_e = std::get<3>(q) + k; // h.query_s + read_length/2;
//            h.is_rc = std::get<4>(q);
//            std::tuple<uint64_t, unsigned int> mer;
//            mer = mers_index[mer_hashv];
//            unsigned int offset = std::get<0>(mer);
//            unsigned int count = std::get<1>(mer);
//
//            for (size_t j = offset; j < offset + count; ++j) {
//                auto r = ref_mers[j];
//                unsigned int ref_id = std::get<0>(r);
//                unsigned int ref_s = std::get<1>(r);
//                unsigned int ref_e = std::get<2>(r) + k; //ref_s + read_length/2;
//
//                h.ref_id = ref_id;
//                h.ref_s = ref_s;
//                h.ref_e = ref_e;
//                all_hits.push_back(h);
//            }
//
//        }
//    }
//
//    std::sort(all_hits.begin(), all_hits.end(), sort_hits);
//
//    std::vector<nam> final_nams; // [ref_id] -> vector(struct nam)
//    nam o;
//    if (all_hits.size() == 0){
//        return final_nams;
//    }
//    else{
//        hit h = all_hits[0];
//        o.ref_id = h.ref_id;
//        o.query_s = h.query_s;
//        o.query_e = h.query_e;
//        o.ref_s = h.ref_s;
//        o.ref_e = h.ref_e;
////        o.previous_query_start = h.query_s;
////        o.previous_ref_start = h.ref_s;
//        o.query_prev_hit_startpos = h.query_s;
//        o.ref_prev_hit_startpos = h.ref_s;
//        o.n_hits = 1;
//        o.is_rc = h.is_rc;
//    }
//
//    hit h;
//    for(size_t i = 1; i < all_hits.size(); ++i) // all but first element
//    {
//        h = all_hits[i];
//        if ( (o.ref_id == h.ref_id) && ( o.is_rc == h.is_rc) && ( o.query_prev_hit_startpos < h.query_s) && (h.query_s <= o.query_e ) && ( o.ref_prev_hit_startpos < h.ref_s) && (h.ref_s <= o.ref_e)){
//            if ( (h.query_e > o.query_e) && (h.ref_e > o.ref_e) ) {
//                o.query_e = h.query_e;
//                o.ref_e = h.ref_e;
////                o.previous_query_start = h.query_s;
////                o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
//                o.query_prev_hit_startpos = h.query_s; // log the last strobemer hit in case of outputting paf
//                o.ref_prev_hit_startpos = h.ref_s; // log the last strobemer hit in case of outputting paf
//                o.n_hits ++;
//            }
//            else if ((h.query_e <= o.query_e) && (h.ref_e <= o.ref_e)) {
////                o.previous_query_start = h.query_s;
////                o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
//                o.query_prev_hit_startpos = h.query_s; // log the last strobemer hit in case of outputting paf
//                o.ref_prev_hit_startpos = h.ref_s; // log the last strobemer hit in case of outputting paf
//                o.n_hits ++;
//            }
//        } else {
//            final_nams.push_back(o);
//            o.ref_id = h.ref_id;
//            o.query_s = h.query_s;
//            o.query_e = h.query_e;
//            o.ref_s = h.ref_s;
//            o.ref_e = h.ref_e;
////            o.previous_query_start = h.query_s;
////            o.previous_ref_start = h.ref_s;
//            o.query_prev_hit_startpos = h.query_s;
//            o.ref_prev_hit_startpos = h.ref_s;
//            o.n_hits = 1;
//            o.is_rc = h.is_rc;
//        }
//    }
//
//
//    final_nams.push_back(o);
//
//
////        for (auto &n : final_nams){
////        std::cerr << "NAM ALT: " << n.ref_id << ": (" << n.n_hits << ", " << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e  << ")" << std::endl;
////    }
////    std::cerr << " " << std::endl;
//
//    return final_nams;
//}


static inline bool sort_hits(const hit &a, const hit &b)
{
    // first sort on query starts, then on reference starts
    return (a.query_s < b.query_s) || ( (a.query_s == b.query_s) && (a.ref_s < b.ref_s) );
}

static inline void find_nams_rescue(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> hits_fw, std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> hits_rc, std::vector<nam> &final_nams, robin_hood::unordered_map< unsigned int, std::vector<hit>> &hits_per_ref, const mers_vector_read &query_mers, const mers_vector &ref_mers, kmer_lookup &mers_index, int k, const std::vector<std::string> &ref_seqs, const std::string &read, unsigned int filter_cutoff ){
//    std::pair<float,int> info (0,0); // (nr_nonrepetitive_hits/total_hits, max_nam_n_hits)
    int nr_good_hits = 0, total_hits = 0;
    bool is_rc = true, no_rep_fw = true, no_rep_rc = true;
//    std::pair<int, int> repeat_fw(0,0), repeat_rc(0,0);
//    std::vector<std::pair<int, int>> repetitive_fw, repetitive_rc;
    for (auto &q : query_mers)
    {
        auto mer_hashv = std::get<0>(q);
        if (mers_index.find(mer_hashv) != mers_index.end()){ //  In  index
            total_hits ++;
            auto ref_hit = mers_index[mer_hashv];
            auto offset = std::get<0>(ref_hit);
            auto count = std::get<1>(ref_hit);
            auto query_s = std::get<2>(q);
            auto query_e = query_s + std::get<3>(q) + k;
            is_rc = std::get<4>(q);
            if (is_rc){
                std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool> s(count, offset, query_s, query_e, is_rc);
                hits_rc.push_back(s);
//                if (count > filter_cutoff){
//                    if (no_rep_rc){ //initialize
//                        repeat_rc.first = query_s;
//                        repeat_rc.second = query_e;
//                        no_rep_rc = false;
//                    }
//                    else if (query_s >= repeat_rc.second){
//                        repetitive_rc.push_back(repeat_rc);
//                        repeat_rc.first = query_s;
//                        repeat_rc.second = query_e;
//                    } else{
//                        repeat_rc.second = std::max(repeat_rc.second, query_e);
//                    }
//                } else{
//                    nr_good_hits ++;
//                }
            } else{
                std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool> s(count, offset, query_s, query_e, is_rc);
                hits_fw.push_back(s);
//                if (count > filter_cutoff){
//                    if (no_rep_fw){ //initialize
//                        repeat_fw.first = query_s;
//                        repeat_fw.second = query_e;
//                        no_rep_fw = false;
//                    }
//                    else if (query_s >= repeat_fw.second ){
//                        repetitive_fw.push_back(repeat_fw);
//                        repeat_fw.first = query_s;
//                        repeat_fw.second = query_e;
//                    } else{
//                        repeat_fw.second = std::max(repeat_fw.second, query_e);
//                    }
//                } else{
//                    nr_good_hits ++;
//                }
            }


        }
    }
//    if (!no_rep_fw) {
//        repetitive_fw.push_back(repeat_fw);
//    }
//    if (!no_rep_rc) {
//        repetitive_rc.push_back(repeat_rc);
//    }
    std::sort(hits_fw.begin(), hits_fw.end());
    std::sort(hits_rc.begin(), hits_rc.end());

//    for (auto &rf : repetitive_fw){
//        std::cerr << "REPEAT MASKED FW: (" << rf.first << " " << rf.second << ") " << std::endl;
//    }
//    for (auto &rc : repetitive_rc){
//        std::cerr << "REPEAT MASKED RC: (" << rc.first << " " << rc.second << ") " << std::endl;
//    }

    hit h;
    int cnt = 0;
    for (auto &q : hits_fw)
    {
//        std::cerr << "Q " << h.query_s << " " << h.query_e << " read length:" << read_length << std::endl;
        auto count = std::get<0>(q);
        auto offset = std::get<1>(q);
        h.query_s = std::get<2>(q);
        h.query_e = std::get<3>(q); // h.query_s + read_length/2;
        h.is_rc = std::get<4>(q);


        if ( ((count <= filter_cutoff) || (cnt < 5)) && (count <= 1000) ){
//            std::cerr << "Found FORWARD: " << count << ", q_start: " <<  h.query_s << ", q_end: " << h.query_e << std::endl;
            int min_diff = 1000;
//            int ref_d;
//            for(size_t j = offset; j < offset+count; ++j) {
//                auto r = ref_mers[j];
//                ref_d = std::get<2>(r) + k - std::get<1>(r); //I changed this code from std::get<3>(r) + k - std::get<2>(r); - but it is probably old, 3 should not be possible since we only had 3 members in the tuple
//                int diff = (h.query_e - h.query_s) - ref_d > 0 ? (h.query_e - h.query_s) -  ref_d : ref_d - (h.query_e - h.query_s);
//                if (diff <= min_diff ){
//                    min_diff = diff;
//                }
//            }

            for(size_t j = offset; j < offset+count; ++j)
            {
                auto r = ref_mers[j];
                h.ref_s = std::get<0>(r);
                auto p = std::get<1>(r);
                int bit_alloc = 8;
                int r_id = (p >> bit_alloc);
                int mask=(1<<bit_alloc) - 1;
                int offset = (p & mask);
                h.ref_e = h.ref_s + offset + k;
//                h.count = count;
//                hits_per_ref[std::get<0>(r)].push_back(h);

                int diff = std::abs((h.query_e - h.query_s) - (h.ref_e - h.ref_s));
                if (diff <= min_diff ){
                    hits_per_ref[r_id].push_back(h);
                    min_diff = diff;
                }
            }
            cnt ++;
        }
        else{
            break;
//            std::cerr << "Found repetitive count FORWARD: " << count << ", q_start: " <<  h.query_s << ", q_end: " << h.query_e << std::endl;

        }

    }

    cnt = 0;
    for (auto &q : hits_rc)
    {
//        std::cerr << "Q " << h.query_s << " " << h.query_e << " read length:" << read_length << std::endl;
        auto count = std::get<0>(q);
        auto offset = std::get<1>(q);
        h.query_s = std::get<2>(q);
        h.query_e = std::get<3>(q); // h.query_s + read_length/2;
        h.is_rc = std::get<4>(q);

        if ( ((count <= filter_cutoff) || (cnt < 5)) && (count <= 1000) ){
//            std::cerr << "Found REVERSE: " << count << ", q_start: " <<  h.query_s << ", q_end: " << h.query_e << std::endl;
            int min_diff = 1000;
//            int ref_d;
//            for(size_t j = offset; j < offset+count; ++j) {
//                auto r = ref_mers[j];
//                ref_d = std::get<2>(r) + k - std::get<1>(r); //same problem here, I changed from 3 to 2 etc., but it doesn't seem right
//                int diff = (h.query_e - h.query_s) - ref_d > 0 ? (h.query_e - h.query_s) -  ref_d : ref_d - (h.query_e - h.query_s);
//                if (diff <= min_diff ){
//                    min_diff = diff;
//                }
//            }

            for(size_t j = offset; j < offset+count; ++j)
            {
                auto r = ref_mers[j];
                h.ref_s = std::get<0>(r);
                auto p = std::get<1>(r);
                int bit_alloc = 8;
                int r_id = (p >> bit_alloc);
                int mask=(1<<bit_alloc) - 1;
                int offset = (p & mask);
                h.ref_e = h.ref_s + offset + k;
//                h.count = count;
//                hits_per_ref[std::get<1>(r)].push_back(h);
                int diff = std::abs((h.query_e - h.query_s) - (h.ref_e - h.ref_s));
                if (diff <= min_diff ){
                    hits_per_ref[r_id].push_back(h);
                    min_diff = diff;
                }
            }
            cnt ++;
        }
        else{
            break;
//            std::cerr << "Found repetitive count REVERSE: " << count << ", q_start: " <<  h.query_s << ", q_end: " << h.query_e << std::endl;

        }

    }

//    std::cerr << "NUMBER OF HITS GENERATED: " << hit_count_all << std::endl;
//    info.first = total_hits > 0 ? ((float) nr_good_hits) / ((float) total_hits) : 1.0;
    int max_nam_n_hits = 0;
    std::vector<nam> open_nams;
    int nam_id_cnt = 0;
//    std::vector<nam> final_nams; // [ref_id] -> vector(struct nam)

    for (auto &it : hits_per_ref)
    {
        auto ref_id = it.first;
        std::vector<hit> hits = it.second;
        std::sort(hits.begin(), hits.end(), sort_hits);
        open_nams = std::vector<nam> (); // Initialize vector
        unsigned int prev_q_start = 0;
        for (auto &h : hits){
            bool is_added = false;
//            std::cerr << "HIT " << h.is_rc << " " << h.query_s <<  ", " << h.query_e << ", " << h.ref_s <<  ", " << h.ref_e << std::endl;
//            bool local_repeat_worse_fit = false;
            for (auto & o : open_nams) {

                // Extend NAM
                if (( o.is_rc == h.is_rc) && (o.query_prev_hit_startpos < h.query_s) && (h.query_s <= o.query_e ) && (o.ref_prev_hit_startpos < h.ref_s) && (h.ref_s <= o.ref_e) ){
                    if ( (h.query_e > o.query_e) && (h.ref_e > o.ref_e) ) {
                        o.query_e = h.query_e;
                        o.ref_e = h.ref_e;
//                        o.previous_query_start = h.query_s;
//                        o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
                        o.query_prev_hit_startpos = h.query_s; // log the last strobemer hit in case of outputting paf
                        o.ref_prev_hit_startpos = h.ref_s; // log the last strobemer hit in case of outputting paf
                        o.n_hits ++;
//                        o.score += (float)1/ (float)h.count;
                        is_added = true;
                        break;
                    }
                    else if ((h.query_e <= o.query_e) && (h.ref_e <= o.ref_e)) {
//                        o.previous_query_start = h.query_s;
//                        o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
                        o.query_prev_hit_startpos = h.query_s; // log the last strobemer hit in case of outputting paf
                        o.ref_prev_hit_startpos = h.ref_s; // log the last strobemer hit in case of outputting paf
                        o.n_hits ++;
//                        o.score += (float)1/ (float)h.count;
                        is_added = true;
                        break;
                    }
//                    else if ( (o.query_e - o.query_s) - (o.ref_e - o.ref_s) > (h.query_e - h.query_s) - (h.ref_e - h.ref_s)  ){
//                        local_repeat_worse_fit = true;
//                    }

                }


            }
//            if (local_repeat_worse_fit){
//                continue;
//            }
            // Add the hit to open matches
            if (!is_added){
                nam n;
                n.nam_id = nam_id_cnt;
                nam_id_cnt ++;
                n.query_s = h.query_s;
                n.query_e = h.query_e;
                n.ref_s = h.ref_s;
                n.ref_e = h.ref_e;
                n.ref_id = ref_id;
//                n.previous_query_start = h.query_s;
//                n.previous_ref_start = h.ref_s;
                n.query_prev_hit_startpos = h.query_s;
                n.ref_prev_hit_startpos = h.ref_s;
                n.n_hits = 1;
                n.is_rc = h.is_rc;
//                n.score += (float)1/ (float)h.count;
                open_nams.push_back(n);
            }


            // Only filter if we have advanced at least k nucleotides
            if (h.query_s > prev_q_start + k) {

                // Output all NAMs from open_matches to final_nams that the current hit have passed
                for (auto &n : open_nams) {
                    if (n.query_e < h.query_s) {
                        int n_max_span = std::max(n.query_e - n.query_s, n.ref_e - n.ref_s);
                        int n_min_span = std::max(n.query_e - n.query_s, n.ref_e - n.ref_s);
                        float n_score;
                        n_score = ( 2*n_min_span -  n_max_span) > 0 ? (float) (n.n_hits * ( 2*n_min_span -  n_max_span) ) : 1;   // this is really just n_hits * ( min_span - (offset_in_span) ) );
//                        n_score = n.n_hits * (n.query_e - n.query_s);
                        n.score = n_score;
                        final_nams.push_back(n);
                        max_nam_n_hits = std::max(n.n_hits, max_nam_n_hits);
                    }
                }

                // Remove all NAMs from open_matches that the current hit have passed
                unsigned int c = h.query_s;
                auto predicate = [c](decltype(open_nams)::value_type const &nam) { return nam.query_e < c; };
                open_nams.erase(std::remove_if(open_nams.begin(), open_nams.end(), predicate), open_nams.end());
                prev_q_start = h.query_s;
            }


        }

        // Add all current open_matches to final NAMs
        for (auto &n : open_nams){
            int n_max_span = std::max(n.query_e - n.query_s, n.ref_e - n.ref_s);
            int n_min_span = std::min(n.query_e - n.query_s, n.ref_e - n.ref_s);
            float n_score;
            n_score = ( 2*n_min_span -  n_max_span) > 0 ? (float) (n.n_hits * ( 2*n_min_span -  n_max_span) ) : 1;   // this is really just n_hits * ( min_span - (offset_in_span) ) );
//            n_score = n.n_hits * (n.query_e - n.query_s);
            n.score = n_score;
            final_nams.push_back(n);
            max_nam_n_hits = std::max(n.n_hits, max_nam_n_hits);
        }
    }

//    for (auto &n : final_nams){
//        std::cerr << "RESCUE NAM: " << n.ref_id << ": (" << n.score << ", " << n.n_hits << ", " << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e  << ")" << " " <<  n.is_rc << std::endl;
//    }
//    info.second = max_nam_n_hits;
    return ;

}



static inline std::pair<float,int> find_nams(std::vector<nam> &final_nams, robin_hood::unordered_map< unsigned int, std::vector<hit>> &hits_per_ref, mers_vector_read &query_mers, mers_vector &ref_mers, kmer_lookup &mers_index, int k, std::vector<std::string> &ref_seqs, std::string &read, unsigned int filter_cutoff ){
//    std::cerr << "ENTER FIND NAMS " <<  std::endl;
//    robin_hood::unordered_map< unsigned int, std::vector<hit>> hits_per_ref; // [ref_id] -> vector( struct hit)
//    std::vector<std::vector<hit>> hits_per_ref(10);
//    int read_length = read.length();
//    std::cerr << " "  <<  std::endl;
    std::pair<float,int> info (0.0f,0); // (nr_nonrepetitive_hits/total_hits, max_nam_n_hits)
    int nr_good_hits = 0, total_hits = 0;
    hit h;
    for (auto &q : query_mers)
//    for (size_t i = 0; i < query_mers.size(); ++i)
    {
//        std::cerr << "Q " << h.query_s << " " << h.query_e << " read length:" << read_length << std::endl;
        auto mer_hashv = std::get<0>(q);
        if (mers_index.find(mer_hashv) != mers_index.end()){ //  In  index
            total_hits ++;
            h.query_s = std::get<2>(q);
            h.query_e = h.query_s + std::get<3>(q) + k; // h.query_s + read_length/2;
            h.is_rc = std::get<4>(q);
            auto mer = mers_index[mer_hashv];
            auto offset = std::get<0>(mer);
            auto count = std::get<1>(mer);
//            if (count == 1){
//                auto r = ref_mers[offset];
//                unsigned int ref_id = std::get<0>(r); //The indexes in this code are not fixed after removal of the 64-bit hash
//                unsigned int ref_s = std::get<1>(r);
//                unsigned int ref_e = std::get<2>(r) + k; //ref_s + read_length/2;
//
//                h.ref_s = ref_s;
//                h.ref_e = ref_e;
//                hits_per_ref[ref_id].push_back(h);
//                h.hit_count = count;
//                hit_count_all ++;
//            } else
            if (count <= filter_cutoff){
                nr_good_hits ++;
//                bool start_log = false;
                int min_diff = 100000;
//                int tries = 0;
//                int ref_d;
//                for(size_t j = offset; j < offset+count; ++j) {
//                    auto r = ref_mers[j];
//                    ref_d = std::get<3>(r) + k - std::get<2>(r);//The indexes in this code are not fixed after removal of the 64-bit hash
//                    int diff = (h.query_e - h.query_s) - ref_d > 0 ? (h.query_e - h.query_s) -  ref_d : ref_d - (h.query_e - h.query_s);
//                    if (diff <= min_diff ){
//                        min_diff = diff;
//                    }
//                }
//                std::cerr << "Found good count: " << count << ", q_start: " <<  h.query_s << ", q_end: " << h.query_e << std::endl;
                for(size_t j = offset; j < offset+count; ++j)
//                for(auto r = begin(ref_mers) + offset; r != begin(ref_mers) + offset + count; ++r)
                {
                    auto r = ref_mers[j];
//                    unsigned int  ref_id,ref_s,ref_e; std::tie(ref_id,ref_s,ref_e) = r;
//                    unsigned int ref_id = std::get<0>(r);//The indexes in this code are not fixed after removal of the 64-bit hash
//                    unsigned int ref_s = std::get<1>(r);
//                    unsigned int ref_e = std::get<2>(r) + k; //ref_s + read_length/2;
                    h.ref_s = std::get<0>(r);
                    auto p = std::get<1>(r);
                    int bit_alloc = 8;
                    int r_id = (p >> bit_alloc);
                    int mask=(1<<bit_alloc) - 1;
                    int offset = (p & mask);
                    h.ref_e = h.ref_s + offset + k;
//                    h.count = count;
//                    hits_per_ref[std::get<1>(r)].push_back(h);
//                    hits_per_ref[std::get<0>(r)].push_back(h);


//                    h.ref_s = ref_s;
//                    h.ref_e = ref_e;
//                    hits_per_ref[ref_id].push_back(h);
                    int diff = std::abs((h.query_e - h.query_s) - (h.ref_e - h.ref_s));
//                    if ((diff > 0) || start_log ){
//                        std::cerr << "Found: " <<  count << " " << diff << " " << h.query_e - h.query_s << " " <<  (h.ref_e - h.ref_s) << std::endl;
//                        start_log = true;
//                    }
                    if (diff <= min_diff ){
                        hits_per_ref[r_id].push_back(h);
                        min_diff = diff;
//                        std::cerr << "Found: query: " <<  h.query_s << " " << h.query_e << " ref: " <<  h.ref_s << " " << h.ref_e << " " << h.is_rc << " diff " << diff << std::endl;
//                        tries ++;
                    }
//                    if (tries > filter_cutoff){
//                        break;
//                    }
//                    h.hit_count = count;
//                    if (count > 1){
//                        int diff = (h.query_e - h.query_s) - (h.ref_e - h.ref_s);
//                        std::cerr << "Found: " <<  h.query_s << " " << h.query_e << " ref: " <<  h.ref_s << " " << h.ref_e << " " << h.is_rc << " diff " << diff << std::endl;
//                    }
//                    hit_count_all ++;

                }
            }

//            else{
//                std::cerr << "Found repetitive count: " << count << ", q_start: " <<  h.query_s << ", q_end: " << h.query_e << std::endl;
//
//            }

        }
    }

//    std::cerr << "NUMBER OF HITS GENERATED: " << hit_count_all << std::endl;
    info.first = total_hits > 0 ? ((float) nr_good_hits) / ((float) total_hits) : 1.0;
    int max_nam_n_hits = 0;
    int nam_id_cnt = 0;
    std::vector<nam> open_nams;
//    std::vector<nam> final_nams; // [ref_id] -> vector(struct nam)

    for (auto &it : hits_per_ref)
    {
        auto ref_id = it.first;
        std::vector<hit> hits = it.second;

//    for(size_t i = 0; i < hits_per_ref.size(); ++i){
//        unsigned int ref_id = i;
//        auto hits = hits_per_ref[i];
        open_nams = std::vector<nam> (); // Initialize vector
        unsigned int prev_q_start = 0;
        for (auto &h : hits){
            bool is_added = false;
//            std::cerr << "HIT " << h.is_rc << " " << h.query_s <<  ", " << h.query_e << ", " << h.ref_s <<  ", " << h.ref_e << std::endl;
//            bool local_repeat_worse_fit = false;
            for (auto & o : open_nams) {

                // Extend NAM
                if (( o.is_rc == h.is_rc) && (o.query_prev_hit_startpos < h.query_s) && (h.query_s <= o.query_e ) && (o.ref_prev_hit_startpos < h.ref_s) && (h.ref_s <= o.ref_e) ){
                    if ( (h.query_e > o.query_e) && (h.ref_e > o.ref_e) ) {
                        o.query_e = h.query_e;
                        o.ref_e = h.ref_e;
//                        o.previous_query_start = h.query_s;
//                        o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
                        o.query_prev_hit_startpos = h.query_s; // log the last strobemer hit in case of outputting paf
                        o.ref_prev_hit_startpos = h.ref_s; // log the last strobemer hit in case of outputting paf
                        o.n_hits ++;
//                        o.score += (float)1/ (float)h.count;
                        is_added = true;
                        break;
                    }
                    else if ((h.query_e <= o.query_e) && (h.ref_e <= o.ref_e)) {
//                        o.previous_query_start = h.query_s;
//                        o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
                        o.query_prev_hit_startpos = h.query_s; // log the last strobemer hit in case of outputting paf
                        o.ref_prev_hit_startpos = h.ref_s; // log the last strobemer hit in case of outputting paf
                        o.n_hits ++;
//                        o.score += (float)1/ (float)h.count;
                        is_added = true;
                        break;
                    }
//                    else if ( (o.query_e - o.query_s) - (o.ref_e - o.ref_s) > (h.query_e - h.query_s) - (h.ref_e - h.ref_s)  ){
//                        local_repeat_worse_fit = true;
//                    }

                }

//                // CHECK IF FALSE REVERSE HITS FROM SYM HASHES
//                if (( o.is_rc == h.is_rc) && (o.query_prev_hit_startpos < h.query_s) && (h.query_s <= o.query_e ) && (o.ref_prev_hit_startpos <= h.ref_e) && (h.ref_e < o.ref_e) ){
//                    if ( (h.query_e > o.query_e) && (h.ref_s < o.ref_s) ) {
//                        o.query_e = h.query_e;
//                        o.ref_s = h.ref_s;
////                        o.previous_query_start = h.query_s;
////                        o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
//                        o.query_prev_hit_startpos = h.query_s; // log the last strobemer hit in case of outputting paf
//                        o.ref_prev_hit_startpos = h.ref_s; // log the last strobemer hit in case of outputting paf
//                        o.n_hits ++;
////                        o.score += (float)1/ (float)h.count;
//                        is_added = true;
//                        break;
//                    }
//                    else if ((h.query_e <= o.query_e) && (h.ref_s >= o.ref_s)) {
//                        o.query_prev_hit_startpos = h.query_s; // log the last strobemer hit in case of outputting paf
//                        o.ref_prev_hit_startpos = h.ref_s; // log the last strobemer hit in case of outputting paf
//                        o.n_hits ++;
//                        is_added = true;
//                        break;
//                    }
//                }

            }
//            if (local_repeat_worse_fit){
//                continue;
//            }
            // Add the hit to open matches
            if (!is_added){
                nam n;
                n.nam_id = nam_id_cnt;
                nam_id_cnt ++;
                n.query_s = h.query_s;
                n.query_e = h.query_e;
                n.ref_s = h.ref_s;
                n.ref_e = h.ref_e;
                n.ref_id = ref_id;
//                n.previous_query_start = h.query_s;
//                n.previous_ref_start = h.ref_s;
                n.query_prev_hit_startpos = h.query_s;
                n.ref_prev_hit_startpos = h.ref_s;
                n.n_hits = 1;
                n.is_rc = h.is_rc;
//                n.score += (float)1 / (float)h.count;
                open_nams.push_back(n);
            }


            // Only filter if we have advanced at least k nucleotides
            if (h.query_s > prev_q_start + k) {

                // Output all NAMs from open_matches to final_nams that the current hit have passed
                for (auto &n : open_nams) {
                    if (n.query_e < h.query_s) {
                        int n_max_span = std::max(n.query_e - n.query_s, n.ref_e - n.ref_s);
                        int n_min_span = std::min(n.query_e - n.query_s, n.ref_e - n.ref_s);
                        float n_score;
                        n_score = ( 2*n_min_span -  n_max_span) > 0 ? (float) (n.n_hits * ( 2*n_min_span -  n_max_span) ) : 1;   // this is really just n_hits * ( min_span - (offset_in_span) ) );
//                        n_score = n.n_hits * (n.query_e - n.query_s);
                        n.score = n_score;
                        final_nams.push_back(n);
                        max_nam_n_hits = std::max(n.n_hits, max_nam_n_hits);
                    }
                }

                // Remove all NAMs from open_matches that the current hit have passed
                unsigned int c = h.query_s;
                auto predicate = [c](decltype(open_nams)::value_type const &nam) { return nam.query_e < c; };
                open_nams.erase(std::remove_if(open_nams.begin(), open_nams.end(), predicate), open_nams.end());
                prev_q_start = h.query_s;
            }


        }

        // Add all current open_matches to final NAMs
        for (auto &n : open_nams){
            int n_max_span = std::max(n.query_e - n.query_s, n.ref_e - n.ref_s);
            int n_min_span = std::min(n.query_e - n.query_s, n.ref_e - n.ref_s);
            float n_score;
            n_score = ( 2*n_min_span -  n_max_span) > 0 ? (float) (n.n_hits * ( 2*n_min_span -  n_max_span) ) : 1;   // this is really just n_hits * ( min_span - (offset_in_span) ) );
//            n_score = n.n_hits * (n.query_e - n.query_s);
            n.score = n_score;
            final_nams.push_back(n);
            max_nam_n_hits = std::max(n.n_hits, max_nam_n_hits);
        }
    }
    info.second = max_nam_n_hits;
//    for (auto &n : final_nams){
//        int diff = (n.query_e - n.query_s) - (n.ref_e - n.ref_s);
//        std::cerr << "NAM ORG: nam_id: " << n.nam_id << " ref_id: " << n.ref_id << ": (" << n.score << ", " << n.n_hits << ", " << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e  << ")" << " diff: " << diff << " is_rc: " << n.is_rc << std::endl;
//    }
    return info;

//
//    std::cerr << "DONE" << std::endl;

//    return final_nams;
}




//static inline bool compareByQueryCoord(const hit &a, const hit &b)
//{
//    // first sort on ref ID, then on query, then on reference
//    return (a.ref_id < b.ref_id) ||
//           ( (a.ref_id == b.ref_id) && (a.query_s < b.query_s) ) ||
//           ((a.ref_id == b.ref_id) && (a.query_s == b.query_s ) && (a.ref_s < b.ref_s)) ;
//}
//
//static inline bool compareByQueryLength(const hit &a, const hit &b)
//{
//    return (a.query_e - a.query_s) < ( b.query_e - b.query_s);
//}

//static inline bool compareByNrHitsAndSimilarSpan(const nam &a, const nam &b)
//{
//    // first sort on nr hits, then on diff in span between query and reference, then on reference
//    return (a.n_hits > b.n_hits) ||
//           ( (a.n_hits == b.n_hits) && ( ((a.query_e - a.query_s) - (a.ref_e - a.ref_s)) < ((b.query_e - b.query_s) - (b.ref_e - b.ref_s)) ) );
//}

//static inline bool score(const nam &a, const nam &b)
//{
//    return ( (a.n_hits * (a.query_e - a.query_s)) > (b.n_hits * (b.query_e - b.query_s)) );
//}




inline void output_hits_paf(std::string &paf_output, const std::vector<nam> &all_nams, const std::string& query_acc, const ref_names &acc_map, int k, int read_len, const std::vector<unsigned int> &ref_len_map) {
    // Output results
    if (all_nams.size() == 0) {
        return;
    }
    // Only output single best hit based on: number of randstrobe-matches times span of the merged match.
    std::string o;
    nam n = all_nams[0];
    if (n.is_rc){
        o = "-";
    }
    else{
        o = "+";
    }
//    paf_output << query_acc << "\t" << read_len <<  "\t" << n.query_s << "\t" << n.query_prev_hit_startpos + k << "\t" << o  <<  "\t" << acc_map[n.ref_id] << "\t" << ref_len_map[n.ref_id] << "\t" << n.ref_s << "\t" << n.ref_prev_hit_startpos + k << "\t" << n.n_hits << "\t" << n.ref_prev_hit_startpos + k - n.ref_s << "\t" << "-" << "\n";
    paf_output.append(query_acc);
    paf_output.append("\t");
    paf_output.append(std::to_string(read_len));
    paf_output.append("\t");
    paf_output.append(std::to_string(n.query_s));
    paf_output.append("\t");
    paf_output.append(std::to_string(n.query_prev_hit_startpos + k));
    paf_output.append("\t");
    paf_output.append(o);
    paf_output.append("\t");
    paf_output.append(acc_map[n.ref_id]);
    paf_output.append("\t");
    paf_output.append(std::to_string( ref_len_map[n.ref_id]));
    paf_output.append("\t");
    paf_output.append(std::to_string(n.ref_s));
    paf_output.append("\t");
    paf_output.append(std::to_string(n.ref_prev_hit_startpos + k));
    paf_output.append("\t");
    paf_output.append(std::to_string( n.n_hits));
    paf_output.append("\t");
    paf_output.append(std::to_string(n.ref_prev_hit_startpos + k - n.ref_s));
    paf_output.append("\t255\n");
}

inline void output_hits_paf_PE(std::string &paf_output, nam &n, std::string &query_acc, ref_names &acc_map, int k, int read_len, std::vector<unsigned int> &ref_len_map) {
    // Output results
    std::string o;
    if (n.ref_s >= 0) {
        o = n.is_rc ? "-" : "+";
        paf_output.append(query_acc);
        paf_output.append("\t");
        paf_output.append(std::to_string(read_len));
        paf_output.append("\t");
        paf_output.append(std::to_string(n.query_s));
        paf_output.append("\t");
        paf_output.append(std::to_string(n.query_prev_hit_startpos + k));
        paf_output.append("\t");
        paf_output.append(o);
        paf_output.append("\t");
        paf_output.append(acc_map[n.ref_id]);
        paf_output.append("\t");
        paf_output.append(std::to_string( ref_len_map[n.ref_id]));
        paf_output.append("\t");
        paf_output.append(std::to_string(n.ref_s));
        paf_output.append("\t");
        paf_output.append(std::to_string(n.ref_prev_hit_startpos + k));
        paf_output.append("\t");
        paf_output.append(std::to_string(n.n_hits));
        paf_output.append("\t");
        paf_output.append(std::to_string(n.ref_prev_hit_startpos + k - n.ref_s));
        paf_output.append("\t255\n");
    }
}

static inline std::string reverse_complement(const std::string &read) {
    auto read_rev = read;
    std::reverse(read_rev.begin(), read_rev.end()); // reverse
    for (size_t j = 0; j < read_rev.length(); ++j) { // complement
        if (read_rev[j] == 'A') read_rev[j] = 'T';
        else if (read_rev[j] == 'T') read_rev[j] = 'A';
        else if (read_rev[j] == 'C') read_rev[j] = 'G';
        else if (read_rev[j] == 'G') read_rev[j] = 'C';
    }
    return read_rev;
}


//aln_info parasail_align(std::string &ref, int tlen, std::string &query, int qlen, int sc_mch, int sc_mis, int gapo, int gape) {
//    const char *ref_ptr = ref.c_str();
//    const char *read_ptr = query.c_str();
//    parasail_matrix_t *user_matrix = NULL;
//    user_matrix = parasail_matrix_create("ACGT", sc_mch, -sc_mis);
//    aln_info aln;
////    const std::string ref   = "AGTATCTGGAACTGGACTTTTGGAGCGCTTTCAGGGATAAGGTGAAAAAGGAAATATCTTCCCATAAAAACTGGACAGAAGCATTCTCAGAAACTTATTTGAGATGTGTGTACTCAACTAAGAGAATTGAACCACCGTTTTGAAGGAGCAGTTTTGAAACTCTCTTTTTCTGGAATCTGCAAGTGGATATTTGGCTAGCTTTGGGGATTTCGCTGGAAGCGGGAATACATATAAAAAGCACACAGCAGCGTTCTGAGAAACTGCTTTCTGATGTTTGCATTCAAGTCAAAAGTTGAACACTCCCTTTCATAGAGCAGTCTTGAAACACCCCTTTTGTAGTATCTGGAACTGGACTTTTGGAGCGATTTCAGGGCTAAGGTGAAAAAGGAAATATCTTCCCATAAAAACTGGACAGAAGCATTCTCAGAAACTTGGTTATGCTGTATCTACTCAACTAACAAAGTTGAACCTTTCTTTTGATAGAGCAGTTTTGAAATGGTCTTTTTGTGGAATCTGCAAGTGGATATTTGGCTAGTTTTGAGGATTTCGTTGGAAGCGGGAATTCATACAAATTGCAGACTGCAGCGTTCTGAGAAACATCTTTGTGATGTTTGTATTCAGGACAGAGAGTTGAACATTCCCTATCATAGAGCAGGTTGGAATCACTCCTTTTGTAGTATCTGGAAGTGGACATTTGGAGCGCTTTCAGGCCTATTTTGGAAAGGGAAATATCTTCCCGTAACAACTATGCAGAAGCATTCTCAGAAACTTGTTTGTGATGTGTGCCCTCTACTGACAGAGTTGAACCTTTCTTTTCATAGAGCAGTTTTGAAACACTCTTTTTGTAGAA";
////    const std::string query = "CGGGAATACATATAAAAAGCACACAGCAGCGTTCTGAGAAACTGCTTTCTGATGTTTGCATTAAAGTCAAAAGTTGAACACTCCCTTTCATAGAGCAGTC";
//
//
////    parasail_result_t *result = NULL;
////    result = parasail_sw_trace_striped_sat(read_ptr, qlen, ref_ptr, tlen,  gapo, gape, user_matrix);
////    parasail_result_free(result);
//
//        parasail_result_ssw *result = NULL;
//        result = parasail_ssw(read_ptr, qlen, ref_ptr, tlen, gapo, gape, user_matrix );
//
//    // TODO: Fix cigarstring to use M instead of =/X
////    aln.ed = result->;
//    aln.ref_offset = result->ref_begin1;
//    aln.sw_score = result->score1; //(alignment.query_end - alignment.query_begin) - 4*alignment.mismatches; //approximate for ssw until I implement a cigar parser
//
//    std::stringstream cigar_string;
//    int edit_distance = 0;
//    int sw_score = 0;
//    unsigned ref_pos = 0, read_pos = 0;
//
//    for (int i = 0; i < result->cigarLen; i++) {
//        int count = result->cigar[i] >> 4;
//        char op = "MID"[result->cigar[i] & 0xf];
////        std::cerr << "count: " << count << " op:" << op << std::endl;
//        if ( (i==0) && op == 'D'){
//            ref_pos += count;
////            std::cerr << "First deletion " << i << " " << count << std::endl;
//            continue;
//        }
//        if ( (i==result->cigarLen-1) && op == 'D'){
//            ref_pos += count;
////            std::cerr << "Last deletion " << i << " " << count << std::endl;
//            continue;
//        }
//        cigar_string << count << op;
//        switch (op) {
//            case 'M':
//                for (int j = 0; j < count; j++, ref_pos++, read_pos++) {
//                    if (ref_ptr[ref_pos] != read_ptr[read_pos]) {
//                        edit_distance++;
//                        sw_score -= sc_mis;
//                    } else{
//                        sw_score += sc_mch;
//                    }
//                }
//                break;
//            case 'D':edit_distance += count;
//                ref_pos += count;
//                sw_score -= (gapo + (count-1));
//                break;
//            case 'I':edit_distance += count;
//                read_pos += count;
//                sw_score -= (gapo + (count-1));
//                break;
//            default:assert(0);
//        }
////        std::cerr << "ED " << edit_distance << std::endl;
//    }
//
//    aln.cigar =  cigar_string.str();
//    aln.ed =  edit_distance;
//
//    parasail_result_ssw_free(result);
//
//    return aln;
//}


//inline aln_info ssw_align(std::string &ref, std::string &query, int read_len, int match_score, int mismatch_penalty, int gap_opening_penalty, int gap_extending_penalty) {
//
//    aln_info aln;
//    int32_t maskLen = strlen(query.c_str())/2;
//    maskLen = std::max(maskLen, 15);
//    if (ref.length() > 2000){
////        std::cerr << "ALIGNMENT TO REF LONGER THAN 2000bp - REPORT TO DEVELOPER. Happened for read: " <<  query << " ref len:" << ref.length() << std::endl;
//        aln.global_ed = 100000;
//        aln.ed = 100000;
//        aln.ref_offset = 0;
//        aln.cigar = "*";
//        aln.sw_score = -1000000;
//        return aln;
//    }
//
//    StripedSmithWaterman::Aligner aligner(match_score, mismatch_penalty, gap_opening_penalty, gap_extending_penalty);
////    StripedSmithWaterman::Aligner aligner;
//    // Declares a default filter
//    StripedSmithWaterman::Filter filter;
//    // Declares an alignment that stores the result
//    StripedSmithWaterman::Alignment alignment_ssw;
//    aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment_ssw, maskLen, 1);
//        // Have to give up this optimization untill the 'Command terminated abnormally' bug is fixed in ssw library
////     if (read_len*match_score < 255){
////         std::cerr << "Here: "  << read_len*match_score << " " << ref.length() << std::endl;
////         try
////         {
////             aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment_ssw, maskLen, 0);
////         }
////         catch (...)
////         {
////             aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment_ssw, maskLen, 1);
////         }
////
////     } else {
////            aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment_ssw, maskLen, 1);
////     }
////    std::cerr << passed << std::endl;
////    if(!passed){
////        std::cerr << "Failed" << std::endl;
////        std::cerr << "read: " << query << std::endl;
////        std::cerr << "ref: "  << ref << std::endl;
////    }
//
//
////    std::cerr << "===== SSW result =====" << std::endl;
////    std::cerr << "Best Smith-Waterman score:\t" << alignment_ssw.sw_score << std::endl
////         << "Next-best Smith-Waterman score:\t" << alignment_ssw.sw_score_next_best << std::endl
////         << "Reference start:\t" << alignment_ssw.ref_begin << std::endl
////         << "Reference end:\t" << alignment_ssw.ref_end << std::endl
////         << "Query start:\t" << alignment_ssw.query_begin << std::endl
////         << "Query end:\t" << alignment_ssw.query_end << std::endl
////         << "Next-best reference end:\t" << alignment_ssw.ref_end_next_best << std::endl
////         << "Number of mismatches:\t" << alignment_ssw.mismatches << std::endl
////         << "Cigar: " << alignment_ssw.cigar_string << std::endl;
//
//    aln.global_ed = alignment_ssw.global_ed;
//    aln.ed = alignment_ssw.mismatches;
//    aln.ref_offset = alignment_ssw.ref_begin;
//    aln.cigar = alignment_ssw.cigar_string;
//    aln.sw_score = alignment_ssw.sw_score;
//    aln.length = alignment_ssw.ref_end - alignment_ssw.ref_begin;
//    return aln;
//}


//inline aln_info ksw_align(const char *tseq, int tlen, const char *qseq, int qlen,
//                          int sc_mch, int sc_mis, int gapo, int gape, ksw_extz_t &ez) {
//    int8_t a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
//    int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a, b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
//    const uint8_t *ts = reinterpret_cast<const uint8_t *>(tseq);
//    const uint8_t *qs = reinterpret_cast<const uint8_t *>(qseq);
//    memset(&ez, 0, sizeof(ksw_extz_t));
//    ksw_extz2_sse(0, qlen, qs, tlen, ts, 5, mat, gapo, gape, -1, -1, 10000, KSW_EZ_EXTZ_ONLY, &ez);
//
//    aln_info aln;
////    std::string cigar_mod;
////    cigar_mod.reserve(5*ez.n_cigar);
//    unsigned int tstart_offset = 0;
//    int eqx_len, switch_ind;
//    std::stringstream cigar_string;
//    int edit_distance = 0;
//    int sw_score = 0;
//    unsigned ref_pos = 0, read_pos = 0;
//    for (int i = 0; i < ez.n_cigar; i++) {
//        int count = ez.cigar[i] >> 4;
//        char op = "MID"[ez.cigar[i] & 0xf];
////        std::cerr << "count: " << count << " op:" << op << std::endl;
//        if ( (i==0) && op == 'D'){
//            ref_pos += count;
//            tstart_offset = ref_pos;
////            std::cerr << "First deletion " << i << " " << count << std::endl;
//            continue;
//        }
//        if ( (i==ez.n_cigar-1) && op == 'D'){
//            ref_pos += count;
////            std::cerr << "Last deletion " << i << " " << count << std::endl;
//            continue;
//        }
//        cigar_string << count << op;
//        switch (op) {
//            case 'M': {
////                eqx_len = 0;
////                switch_ind = 0; // switch_ind 0 if prev was match, 1 if mismatch
////                char o = '=';
//                for (int j = 0; j < count; j++, ref_pos++, read_pos++) {
//                    if (tseq[ref_pos] != qseq[read_pos]) {
//                        edit_distance++;
//                        sw_score -= -b;
////                        if ((switch_ind == 0) && (j > 0)) { // prev was match
////                            cigar_string << eqx_len << '=';
////                            eqx_len = 0;
////                        }
////                        switch_ind = 1;
////                        o = 'X';
////                        eqx_len++;
//                    } else{
//                        sw_score += sc_mch;
////                        if (switch_ind == 1) { // prev was mismatch
////                            cigar_string << eqx_len << 'X';
////                            eqx_len = 0;
////                            o = '=';
////                            switch_ind = 0;
////                        }
////                        eqx_len++;
//                    }
//                }
////                cigar_string << eqx_len << o;
//                break;
//            }
//            case 'D': {
//                edit_distance += count;
//                ref_pos += count;
//                sw_score -= (gapo + (count - 1));
////                cigar_string << count << op;
//                break;
//            }
//            case 'I': {
//                edit_distance += count;
//                read_pos += count;
//                sw_score -= (gapo + (count - 1));
////                cigar_string << count << op;
//                break;
//            }
//            default:assert(0);
//        }
////        std::cerr << "ED " << edit_distance << std::endl;
//    }
//    aln.ed = edit_distance;
//    aln.sw_score = sw_score;
//    aln.ref_offset = tstart_offset;
//    aln.cigar = cigar_string.str();
//    free(ez.cigar); //free(ts); free(qs);
//    return aln;
//}

inline int HammingDistance(const std::string &One, const std::string &Two)
{
    if (One.length() != Two.length()){
        return -1;
    }

    int counter = 0;

    for(int i=0; i<One.length(); i++) {
        if (One[i] != Two[i]) counter++;
    }

    return counter;
}


inline int HammingToCigarEQX2(const std::string &One, const std::string &Two, std::stringstream &cigar, int match, int mismatch, int &aln_score, int &soft_left, int &soft_right)
{
    if (One.length() != Two.length()){
        return -1;
    }

    // Decide softclipps
    int peak_score = 0;
//    int peak_score_pos = 0;
    int curr_score = 0;
    int end_softclipp = 0;
    for(int i=0; i<One.length(); i++) {
        if (One[i] == Two[i]){
            curr_score += match;
        } else {
            curr_score -= mismatch;
        }

        if (curr_score >= peak_score){
            peak_score = curr_score;
//            peak_score_pos = i;
            end_softclipp = i;
        }
    }
//    std::cout << "End softclipp: " << end_softclipp << std::endl;

    peak_score = 0;
//    peak_score_pos = 0;
    curr_score = 0;
    int start_softclipp = 0;
    for (int i = One.length() - 1; i >= 0; i--) {
        if (One[i] == Two[i]){
            curr_score += match;
        } else {
            curr_score -= mismatch;
        }

        if (curr_score >= peak_score){
            peak_score = curr_score;
//            peak_score_pos = i;
            start_softclipp = i;
        }
    }

    if (start_softclipp >= end_softclipp){
        aln_score = 0;
        soft_left = 50; //default
        soft_right = 50; //default
        return -1;
    }
//    std::cout << "Start softclipp: " << start_softclipp << std::endl;

    if (start_softclipp > 0){
        cigar << start_softclipp << 'S';
    }

    int counter = 1;
    bool prev_is_match = One[start_softclipp+1] == Two[start_softclipp+1];
    int hamming_mod = prev_is_match ? 0 : 1;
    bool curr_match = false;
    for (int i = start_softclipp + 1; i < end_softclipp + 1; i++) {
        curr_match = (One[i] == Two[i]);

        if ( !curr_match && prev_is_match ){
            cigar << counter << '=';
            aln_score += counter * match;
            counter = 0;
        }
        else if ( curr_match && !prev_is_match ){
            cigar << counter << 'X';
            aln_score -= counter * mismatch;
            hamming_mod += counter;
//            std::cout << "Added1: " << counter << " current: " << hamming_mod << std::endl;
//            if (beginning){
//                needs_aln = counter > 2;
//            }
            counter = 0;
        }
        prev_is_match = curr_match;
        counter++;
    }

    // Print last
    if (curr_match) {
        cigar << counter << '=';
        aln_score += counter * match;
    } else {
        cigar << counter << 'X';
        hamming_mod += counter;
//        std::cout << "Added2: " << counter << " current: " << hamming_mod << std::endl;
//        needs_aln = counter > 2;
    }

    if (One.length() - end_softclipp - 1 > 0){
        cigar << One.length() - end_softclipp - 1 << 'S';
    }

//    if (aln_score < 0){
//            std::cout << "NEGATIVE SCORE: " << aln_score << std::endl;
//    }

    soft_left = start_softclipp;
    soft_right = One.length() - end_softclipp - 1 > 0 ? One.length() - end_softclipp - 1 : 0;

    return hamming_mod;
}


inline bool HammingToCigarEQX(const std::string &One, const std::string &Two, std::stringstream &cigar)
{
    if (One.length() != Two.length()){
        return true;
    }

    int counter = 1;
    bool prev_is_match = One[0] == Two[0];
    bool beginning = true;
    bool needs_aln = false; // !prev_is_match
    bool curr_match;
    for(int i=1; i<One.length(); i++) {
        curr_match = (One[i] == Two[i]);

        if ( !curr_match && prev_is_match ){
            cigar << counter << '=';
            counter = 0;
        }
        else if ( curr_match && !prev_is_match ){
            cigar << counter << 'X';
            if (beginning){
                needs_aln = counter > 2;
            }
            counter = 0;
        }
        prev_is_match = curr_match;
        counter++;
    }

    // Print last
    if ( curr_match  ){
        cigar << counter << '=';
    } else{
        cigar << counter << 'X';
        needs_aln = counter > 2;
    }
    return needs_aln;
}


static inline bool sort_lowest_ed_scores_single(const std::tuple<int, alignment> &a,
                                                const std::tuple<int, alignment> &b)
{
    return (std::get<0>(a) < std::get<0>(b));
}

static inline bool sort_highest_sw_scores_single(const std::tuple<int, alignment> &a,
                                                 const std::tuple<int, alignment> &b)
{
    return (std::get<0>(a) > std::get<0>(b));
}


inline void align_SE(const alignment_params &aln_params, Sam& sam, std::vector<nam> &all_nams, const KSeq& record, int k, const std::vector<unsigned int> &ref_len_map, const std::vector<std::string> &ref_seqs, logging_variables &log_vars, float dropoff, int max_tries) {
    auto query_acc = record.name;
    auto read = record.seq;
    auto qual = record.qual;
    auto read_len = read.size();

    std::string read_rc;
    bool rc_already_comp = false;

    if (all_nams.size() == 0) {
        sam.add_unmapped(record);
        return;
    }

    int cnt = 0;
    float score_dropoff;
    nam n_max = all_nams[0];
    float s1 = n_max.score;
    // old mapq commented out - it is based only on seeds which is worse than calculating them based on base level alignments (even though heuristic)
//    if (all_nams.size() > 1) {
//        nam n_second = all_nams[1];
//        float s2 = n_second.score;
//        float min_matches;
//        min_matches  = (float)n_max.n_hits/10 > 1 ? (float)n_max.n_hits/10 : 1;
//        mapq = 40*(1 - s2/s1)*min_matches*log(s1) < 60 ? 40*(1 - s2/s1)*min_matches*log(s1) : 60 ;
//    }
//    int extra_ref = 50;
    int best_align_dist = ~0U >> 1;
    int best_align_sw_score = -1000;

    bool aln_did_not_fit;
    alignment sam_aln;
    int mapq = 60;
    int min_mapq_diff = best_align_dist;

    for (auto &n : all_nams) {
        aln_did_not_fit = false;
        score_dropoff = (float) n.n_hits / n_max.n_hits;
//        score_dropoff = (float) n.score / n_max.score;

        if ( (cnt >= max_tries) || best_align_dist == 0 || score_dropoff < dropoff){ // only consider top 20 hits as minimap2 and break if alignment is exact match to reference or the match below droppoff cutoff.
            break;
        }

        log_vars.tot_all_tried ++;

        int ref_diff = n.ref_e - n.ref_s;
        int read_diff = n.query_e - n.query_s;
        int min_diff =  read_diff ^ ((ref_diff ^ read_diff) & -(ref_diff < read_diff));
        int max_diff = ref_diff ^ ((ref_diff ^ read_diff) & -(ref_diff < read_diff));
        int diff = max_diff - min_diff;


//        // decide if read should be fw or rc aligned to reference here by checking exact match of first and last strobe in the NAM
//
//        if ( (ref_segm.substr(n.query_s, k) == read.substr(n.query_s, k) ) ) { //&& (ref_segm.substr(n.query_e - k + (ref_diff - read_diff), k) == read.substr(n.query_e - k, k)) ){
//            n.is_rc = false;
//        }
//        else {
//            if (!rc_already_comp){
//                read_rc = reverse_complement(read);
//                rc_already_comp = true;
//            }
//
//            if ((ref_segm.substr(n.query_s, k) == read_rc.substr(n.query_s, k))) { // && (ref_segm.substr(n.query_e - k + (ref_diff - read_diff), k) == read_rc.substr(n.query_e - k, k)) ){
//                n.is_rc = true;
//            } else {
//                did_not_fit++;
//                aln_did_not_fit = true;
//            }
//        }

        // decide if read should be fw or rc aligned to reference here by checking exact match of first and last strobe in the NAM
        bool fits = false;
        std::string ref_start_kmer;
        std::string ref_end_kmer;
        std::string read_start_kmer;
        std::string read_end_kmer;
        std::string read_rc_start_kmer;
        std::string read_rc_end_kmer;
        ref_start_kmer = ref_seqs[n.ref_id].substr(n.ref_s, k);
        ref_end_kmer = ref_seqs[n.ref_id].substr(n.ref_e-k, k);

        if (!n.is_rc) {
            read_start_kmer = read.substr(n.query_s, k);
            read_end_kmer = read.substr(n.query_e-k, k);
            if ((ref_start_kmer == read_start_kmer) && (ref_end_kmer == read_end_kmer)) {
//            n.is_rc = false;
                fits = true;
            } else  {
                //  FALSE FORWARD TAKE CARE OF FALSE HITS HERE - it can be false forwards or false rc because of symmetrical hash values
                //    we need two extra checks for this - hopefully this will remove all the false hits we see (true hash collisions should be very few)

//              std::cerr << " CHECKING1!! " << std::endl;
                // false reverse hit, change coordinates in nam to forward
                if (!rc_already_comp){
                    read_rc = reverse_complement(read);
                    rc_already_comp = true;
                }

                int q_start_tmp = read_len - n.query_e;
                int q_end_tmp = read_len - n.query_s;
                read_start_kmer = read_rc.substr(q_start_tmp, k);
                read_end_kmer = read_rc.substr(q_end_tmp-k, k);
                if ((ref_start_kmer == read_start_kmer) && (ref_end_kmer == read_end_kmer)){
                    fits = true;
                    n.is_rc = true;
                    n.query_s = q_start_tmp;
                    n.query_e = q_end_tmp;
//                std::cerr << " DETECTED FALSE RC FROM SYMM!! " << std::endl;
                }

            }
        } else {
            if (!rc_already_comp){
                read_rc = reverse_complement(read);
                rc_already_comp = true;
            }
            read_rc_start_kmer = read_rc.substr(n.query_s, k);
            read_rc_end_kmer = read_rc.substr(n.query_e-k, k);
            if ( (ref_start_kmer == read_rc_start_kmer) && (ref_end_kmer == read_rc_end_kmer) ) { // && (ref_segm.substr(n.query_e - k + (ref_diff - read_diff), k) == read_rc.substr(n.query_e - k, k)) ){
                n.is_rc = true;
                fits = true;
            } else{
                //  FALSE REVERSE TAKE CARE OF FALSE HITS HERE - it can be false forwards or false rc because of symmetrical hash values
                //    we need two extra checks for this - hopefully this will remove all the false hits we see (true hash collisions should be very few)

                int q_start_tmp = read_len - n.query_e;
                int q_end_tmp = read_len - n.query_s;
                read_start_kmer = read.substr(q_start_tmp, k);
                read_end_kmer = read.substr(q_end_tmp-k, k);
//            std::cerr << " CHECKING2!! " <<   n.query_s << " " <<   n.query_e << " " << std::endl;
//            std::cerr << read_start_kmer  << " " <<  ref_start_kmer << " " <<  read_end_kmer << " " << ref_end_kmer << std::endl;

                if ((ref_start_kmer == read_start_kmer) && (ref_end_kmer == read_end_kmer)){
                    fits = true;
                    n.is_rc = false;
                    n.query_s = q_start_tmp;
                    n.query_e = q_end_tmp;
//                std::cerr << " DETECTED FALSE FW FROM SYMM!! " << std::endl;
                }
            }
        }
        if (!fits){
            log_vars.did_not_fit++;
            aln_did_not_fit = true;
        }

        // deal with any read hanging of ends of reference not to get 'std::out_of_range' what(): basic_string::substr
        int ref_tmp_start = n.ref_s - n.query_s;
        int ref_tmp_segm_size = read_len + diff;
        int ref_len = ref_len_map[n.ref_id];
        int ref_start = std::max(0, ref_tmp_start);
        int ref_segm_size = ref_tmp_segm_size < ref_len - ref_start ? ref_tmp_segm_size : ref_len - 1 - ref_start;
        std::string ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_segm_size);

        int soft_left = 50;
        int soft_right = 50;
        int hamming_mod;
        int hamming_dist = -1;
        std::string r_tmp;
        bool is_rc;
        if (n.is_rc){
            r_tmp = read_rc;
            is_rc = true;
        }else{
            r_tmp = read;
            is_rc = false;
        }
//        std::cout << "DIFF: "  <<  diff << ", " << n.score << ", " << ref_segm.length() << std::endl;

        if (ref_segm.length() == read_len){
            hamming_dist = HammingDistance(r_tmp, ref_segm);
//            std::cout << "Hammingdist: " << n.score << ", "  <<  n.n_hits << ", " << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e  << ") hd:" << hamming_dist << ", best ed so far: " << best_align_dist  << std::endl;
            if ( (hamming_dist >=0)){
                int sw_score = aln_params.match*(read_len-hamming_dist) - aln_params.mismatch*hamming_dist;
                int diff_to_best = std::abs(best_align_sw_score - sw_score);
                min_mapq_diff = std::min(min_mapq_diff, diff_to_best);
                std::stringstream cigar_string;
                int aln_score = 0;
                hamming_mod = HammingToCigarEQX2(r_tmp, ref_segm, cigar_string, aln_params.match, aln_params.mismatch, aln_score, soft_left, soft_right);
//                if (hamming_dist < best_align_dist){
                if (aln_score > best_align_sw_score){
                    min_mapq_diff = std::max(0, sw_score - best_align_sw_score); // new distance to next best match

//                    min_mapq_diff = best_align_dist - hamming_dist; // new distance to next best match
//                    needs_aln = HammingToCigarEQX(r_tmp, ref_segm, cigar_string);

//                    sw_score = aln_params.match*(read_len-hamming_mod) - aln_params.mismatch*hamming_mod;

                    sam_aln.cigar = cigar_string.str();
                    best_align_dist = hamming_dist;
                    //                sam_aln.cigar = std::to_string(read_len) + "M";
                    sam_aln.ed = hamming_mod;
                    sam_aln.ref_start = ref_start + soft_left +1; // +1 because SAM is 1-based!
                    sam_aln.is_rc = is_rc;
                    sam_aln.ref_id = n.ref_id;
                    sam_aln.sw_score = aln_score;
                    best_align_sw_score = aln_score;
                    sam_aln.aln_score = aln_score;
                    sam_aln.aln_length = read_len;
                }
            }
        }

//        std::cout << hamming_dist  << ", " << read_len << ", " << (float) hamming_dist / (float) read_len << std::endl;

        // ((float) sam_aln.ed / read_len) < 0.05  //Hamming distance worked fine, no need to ksw align
        if ( (hamming_dist >=0) && (diff == 0) && (((float) hamming_dist / (float) read_len) < 0.05) ) { // Likely substitutions only (within NAM region) no need to call ksw alingment
            ;
//            if (hamming_dist < best_align_dist){
//                ;
////                best_align_dist = hamming_dist;
//////                sam_aln.cigar = std::to_string(read_len) + "M";
////                sam_aln.ed = hamming_mod;
////                sam_aln.ref_start = ref_start +1; // +1 because SAM is 1-based!
////                sam_aln.is_rc = is_rc;
////                sam_aln.ref_id = n.ref_id;
//            }
//            std::cout << "HERE 1 " << sam_aln.ref_start << " " << hamming_dist  << ", " << read_len << ", " << (float) hamming_dist / (float) read_len << std::endl;

        } else {

//        } else if ( (best_align_dist > 1) || aln_did_not_fit ){
//        } else if ( (best_align_dist > 1) || ( aln_did_not_fit || needs_aln || (((float) hamming_dist / (float) read_len) >= 0.05) ) ){
            int extra_ref_left = std::min(soft_left, 50);
            int extra_ref_right = std::min(soft_right, 50);
//            extra_ref = 50; //(read_diff - ref_diff) > 0 ?  (read_diff - ref_diff) : 0;
            int a = n.ref_s - n.query_s - extra_ref_left;
            int ref_start = std::max(0, a);
            int b = n.ref_e + (read_len - n.query_e)+ extra_ref_right;
            int ref_len = ref_len_map[n.ref_id];
            int ref_end = std::min(ref_len, b);
            ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_end - ref_start);
//            ksw_extz_t ez;
//            const char *ref_ptr = ref_segm.c_str();
//            const char *read_ptr = r_tmp.c_str();
            aln_info info;
//            std::cout << "Extra ref: " << extra_ref << " " << read_diff << " " << ref_diff << " " << ref_start << " " << ref_end << std::endl;
//            info = ksw_align(ref_ptr, ref_segm.size(), read_ptr, r_tmp.size(), 1, 4, 6, 1, ez);
            info = ssw_align(ref_segm, r_tmp, aln_params.match, aln_params.mismatch, aln_params.gap_open, aln_params.gap_extend);
//            info.ed = info.global_ed; // read_len - info.sw_score;
            int diff_to_best = std::abs(best_align_sw_score - info.sw_score);
            min_mapq_diff = std::min(min_mapq_diff, diff_to_best);
            log_vars.tot_ksw_aligned ++;
//            if (info.global_ed <= best_align_dist){
            if (info.sw_score >= best_align_sw_score){
                min_mapq_diff = std::max(0, info.sw_score - best_align_sw_score); // new distance to next best match
                best_align_dist = info.global_ed;
                sam_aln.cigar = info.cigar;
                sam_aln.ed = info.ed;
                sam_aln.ref_start =  ref_start + info.ref_offset +1; // +1 because SAM is 1-based!
                sam_aln.is_rc = is_rc;
                sam_aln.ref_id = n.ref_id;
                sam_aln.sw_score = info.sw_score;
                best_align_sw_score = info.sw_score;
                sam_aln.aln_score = info.sw_score;
                sam_aln.aln_length = info.length;
            }
//            std::cout << "HERE 2 "  << sam_aln.ref_start << " global ed: " << info.global_ed  <<  ", hamming: " << hamming_dist << ", "<< read_len << ", " << (float) hamming_dist / (float) read_len << std::endl;

//            std::cout << "Aligned: " << n.score << ", "  << n.n_hits << ", " << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e  << ") ed:" << info.ed << ", best ed so far: " << best_align_dist  << std::endl;

        }
        cnt ++;
    }

    if (all_nams.size() > 0) {
        sam_aln.mapq = std::min(min_mapq_diff, 60);
        sam.add(sam_aln, record, read_rc);
    }
}


static inline void align_SE_secondary_hits(alignment_params &aln_params, Sam& sam, std::vector<nam> &all_nams, const KSeq& record, int k, std::vector<unsigned int> &ref_len_map, std::vector<std::string> &ref_seqs, logging_variables &log_vars, float dropoff, int max_tries, int max_secondary ) {

    auto query_acc = record.name;
    auto read = record.seq;
    auto qual = record.qual;
    auto read_len = read.size();
    std::string read_rc;
    bool rc_already_comp = false;

    if (all_nams.size() == 0) {
        sam.add_unmapped(record);
        return;
    }

    std::vector<std::tuple<int,alignment>> alignments; // (score, aln)
    int cnt = 0;
    float score_dropoff;
    nam n_max = all_nams[0];
    float s1 = n_max.score;

//    int extra_ref = 50;
    int best_align_dist = ~0U >> 1;
    int best_align_sw_score = -1000;

    int min_mapq_diff = best_align_dist;
    bool aln_did_not_fit;
//    int best_align_sw_score = -1000;
    for (auto &n : all_nams) {
        alignment sam_aln;
        sam_aln.ed = 1000; // init
        aln_did_not_fit = false;
        score_dropoff = (float) n.n_hits / n_max.n_hits;
//        score_dropoff = (float) n.score / n_max.score;

        if ( (cnt >= max_tries) || best_align_dist == 0 || score_dropoff < dropoff){ // only consider top 20 hits as minimap2 and break if alignment is exact match to reference or the match below droppoff cutoff.
            break;
        }

        log_vars.tot_all_tried ++;

        int ref_diff = n.ref_e - n.ref_s;
        int read_diff = n.query_e - n.query_s;
        int min_diff =  read_diff ^ ((ref_diff ^ read_diff) & -(ref_diff < read_diff));
        int max_diff = ref_diff ^ ((ref_diff ^ read_diff) & -(ref_diff < read_diff));
        int diff = max_diff - min_diff;

        // decide if read should be fw or rc aligned to reference here by checking exact match of first and last strobe in the NAM
        bool fits = false;
        std::string ref_start_kmer;
        std::string ref_end_kmer;
        std::string read_start_kmer;
        std::string read_end_kmer;
        std::string read_rc_start_kmer;
        std::string read_rc_end_kmer;
        ref_start_kmer = ref_seqs[n.ref_id].substr(n.ref_s, k);
        ref_end_kmer = ref_seqs[n.ref_id].substr(n.ref_e-k, k);

        if (!n.is_rc) {
            read_start_kmer = read.substr(n.query_s, k);
            read_end_kmer = read.substr(n.query_e-k, k);
            if ((ref_start_kmer == read_start_kmer) && (ref_end_kmer == read_end_kmer)) {
//            n.is_rc = false;
                fits = true;
            } else  {
                //  FALSE FORWARD TAKE CARE OF FALSE HITS HERE - it can be false forwards or false rc because of symmetrical hash values
                //    we need two extra checks for this - hopefully this will remove all the false hits we see (true hash collisions should be very few)

//              std::cerr << " CHECKING1!! " << std::endl;
                // false reverse hit, change coordinates in nam to forward
                if (!rc_already_comp){
                    read_rc = reverse_complement(read);
                    rc_already_comp = true;
                }

                int q_start_tmp = read_len - n.query_e;
                int q_end_tmp = read_len - n.query_s;
                read_start_kmer = read_rc.substr(q_start_tmp, k);
                read_end_kmer = read_rc.substr(q_end_tmp-k, k);
                if ((ref_start_kmer == read_start_kmer) && (ref_end_kmer == read_end_kmer)){
                    fits = true;
                    n.is_rc = true;
                    n.query_s = q_start_tmp;
                    n.query_e = q_end_tmp;
//                std::cerr << " DETECTED FALSE RC FROM SYMM!! " << std::endl;
                }

            }
        } else {
            if (!rc_already_comp){
                read_rc = reverse_complement(read);
                rc_already_comp = true;
            }
            read_rc_start_kmer = read_rc.substr(n.query_s, k);
            read_rc_end_kmer = read_rc.substr(n.query_e-k, k);
            if ( (ref_start_kmer == read_rc_start_kmer) && (ref_end_kmer == read_rc_end_kmer) ) { // && (ref_segm.substr(n.query_e - k + (ref_diff - read_diff), k) == read_rc.substr(n.query_e - k, k)) ){
                n.is_rc = true;
                fits = true;
            } else{
                //  FALSE REVERSE TAKE CARE OF FALSE HITS HERE - it can be false forwards or false rc because of symmetrical hash values
                //    we need two extra checks for this - hopefully this will remove all the false hits we see (true hash collisions should be very few)

                int q_start_tmp = read_len - n.query_e;
                int q_end_tmp = read_len - n.query_s;
                read_start_kmer = read.substr(q_start_tmp, k);
                read_end_kmer = read.substr(q_end_tmp-k, k);
//            std::cerr << " CHECKING2!! " <<   n.query_s << " " <<   n.query_e << " " << std::endl;
//            std::cerr << read_start_kmer  << " " <<  ref_start_kmer << " " <<  read_end_kmer << " " << ref_end_kmer << std::endl;

                if ((ref_start_kmer == read_start_kmer) && (ref_end_kmer == read_end_kmer)){
                    fits = true;
                    n.is_rc = false;
                    n.query_s = q_start_tmp;
                    n.query_e = q_end_tmp;
//                std::cerr << " DETECTED FALSE FW FROM SYMM!! " << std::endl;
                }
            }
        }
        if (!fits){
            log_vars.did_not_fit++;
            aln_did_not_fit = true;
        }

        // deal with any read hanging of ends of reference not to get 'std::out_of_range' what(): basic_string::substr
        int ref_tmp_start = n.ref_s - n.query_s;
        int ref_tmp_segm_size = read_len + diff;
        int ref_len = ref_len_map[n.ref_id];
        int ref_start = std::max(ref_tmp_start, 0);
        int ref_segm_size = ref_tmp_segm_size < ref_len - ref_start ? ref_tmp_segm_size : ref_len - 1 - ref_start;

        std::string ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_segm_size);

//        // decide if read should be fw or rc aligned to reference here by checking exact match of first and last strobe in the NAM
//
//        if ( (ref_segm.substr(n.query_s, k) == read.substr(n.query_s, k) ) ) { //&& (ref_segm.substr(n.query_e - k + (ref_diff - read_diff), k) == read.substr(n.query_e - k, k)) ){
//            n.is_rc = false;
//        }
//        else {
//            if (!rc_already_comp){
//                read_rc = reverse_complement(read);
//                rc_already_comp = true;
//            }
//
//            if ((ref_segm.substr(n.query_s, k) == read_rc.substr(n.query_s, k))) { // && (ref_segm.substr(n.query_e - k + (ref_diff - read_diff), k) == read_rc.substr(n.query_e - k, k)) ){
//                n.is_rc = true;
//            } else {
//                did_not_fit++;
//                aln_did_not_fit = true;
//            }
//        }

        int hamming_mod;
        int soft_left = 50;
        int soft_right = 50;
//        bool needs_aln = false;
        int hamming_dist = -1;
        int sw_score = -999;
        std::string r_tmp;
        bool is_rc;
        if (n.is_rc){
            r_tmp = read_rc;
            is_rc = true;
        }else{
            r_tmp = read;
            is_rc = false;
        }
//        std::cout << "DIFF: "  <<  diff << ", " << n.score << ", " << ref_segm.length() << std::endl;

        if (ref_segm.length() == read_len){
            hamming_dist = HammingDistance(r_tmp, ref_segm);
//            std::cout << "Hammingdist: " << n.score << ", "  <<  n.n_hits << ", " << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e  << ") hd:" << hamming_dist << ", best ed so far: " << best_align_dist  << std::endl;
            if ( (hamming_dist >=0)){
                sw_score =  aln_params.match*(read_len-hamming_dist) - aln_params.mismatch*hamming_dist;
                int diff_to_best = std::abs(best_align_sw_score - sw_score);
                min_mapq_diff = std::min(min_mapq_diff, diff_to_best);
//                if (hamming_dist < best_align_dist) {
//                    min_mapq_diff = best_align_dist - hamming_dist; // new distance to next best match
//                }
                if (sw_score > best_align_sw_score){
                    min_mapq_diff = std::max(0, sw_score - best_align_sw_score); // new distance to next best match
                }
                std::stringstream cigar_string;
//                needs_aln = HammingToCigarEQX(r_tmp, ref_segm, cigar_string);
                int aln_score = 0;
                hamming_mod = HammingToCigarEQX2(r_tmp, ref_segm, cigar_string, aln_params.match, aln_params.mismatch, aln_score, soft_left, soft_right);
//                sw_score =  aln_params.match*(read_len-hamming_mod) - aln_params.mismatch*hamming_mod;

                sam_aln.cigar = cigar_string.str();
                best_align_dist = hamming_dist;
                //                sam_aln.cigar = std::to_string(read_len) + "M";
                sam_aln.global_ed = hamming_dist;
                sam_aln.ed = hamming_mod;
                sam_aln.ref_start = ref_start + soft_left +1; // +1 because SAM is 1-based!
                sam_aln.is_rc = is_rc;
                sam_aln.ref_id = n.ref_id;
                sam_aln.sw_score = sw_score;
                sam_aln.aln_score = aln_score;
                sam_aln.aln_length = read_len;

//                best_align_sw_score = sam_aln.sw_score;

            }
        }
        if ( (hamming_dist >=0) && (diff == 0) && (((float) hamming_dist / (float) read_len) < 0.05) ) { // Likely substitutions only (within NAM region) no need to call ksw alingment
            ;
//            if (hamming_dist < best_align_dist){
//                ;
////                best_align_dist = hamming_dist;
//////                sam_aln.cigar = std::to_string(read_len) + "M";
////                sam_aln.ed = hamming_mod;
////                sam_aln.ref_start = ref_start +1; // +1 because SAM is 1-based!
////                sam_aln.is_rc = is_rc;
////                sam_aln.ref_id = n.ref_id;
////                sam_aln.sw_score = aln_params.match*(read_len-hamming_dist) - aln_params.mismatch*hamming_dist;
//            }
////        } else if ( (best_align_dist > 1) || aln_did_not_fit ){
        } else {
//         if ( (hamming_dist < 0 ) || (ref_segm.length() != read_len) || aln_did_not_fit || needs_aln ){
            int extra_ref_left = std::min(soft_left, 50);
            int extra_ref_right = std::min(soft_right, 50);
            //            extra_ref = 50; //(read_diff - ref_diff) > 0 ?  (read_diff - ref_diff) : 0;
            int a = n.ref_s - n.query_s - extra_ref_left;
            int ref_start = std::max(0, a);
            int b = n.ref_e + (read_len - n.query_e)+ extra_ref_right;
            int ref_len = ref_len_map[n.ref_id];
            int ref_end = std::min(ref_len, b);
            ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_end - ref_start);
//            ksw_extz_t ez;
//            const char *ref_ptr = ref_segm.c_str();
//            const char *read_ptr = r_tmp.c_str();
            aln_info info;
//            std::cout << "Extra ref: " << extra_ref << " " << read_diff << " " << ref_diff << " " << ref_start << " " << ref_end << std::endl;
//            info = ksw_align(ref_ptr, ref_segm.size(), read_ptr, r_tmp.size(), 1, 4, 6, 1, ez);
            info = ssw_align(ref_segm, r_tmp, aln_params.match, aln_params.mismatch, aln_params.gap_open, aln_params.gap_extend);
//            info.ed = info.global_ed; // read_len - info.sw_score;
            sw_score = info.sw_score;
            log_vars.tot_ksw_aligned ++;
            int diff_to_best = sw_score < best_align_sw_score ? best_align_sw_score - sw_score : sw_score - best_align_sw_score;
            min_mapq_diff = min_mapq_diff < diff_to_best ? min_mapq_diff : diff_to_best;
//            if (info.global_ed <= best_align_dist) {
//                min_mapq_diff = best_align_dist - info.global_ed; // new distance to next best match
//            }
            if (sw_score >= best_align_sw_score){
                min_mapq_diff = (sw_score - best_align_sw_score) > 0 ? (sw_score - best_align_sw_score)  : 0 ; // new distance to next best match
            }
            best_align_dist = info.global_ed;
            sam_aln.global_ed = info.global_ed;
            sam_aln.cigar = info.cigar;
            sam_aln.ed = info.ed;
            sam_aln.ref_start =  ref_start + info.ref_offset +1; // +1 because SAM is 1-based!
            sam_aln.is_rc = is_rc;
            sam_aln.ref_id = n.ref_id;
            sam_aln.sw_score = info.sw_score;
            sam_aln.aln_score = info.sw_score;
            sam_aln.aln_length = info.length;
//            std::cout << "Aligned: " << n.score << ", "  << n.n_hits << ", " << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e  << ") ed:" << info.ed << ", best ed so far: " << best_align_dist  << std::endl;

        }
        if (sw_score > best_align_sw_score){
            best_align_sw_score = sw_score;
        }

        std::tuple<int, alignment> t (sam_aln.sw_score, sam_aln);
        alignments.push_back(t);
        cnt ++;
    }

    //
    if (all_nams.size() > 0) {
        std::sort(alignments.begin(), alignments.end(), sort_highest_sw_scores_single); // Sorting by highest sw first
        int max_out = alignments.size() < max_secondary ? alignments.size() : max_secondary;
        for (int i = 0; i < max_out; ++i) {
            auto aln = alignments[i];
            auto sam_aln = std::get<1>(aln);
            if ((sam_aln.sw_score - best_align_sw_score) > (2*aln_params.mismatch + aln_params.gap_open) ){
                break;
            }
            bool is_secondary = false;

            if (i > 0) {
                is_secondary = true;
                sam_aln.mapq = 255;
            } else {
                sam_aln.mapq = std::min(min_mapq_diff, 60);
            }
            sam.add(sam_aln, record, read_rc, is_secondary);
        }
    }
}


static inline void align_segment(alignment_params &aln_params, std::string &read_segm, std::string &ref_segm, int read_segm_len, int ref_segm_len, int ref_start,  int ext_left, int ext_right, bool aln_did_not_fit, bool is_rc, alignment &sam_aln_segm, unsigned int &tot_ksw_aligned) {
    int hamming_dist = -1;
    int soft_left = 50;
    int soft_right = 50;
    int hamming_mod;
    int ref_segm_len_ham = ref_segm_len - ext_left - ext_right; // we send in the already extended ref segment to save time. This is not true in center alignment if merged match have diff length
//    std::cerr << "ref_segm_len_ham: " << ref_segm_len_ham << " read_segm_len: " << read_segm_len << std::endl;
    if ( (ref_segm_len_ham == read_segm_len) && (!aln_did_not_fit) ){
        std::string ref_segm_ham = ref_segm.substr(ext_left, read_segm_len);
//        std::cout << "ref_segm_ham " << ref_segm_ham << std::endl;

        hamming_dist = HammingDistance(read_segm, ref_segm_ham);
//        std::cout << "hamming_dist " << hamming_dist << std::endl;

        if ( (hamming_dist >= 0) && (((float) hamming_dist / read_segm_len) < 0.05) ) { //Hamming distance worked fine, no need to ksw align
            std::stringstream cigar_string;
            int aln_score = 0;
            hamming_mod = HammingToCigarEQX2(read_segm, ref_segm_ham, cigar_string, aln_params.match, aln_params.mismatch, aln_score, soft_left, soft_right);
            sam_aln_segm.cigar = cigar_string.str();
            sam_aln_segm.ed = hamming_mod;
            sam_aln_segm.sw_score = aln_score; // aln_params.match*(read_len-hamming_dist) - aln_params.mismatch*hamming_dist;
            sam_aln_segm.ref_start = ref_start + ext_left + soft_left+1; // +1 because SAM is 1-based!
            sam_aln_segm.is_rc = is_rc;
            sam_aln_segm.is_unaligned = false;
            sam_aln_segm.aln_score = aln_score;
            sam_aln_segm.aln_length = read_segm_len;
//            std::cerr << "HAMMING WORKED " << sam_aln_segm.cigar << " sam_aln_segm.ref_start " << sam_aln_segm.ref_start << std::endl;
            return;
        }
    }

    aln_info info;
    info = ssw_align(ref_segm, read_segm, aln_params.match, aln_params.mismatch, aln_params.gap_open, aln_params.gap_extend);
    tot_ksw_aligned ++;
    sam_aln_segm.cigar = info.cigar;
    sam_aln_segm.ed = info.ed;
//    std::cerr << r_tmp << " " << n.n_hits << " " << n.score << " " <<  diff << " " << sam_aln.ed << " "  <<  n.query_s << " "  << n.query_e << " "<<  n.ref_s << " "  << n.ref_e << " " << n.is_rc << " " << hamming_dist << " " << sam_aln.cigar << " " << info.sw_score << std::endl;
    sam_aln_segm.sw_score = info.sw_score;
    sam_aln_segm.ref_start =  ref_start + info.ref_offset +1; // +1 because SAM is 1-based!
    sam_aln_segm.is_rc = is_rc;
    sam_aln_segm.is_unaligned = false;
    sam_aln_segm.aln_score = info.sw_score;
    sam_aln_segm.aln_length = info.length;
//    std::cerr << " ALIGN SCORE: " << sam_aln_segm.sw_score << " cigar: " << sam_aln_segm.cigar << std::endl;
}

static inline void get_alignment(alignment_params &aln_params, nam &n, const std::vector<unsigned int> &ref_len_map, const std::vector<std::string> &ref_seqs, const std::string &read, std::string &read_rc, alignment &sam_aln, int k, bool &rc_already_comp, unsigned int &did_not_fit, unsigned int &tot_ksw_aligned){
    auto read_len = read.size();
    bool aln_did_not_fit = false;
    int ref_diff = n.ref_e - n.ref_s;
    int read_diff = n.query_e - n.query_s;
    int min_diff =  read_diff ^ ((ref_diff ^ read_diff) & -(ref_diff < read_diff));
    int max_diff = ref_diff ^ ((ref_diff ^ read_diff) & -(ref_diff < read_diff));
    int diff = max_diff - min_diff;
//    int max_allowed_mask = aln_params.gap_open/aln_params.match - 1 > 0 ? aln_params.gap_open/aln_params.match - 1 : 1;

//    std::cerr << "n.ID " << n.nam_id  << " n.n_hits " << n.n_hits << " n.ref_s " <<  n.ref_s <<  " n.ref_e " << n.ref_e << " read " << read << std::endl;

    // decide if read should be fw or rc aligned to reference here by checking exact match of first and last strobe in the NAM
    bool fits = false;
    std::string ref_start_kmer;
    std::string ref_end_kmer;
    std::string read_start_kmer;
    std::string read_end_kmer;
    std::string read_rc_start_kmer;
    std::string read_rc_end_kmer;
    ref_start_kmer = ref_seqs[n.ref_id].substr(n.ref_s, k);
    ref_end_kmer = ref_seqs[n.ref_id].substr(n.ref_e-k, k);

    if (!n.is_rc) {
        read_start_kmer = read.substr(n.query_s, k);
        read_end_kmer = read.substr(n.query_e-k, k);
        if ((ref_start_kmer == read_start_kmer) && (ref_end_kmer == read_end_kmer)) {
//            n.is_rc = false;
            fits = true;
        } else  {
            //  FALSE FORWARD TAKE CARE OF FALSE HITS HERE - it can be false forwards or false rc because of symmetrical hash values
            //    we need two extra checks for this - hopefully this will remove all the false hits we see (true hash collisions should be very few)

//              std::cerr << " CHECKING1!! " << std::endl;
            // false reverse hit, change coordinates in nam to forward
            if (!rc_already_comp){
                read_rc = reverse_complement(read);
                rc_already_comp = true;
            }

            int q_start_tmp = read_len - n.query_e;
            int q_end_tmp = read_len - n.query_s;
            read_start_kmer = read_rc.substr(q_start_tmp, k);
            read_end_kmer = read_rc.substr(q_end_tmp-k, k);
            if ((ref_start_kmer == read_start_kmer) && (ref_end_kmer == read_end_kmer)){
                fits = true;
                n.is_rc = true;
                n.query_s = q_start_tmp;
                n.query_e = q_end_tmp;
//                std::cerr << " DETECTED FALSE RC FROM SYMM!! " << std::endl;
            }

        }
    } else {
        if (!rc_already_comp){
            read_rc = reverse_complement(read);
            rc_already_comp = true;
        }
        read_rc_start_kmer = read_rc.substr(n.query_s, k);
        read_rc_end_kmer = read_rc.substr(n.query_e-k, k);
        if ( (ref_start_kmer == read_rc_start_kmer) && (ref_end_kmer == read_rc_end_kmer) ) { // && (ref_segm.substr(n.query_e - k + (ref_diff - read_diff), k) == read_rc.substr(n.query_e - k, k)) ){
            n.is_rc = true;
            fits = true;
        } else{
            //  FALSE REVERSE TAKE CARE OF FALSE HITS HERE - it can be false forwards or false rc because of symmetrical hash values
            //    we need two extra checks for this - hopefully this will remove all the false hits we see (true hash collisions should be very few)

            int q_start_tmp = read_len - n.query_e;
            int q_end_tmp = read_len - n.query_s;
            read_start_kmer = read.substr(q_start_tmp, k);
            read_end_kmer = read.substr(q_end_tmp-k, k);
//            std::cerr << " CHECKING2!! " <<   n.query_s << " " <<   n.query_e << " " << std::endl;
//            std::cerr << read_start_kmer  << " " <<  ref_start_kmer << " " <<  read_end_kmer << " " << ref_end_kmer << std::endl;

            if ((ref_start_kmer == read_start_kmer) && (ref_end_kmer == read_end_kmer)){
                fits = true;
                n.is_rc = false;
                n.query_s = q_start_tmp;
                n.query_e = q_end_tmp;
//                std::cerr << " DETECTED FALSE FW FROM SYMM!! " << std::endl;
            }
        }
    }

    if (!fits) {
        did_not_fit++;
        aln_did_not_fit = true;
        sam_aln.not_proper = true;
    }

    std::string r_tmp;
    bool is_rc;
    if (n.is_rc){
        r_tmp = read_rc;
        is_rc = true;
    }else{
        r_tmp = read;
        is_rc = false;
    }

    int ext_left;
    int ext_right;
    int ref_tmp_segm_size;
    int ref_len = ref_len_map[n.ref_id];
    size_t ref_segm_size;
    int ref_tmp_start;
    int ref_start;
    std::string ref_segm;
    std::string read_segm;

    if (true){ // currently deactivate partial SSW implementation.. (!fits){ // full alignment
        ref_tmp_start = n.ref_s - n.query_s > 0 ? n.ref_s - n.query_s : 0;
        ext_left = ref_tmp_start < 50 ? ref_tmp_start : 50;
        ref_start = ref_tmp_start - ext_left;

        ref_tmp_segm_size = read_len + diff;
        ext_right = ref_len - (n.ref_e +1) < 50 ? ref_len - (n.ref_e +1) : 50;

        ref_segm_size = ref_tmp_segm_size + ext_left + ext_right;
//        if (ref_start < 0  ){
//            std::cerr << "Get_alignment Bug1! ref start: " << ref_start << " ref_segm_size: " << ref_segm_size  << " ref len:  " << ref_seqs[n.ref_id].length() << std::endl;
//        }
//        if (ref_start + ref_segm_size >= ref_seqs[n.ref_id].length() ){
//            std::cerr << "Get_alignment Bug2! ref start: " << ref_start << " ref_segm_size: " << ref_segm_size << " ref len:  " << ref_seqs[n.ref_id].length() << std::endl;
//        }
//        if (ref_segm_size <= 20){
//            std::cerr << "Get_alignment Bug3! ref start: " << ref_start << " ref_segm_size: " << ref_segm_size << " ref len:  " << ref_seqs[n.ref_id].length() << std::endl;
//        }
        ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_segm_size);
//        std::cerr << " ref_tmp_start " << ref_tmp_start << " ext left " << ext_left << " ext right " << ext_right << " ref_tmp_segm_size " << ref_tmp_segm_size << " ref_segm_size " << ref_segm_size << " ref_segm " << ref_segm << std::endl;
        sam_aln.ref_id = n.ref_id;
        align_segment(aln_params, r_tmp, ref_segm, read_len, ref_segm_size, ref_start, ext_left, ext_right, aln_did_not_fit, is_rc, sam_aln, tot_ksw_aligned);
    } else{
        // test full hamming based alignment first
        ref_tmp_start = n.ref_s - n.query_s > 0 ? n.ref_s - n.query_s : 0;
        int ref_start = ref_tmp_start > 0 ? ref_tmp_start : 0;
        ref_tmp_segm_size = read_len + diff;
        ref_segm_size = ref_tmp_segm_size < ref_len - ref_start + 1 ? ref_tmp_segm_size : ref_len - ref_start + 1;

        std::string ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_segm_size);
        if ( (ref_segm_size == read_len) && fits ){

            int hamming_dist = HammingDistance(r_tmp, ref_segm);

            if ( (hamming_dist >= 0) && (((float) hamming_dist / ref_segm_size) < 0.05) ) { //Hamming distance worked fine, no need to ksw align
                std::stringstream cigar_string;
                int aln_score = 0;
                int soft_left = 0;
                int soft_right = 0;
                int hamming_mod = HammingToCigarEQX2(r_tmp, ref_segm, cigar_string, aln_params.match, aln_params.mismatch, aln_score, soft_left, soft_right);
                sam_aln.cigar = cigar_string.str();
                sam_aln.ed = hamming_mod;
                sam_aln.sw_score = aln_score; // aln_params.match*(read_len-hamming_dist) - aln_params.mismatch*hamming_dist;
                sam_aln.ref_start = ref_start + ext_left + soft_left+1; // +1 because SAM is 1-based!
                sam_aln.is_rc = is_rc;
                sam_aln.is_unaligned = false;
                sam_aln.aln_score = aln_score;
//                std::cerr << "FULL HAMMING , returning " << sam_aln.cigar << " sam_aln_segm.ref_start " << sam_aln.ref_start << std::endl;
                return;
            }
        }
//        std::cerr << " GOING HERE!!!!!" << std::endl;
//        std::cerr << " " << std::endl;

//        std::cerr << diff << " " << ref_segm_size << "  " << read_len << std::endl;

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
            int left_ref_end_bp = -1;
            int right_ref_start_bp = -1;
            if (break_point == n.query_s){
                left_ref_end_bp = n.ref_s +k;
                right_ref_start_bp = n.ref_s;
            } else if (break_point == (n.query_e - k)) {
                left_ref_end_bp = n.ref_e;
                right_ref_start_bp = n.ref_e-k;
            } else  {
                std::cerr << "BUUUUUUG " << std::endl;
            }


//            int ref_diff = n.ref_e - n.ref_s;
//            int read_diff = n.query_e - n.query_s;
//            int min_diff =  read_diff ^ ((ref_diff ^ read_diff) & -(ref_diff < read_diff));
//            int max_diff = ref_diff ^ ((ref_diff ^ read_diff) & -(ref_diff < read_diff));
//            int diff = max_diff - min_diff;

            // Left region align
            alignment sam_aln_segm_left;
            read_segm = r_tmp.substr(0, left_region_bp);
            ref_tmp_start = n.ref_s - n.query_s > 0 ? n.ref_s - n.query_s : 0;
            ext_left = ref_tmp_start < 50 ? ref_tmp_start : 50;
            ext_right = 0;
            ref_start = ref_tmp_start - ext_left;
            ref_segm_size = left_region_bp + ext_left + diff;
            ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_segm_size);
//            std::cerr << " "  << std::endl;
//            std::cerr << "GOING IN LEFT: " << " read segm len " << read_segm.length() << " ref segm len " << ref_segm_size  << " ext_left: " << ext_left << std::endl;
//            std::cerr << diff << std::endl;
//            std::cerr << read_segm << std::endl;
//            std::cerr << ref_segm << std::endl;

            align_segment(aln_params, read_segm, ref_segm, read_segm.length(), ref_segm_size, ref_start, ext_left, ext_right, aln_did_not_fit, is_rc, sam_aln_segm_left, tot_ksw_aligned);
//            std::cerr << "LEFT CIGAR: " << sam_aln_segm_left.cigar << std::endl;

            //Right region align
            alignment sam_aln_segm_right;
            read_segm = r_tmp.substr(right_region_bp, read_len - right_region_bp );
            ref_tmp_segm_size = right_ref_start_bp + (read_len + diff - right_region_bp) < ref_len ? (read_len + diff - right_region_bp) : ref_len - right_ref_start_bp;
            ext_left = 0;
            ext_right = std::min(ref_len - (right_ref_start_bp + ref_tmp_segm_size), 50);
            ref_segm_size = ref_tmp_segm_size + ext_right;
            ref_segm = ref_seqs[n.ref_id].substr(right_ref_start_bp, ref_segm_size);
//            std::cerr << " "  << std::endl;
//            std::cerr << "GOING IN RIGHT: " << " read segm len " << read_segm.length() << " ref segm len " << ref_segm_size  << " ext_right: " << ext_right << std::endl;
//            std::cerr << diff << std::endl;
//            std::cerr << read_segm << std::endl;
//            std::cerr << ref_segm << std::endl;
//            std::cerr << "read_segm.length(): " << read_segm.length() << " ref_segm_size " << ref_segm_size << std::endl;
            align_segment(aln_params, read_segm, ref_segm, read_segm.length(), ref_segm_size, ref_start, ext_left, ext_right, aln_did_not_fit, is_rc, sam_aln_segm_right, tot_ksw_aligned);
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
//            std::cerr << "JOINT CIGAR: " << sam_aln.cigar << std::endl;
//            std::cerr << " " << std::endl;


        } else{
//            std::cerr << "NOOOO MAX BREAKPOINT " << break_point << " candidates: "  <<  n.query_s  << " " << n.query_e - k << std::endl;
            // full align
            ref_tmp_start = n.ref_s - n.query_s > 0 ? n.ref_s - n.query_s : 0;
            ext_left = ref_tmp_start < 50 ? ref_tmp_start : 50;
            ref_start = ref_tmp_start - ext_left;

            ref_tmp_segm_size = read_len + diff;
            ext_right = ref_len - (n.ref_e +1) < 50 ? ref_len - (n.ref_e +1) : 50;

            ref_segm_size = ref_tmp_segm_size + ext_left + ext_right;
            ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_segm_size);
//        std::cerr << " ref_tmp_start " << ref_tmp_start << " ext left " << ext_left << " ext right " << ext_right << " ref_tmp_segm_size " << ref_tmp_segm_size << " ref_segm_size " << ref_segm_size << " ref_segm " << ref_segm << std::endl;
            sam_aln.ref_id = n.ref_id;
            align_segment(aln_params, r_tmp, ref_segm, read_len, ref_segm_size, ref_start, ext_left, ext_right, aln_did_not_fit, is_rc, sam_aln, tot_ksw_aligned);
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

    }

}


//
//
//static inline void get_alignment(alignment_params &aln_params, nam &n, std::vector<unsigned int> &ref_len_map, std::vector<std::string> &ref_seqs, std::string &read, std::string &read_rc, int read_len, alignment &sam_aln, int k, int cnt, bool &rc_already_comp, unsigned int &did_not_fit, unsigned int &tot_ksw_aligned){
//    bool aln_did_not_fit = false;
//    int ref_diff = n.ref_e - n.ref_s;
//    int read_diff = n.query_e - n.query_s;
//    int min_diff =  read_diff ^ ((ref_diff ^ read_diff) & -(ref_diff < read_diff));
//    int max_diff = ref_diff ^ ((ref_diff ^ read_diff) & -(ref_diff < read_diff));
//    int diff = max_diff - min_diff;
////    int max_allowed_mask = aln_params.gap_open/aln_params.match - 1 > 0 ? aln_params.gap_open/aln_params.match - 1 : 1;
//
//    // deal with any read hanging of ends of reference not to get 'std::out_of_range' what(): basic_string::substr
//    int ref_tmp_start = n.ref_s - n.query_s;
//    int ref_tmp_segm_size = read_len + diff;
//    int ref_len = ref_len_map[n.ref_id];
//    int ref_start = ref_tmp_start > 0 ? ref_tmp_start : 0;
//    int ref_segm_size = ref_tmp_segm_size < ref_len - ref_start ? ref_tmp_segm_size : ref_len - 1 - ref_start;
//
//    std::string ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_segm_size);
//
//    // decide if read should be fw or rc aligned to reference here by checking exact match of first and last strobe in the NAM
//    bool fits = false;
//    std::string ref_start_kmer;
//    std::string ref_end_kmer;
//    std::string read_start_kmer;
//    std::string read_end_kmer;
//    std::string read_rc_start_kmer;
//    std::string read_rc_end_kmer;
//    ref_start_kmer = ref_seqs[n.ref_id].substr(n.ref_s, k);
//    ref_end_kmer = ref_seqs[n.ref_id].substr(n.ref_e-k, k);
//
//    if (!n.is_rc) {
//        read_start_kmer = read.substr(n.query_s, k);
//        read_end_kmer = read.substr(n.query_e-k, k);
//        if ((ref_start_kmer == read_start_kmer) && (ref_end_kmer == read_end_kmer)) {
////            n.is_rc = false;
//            fits = true;
//        } else  {
//            //  FALSE FORWARD TAKE CARE OF FALSE HITS HERE - it can be false forwards or false rc because of symmetrical hash values
//            //    we need two extra checks for this - hopefully this will remove all the false hits we see (true hash collisions should be very few)
//
////              std::cerr << " CHECKING1!! " << std::endl;
//            // false reverse hit, change coordinates in nam to forward
//            if (!rc_already_comp){
//                read_rc = reverse_complement(read);
//                rc_already_comp = true;
//            }
//
//            int q_start_tmp = read_len - n.query_e;
//            int q_end_tmp = read_len - n.query_s;
//            read_start_kmer = read_rc.substr(q_start_tmp, k);
//            read_end_kmer = read_rc.substr(q_end_tmp-k, k);
//            if ((ref_start_kmer == read_start_kmer) && (ref_end_kmer == read_end_kmer)){
//                fits = true;
//                n.is_rc = true;
//                n.query_s = q_start_tmp;
//                n.query_e = q_end_tmp;
////                std::cerr << " DETECTED FALSE RC FROM SYMM!! " << std::endl;
//            }
//
//        }
//    } else {
//        if (!rc_already_comp){
//            read_rc = reverse_complement(read);
//            rc_already_comp = true;
//        }
//        read_rc_start_kmer = read_rc.substr(n.query_s, k);
//        read_rc_end_kmer = read_rc.substr(n.query_e-k, k);
//        if ( (ref_start_kmer == read_rc_start_kmer) && (ref_end_kmer == read_rc_end_kmer) ) { // && (ref_segm.substr(n.query_e - k + (ref_diff - read_diff), k) == read_rc.substr(n.query_e - k, k)) ){
//            n.is_rc = true;
//            fits = true;
//        } else{
//            //  FALSE REVERSE TAKE CARE OF FALSE HITS HERE - it can be false forwards or false rc because of symmetrical hash values
//            //    we need two extra checks for this - hopefully this will remove all the false hits we see (true hash collisions should be very few)
//
//            int q_start_tmp = read_len - n.query_e;
//            int q_end_tmp = read_len - n.query_s;
//            read_start_kmer = read.substr(q_start_tmp, k);
//            read_end_kmer = read.substr(q_end_tmp-k, k);
////            std::cerr << " CHECKING2!! " <<   n.query_s << " " <<   n.query_e << " " << std::endl;
////            std::cerr << read_start_kmer  << " " <<  ref_start_kmer << " " <<  read_end_kmer << " " << ref_end_kmer << std::endl;
//
//            if ((ref_start_kmer == read_start_kmer) && (ref_end_kmer == read_end_kmer)){
//                fits = true;
//                n.is_rc = false;
//                n.query_s = q_start_tmp;
//                n.query_e = q_end_tmp;
////                std::cerr << " DETECTED FALSE FW FROM SYMM!! " << std::endl;
//            }
//        }
//    }
//
//    if (!fits) {
//        did_not_fit++;
//        aln_did_not_fit = true;
//        sam_aln.not_proper = true;
//    }
//
//    int hamming_dist = -1;
//    std::string r_tmp;
//    bool is_rc;
//    if (n.is_rc){
//        r_tmp = read_rc;
//        is_rc = true;
//    }else{
//        r_tmp = read;
//        is_rc = false;
//    }
//
////    std::cerr<< r_tmp << std::endl;
////    std::cerr<< ref_segm << std::endl;
////    std::cerr<< diff << std::endl;
//    int soft_left = 50;
//    int soft_right = 50;
//    int hamming_mod;
////    bool needs_aln = false;
//    if ( (ref_segm_size == read_len) && (!aln_did_not_fit) ){
//        hamming_dist = HammingDistance(r_tmp, ref_segm);
////        std::cerr<< "Here " << hamming_dist << std::endl;
////        std::cerr<< aln_params.gap_open/aln_params.match  << std::endl;
//        if ( (hamming_dist >= 0) && (((float) hamming_dist / read_len) < 0.05) ) { //Hamming distance worked fine, no need to ksw align
//            std::stringstream cigar_string;
////            needs_aln = HammingToCigarEQX(r_tmp, ref_segm, cigar_string);
//            int aln_score = 0;
//            hamming_mod = HammingToCigarEQX2(r_tmp, ref_segm, cigar_string, aln_params.match, aln_params.mismatch, aln_score, soft_left, soft_right);
//
////            needs_aln = false;
//            sam_aln.cigar = cigar_string.str();
////            sam_aln.cigar = std::to_string(read_len) + "M";
////            std::cerr<< "Here ham dist: " << hamming_dist << " ham mod: " << hamming_mod << " " << r_tmp.size() << " " << ref_segm.size()  << std::endl;
//            sam_aln.ed = hamming_mod;
////            sam_aln.sw_score = aln_score;
//            sam_aln.sw_score = aln_score; // aln_params.match*(read_len-hamming_dist) - aln_params.mismatch*hamming_dist;
//            sam_aln.ref_start = ref_start + soft_left+1; // +1 because SAM is 1-based!
//            sam_aln.is_rc = is_rc;
//            sam_aln.ref_id = n.ref_id;
//            sam_aln.is_unaligned = false;
//            sam_aln.aln_score = aln_score;
//            return;
////            if (hamming_mod == hamming_dist ){ // masked only what is justified by alingment parameters max_allowed_mask
////                return;
////            }
//        }
//        //TODO: Only do ksw of the ends outside the NAM to increase speed here
////        else{ // Segment(s) of read outside the NAM span is not fitting to reference, align the segments
////            std::cerr<< sam_aln.ed << " " << sam_aln.sw_score << " " <<   n.query_s << " " << n.query_e << std::endl;
////            std::cerr<< r_tmp << std::endl;
////            std::cerr<< ref_segm.substr(0,read_len) << std::endl;
////
////        }
//    }
//
//    // We didn't get away with hamming distance, do full ksw alignment
////    else {
////    std::cerr<< "3" << std::endl;
//
//    int extra_ref_left = soft_left <= 50 ? soft_left : 50;
//    int extra_ref_right = soft_right <= 50 ? soft_right: 50;
//    int a = n.ref_s - n.query_s - extra_ref_left;
//    ref_start = std::max(0, a);
//    int b = n.ref_e + (read_len - n.query_e)+ extra_ref_right;
//    int ref_end = std::min(ref_len, b);
//    ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_end - ref_start);
////    ksw_extz_t ez;
////    const char *ref_ptr = ref_segm.c_str();
////    const char *read_ptr = r_tmp.c_str();
//    aln_info info;
////    std::cerr<< "4" << std::endl;
////    info = ksw_align(ref_ptr, ref_segm.size(), read_ptr, r_tmp.size(), 1, 4, 6, 1, ez);
//    info = ssw_align(ref_segm, r_tmp, read_len, aln_params.match, aln_params.mismatch, aln_params.gap_open, aln_params.gap_extend);
////    if (info.ed == 100000){
////        std::cerr<< "________________________________________" << std::endl;
////        std::cerr<< "NORMAL MODE" << std::endl;
////        std::cerr<< read << "   " << read_rc << std::endl;
////        std::cerr << r_tmp << " " << n.n_hits << " " << n.score << " " <<  diff << " " << sam_aln.ed << " "  <<  n.query_s << " "  << n.query_e << " "<<  n.ref_s << " "  << n.ref_e << " " << n.is_rc << " " << hamming_dist << " " << sam_aln.cigar << " " << info.sw_score << std::endl;
////        std::cerr << "a " << a << " b " << b << " ref_start " <<  ref_start << " ref_end " << ref_end << "  ref_end - ref_start "  <<  ref_end - ref_start << " extra_ref_left "  << extra_ref_left << " extra_ref_right "<<  extra_ref_right << "  n.is_flipped " <<  n.is_flipped << std::endl;
////        std::cerr<< "________________________________________" << std::endl;
////    }
//
////    std::cerr<< "5" << std::endl;
//    sam_aln.cigar = info.cigar;
//    sam_aln.ed = info.ed;
////    std::cerr << r_tmp << " " << n.n_hits << " " << n.score << " " <<  diff << " " << sam_aln.ed << " "  <<  n.query_s << " "  << n.query_e << " "<<  n.ref_s << " "  << n.ref_e << " " << n.is_rc << " " << hamming_dist << " " << sam_aln.cigar << " " << info.sw_score << std::endl;
//    sam_aln.sw_score = info.sw_score;
//    sam_aln.ref_start =  ref_start + info.ref_offset +1; // +1 because SAM is 1-based!
//    sam_aln.is_rc = is_rc;
//    sam_aln.ref_id = n.ref_id;
//    sam_aln.is_unaligned = false;
//    sam_aln.aln_score = info.sw_score;
//    tot_ksw_aligned ++;
////    }
//}

static inline void get_MAPQ(std::vector<nam> &all_nams, nam &n_max, int &mapq){
    float s1 = n_max.score;
    mapq = 60; // MAPQ = 40(1−s2/s1) ·min{1,|M|/10} · log s1
    if (all_nams.size() > 1) {
        nam n_second = all_nams[1];
        float s2 = n_second.score;
        float min_matches;
        min_matches  = (float)n_max.n_hits/10 > 1 ? (float)n_max.n_hits/10 : 1;
        mapq = 40*(1 - s2/s1)*min_matches*log(s1) < 60 ? 40*(1 - s2/s1)*min_matches*log(s1) : 60 ;
    }
}

static inline void get_joint_MAPQ(float s1, float s2, int joint_n_matches, int &mapq1, int &mapq2){
    // MAPQ = 40(1−s2/s1) ·min{1,|M|/10} · log s1     // from minimap paper
    float min_matches;
    min_matches  = (float)joint_n_matches/10 > 1 ? (float)joint_n_matches/10 : 1;
    mapq1 = 40*(1 - s2/s1)*min_matches*log(s1) < 60 ? 40*(1 - s2/s1)*min_matches*log(s1) : 60 ;
    if (mapq1 < 0){
        mapq1 = 0;
    }
    mapq2 = mapq1;
}

static inline void get_joint_MAPQ_from_alingments(float S1, float S2, int &mapq1, int &mapq2){
    if (S1 == S2){ // At least two identical placements
        mapq1 = 0;
        mapq2 = mapq1;
    } else {
        int diff = S1 - S2; // (1.0 - (S1 - S2) / S1);
//        float log10_p = diff > 6 ? -6.0 : -diff; // Corresponds to: p_error= 0.1^diff // change in sw score times rough illumina error rate. This is highly heauristic, but so seem most computations of mapq scores
        if ((S1 > 0) && (S2 > 0)) {
            mapq1 = diff <= 60 ? diff : 60;
//            mapq1 = -10 * log10_p < 60 ? -10 * log10_p : 60;
            mapq2 = mapq1;
        } else if ((S1 > 0) && (S2 <= 0)) {
            mapq1 = 60;
            mapq2 = mapq1;
        } else { // both negative SW one is better
            mapq1 = 1;
            mapq2 = mapq1;
        }
    }
}


static inline float normal_pdf(float x, float m, float s)
{
    static const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

static inline bool score_sw(const alignment &a, const alignment &b)
{
    return ( a.sw_score > b.sw_score );
}

static inline bool sort_scores(const std::tuple<double, alignment, alignment> &a,
                               const std::tuple<double, alignment, alignment> &b)
{
    return (std::get<0>(a) > std::get<0>(b));
}

static inline bool sort_joint_hits(const std::tuple<int, nam, nam> &a,
                                   const std::tuple<int, nam, nam> &b)
{
    return (std::get<0>(a) > std::get<0>(b));
}




static inline void get_best_scoring_pair(std::vector<alignment> &aln_scores1, std::vector<alignment> &aln_scores2, std::vector<std::tuple<double,alignment,alignment>> &high_scores, float mu, float sigma)
{
    double S;
    float x;
//            std::cerr << "Scoring" << std::endl;
    for (auto &a1 : aln_scores1) {
        for (auto &a2 : aln_scores2) {
            x = a1.ref_start > a2.ref_start ? (float) (a1.ref_start - a2.ref_start) : (float)(a2.ref_start - a1.ref_start);
//                    std::cerr << x << " " << (a1.ref_start - a2.ref_start) << " " << (a2.ref_start - a1.ref_start) << std::endl;
            if ( (a1.is_rc ^ a2.is_rc) && (x < mu+4*sigma)  ){
                // r1.sw_score + r2.sw_score - log P(d(r1,r2)) if -log P(d(r1,r2)) < 3,

                S = (double)a1.sw_score + (double)a2.sw_score + log( normal_pdf(x, mu, sigma ) );  //* (1 - s2 / s1) * min_matches * log(s1);
//                        std::cerr << S << " " << x << " " << log(normal_pdf(x, mu, sigma )) << " " << normal_pdf(x, mu, sigma ) << std::endl;
                std::tuple<double, alignment, alignment> t (S, a1, a2);
                high_scores.push_back(t);
            }
            else{ // individual score
                S = (double)a1.sw_score + (double)a2.sw_score - 10; // 10 corresponds to  a value of log( normal_pdf(x, mu, sigma ) ) of more than 4 stddevs away
//                        std::cerr << S << " individual score " << x << " " << std::endl;
                std::tuple<double, alignment, alignment> t (S, a1, a2);
                high_scores.push_back(t);
            }
        }
    }
    std::sort(high_scores.begin(), high_scores.end(), sort_scores); // Sorting by highest score first
}

static inline void get_best_scoring_NAM_locations(std::vector<nam> &all_nams1, std::vector<nam> &all_nams2, std::vector<std::tuple<int,nam,nam>> &joint_NAM_scores, float mu, float sigma, robin_hood::unordered_set<int> &added_n1, robin_hood::unordered_set<int> &added_n2)
{
    if ( all_nams1.empty() && all_nams2.empty() ){
        return;
    }
    int joint_hits;
    nam n;  //dummy nam
    n.ref_s = -1;
    int hjss = 0; // highest joint score seen
//            std::cerr << "Scoring" << std::endl;
    int a,b;
    for (auto &n1 : all_nams1) {
        for (auto &n2 : all_nams2) {
            if ((n1.n_hits + n2.n_hits) < hjss/2){
                break;
            }

            if ( (n1.is_rc ^ n2.is_rc) && (n1.ref_id == n2.ref_id) ){
                a = n1.ref_s - n1.query_s  > 0 ? n1.ref_s - n1.query_s : 0;
                b = n2.ref_s - n2.query_s  > 0 ? n2.ref_s - n2.query_s : 0;
//                std::cerr << "LEWL " <<  a << " " << b << std::endl;
                bool r1_r2 = n2.is_rc && (a < b) && ((b-a) < mu+10*sigma); // r1 ---> <---- r2
                bool r2_r1 = n1.is_rc && (b < a) && ((a-b) < mu+10*sigma); // r2 ---> <---- r1
//                std::cerr << r1_r2 << " " << r2_r1 << std::endl;

                if ( r1_r2 || r2_r1 ){
//                    int diff1 = (n1.query_e - n1.query_s) - (n1.ref_e - n1.ref_s);
//                    int  n1_penalty = diff1 > 0 ? diff1 : - diff1;
//                    int diff2 = (n2.query_e - n2.query_s) - (n2.ref_e - n2.ref_s);
//                    int  n2_penalty = diff2 > 0 ? diff2 : - diff2;
                    joint_hits = n1.n_hits + n2.n_hits; // - n1_penalty - n2_penalty; // trying out idea about penalty but it needs to be on the individual seed level - to late on merged match level.
                    std::tuple<int, nam, nam> t (joint_hits, n1, n2);
                    joint_NAM_scores.push_back(t);
                    added_n1.insert(n1.nam_id);
                    added_n2.insert(n2.nam_id);
                    if (joint_hits > hjss) {
                        hjss = joint_hits;
                    }
                }
            }

//            x = n1.ref_s > n2.ref_s ? (float) (n1.ref_s - n2.ref_s) : (float)(n2.ref_s - n1.ref_s);
////                    std::cerr << x << " " << (n1.ref_s - n2.ref_s) << " " << (n2.ref_s - n1.ref_s) << std::endl;
//            if ( (n1.is_rc ^ n2.is_rc) && (x < mu+10*sigma) && (n1.ref_id == n2.ref_id) ){
//                joint_hits = n1.n_hits + n2.n_hits;
//
////                        std::cerr << S << " " << x << " " << log(normal_pdf(x, mu, sigma )) << " " << normal_pdf(x, mu, sigma ) << std::endl;
//                std::tuple<int, nam, nam> t (joint_hits, n1, n2);
//                joint_NAM_scores.push_back(t);
//                added_n1.insert(n1.ref_s);
//                added_n2.insert(n2.ref_s);
//                if (joint_hits > hjss) {
//                    hjss = joint_hits;
//                }
//            }
        }
    }
//    std::cerr << "ADDED " << added_n1.size() << " " <<  added_n2.size() << std::endl;
//    for (auto z : added_n1){
//        std::cerr << z  << std::endl;
//    }

    if ( !all_nams1.empty() ){
        int hjss1 = hjss > 0 ? hjss : all_nams1[0].n_hits;
//        std::cerr <<" hjss1 " << hjss1 << " " << std::endl;
        //    int hjss1 = all_nams1[0].n_hits;
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
            //                        std::cerr << S << " individual score " << x << " " << std::endl;
            std::tuple<int, nam, nam> t (joint_hits, n1, n);
            joint_NAM_scores.push_back(t);
        }
    }

    if ( !all_nams2.empty() ){
        int hjss2 = hjss  > 0 ? hjss : all_nams2[0].n_hits;
//        std::cerr <<" hjss2 " << hjss2 << " " << std::endl;
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
            std::tuple<int, nam, nam> t (joint_hits, n, n2);
            joint_NAM_scores.push_back(t);
        }
    }

//    std::cerr << " All scores " << joint_NAM_scores.size() << std::endl;
    added_n1.clear();
    added_n2.clear();
    std::sort(joint_NAM_scores.begin(), joint_NAM_scores.end(), sort_joint_hits); // Sorting by highest score first

//    for (auto zz : joint_NAM_scores){
//        auto score_ = std::get<0>(zz);
//        auto n1_tmp = std::get<1>(zz);
//        auto n2_tmp = std::get<2>(zz);
//        std::cerr << "joint_NAM_score: " << score_ << " NAM ids: "  << n1_tmp.nam_id  << " " << n2_tmp.nam_id  << " Is RC: "  << n1_tmp.is_rc  << " " << n2_tmp.is_rc  << " Nhits: " << n1_tmp.n_hits  << " " << n2_tmp.n_hits  << " scores: " << n1_tmp.score  << " " << n2_tmp.score  << " ref_starts:" << n1_tmp.ref_s  << " " << n2_tmp.ref_s << " query_starts:" << n1_tmp.query_s  << " " << n2_tmp.query_e  << std::endl;
//    }
}


static inline void rescue_mate(alignment_params &aln_params , nam &n, const std::vector<unsigned int> &ref_len_map, const std::vector<std::string> &ref_seqs, std::string &guide_read, std::string &guide_read_rc, bool &guide_rc_already_comp, std::string &read, std::string &read_rc, alignment &sam_aln, bool &rc_already_comp, unsigned int &tot_ksw_aligned, float &mu, float &sigma, unsigned int &tot_rescued, int k
) {
    int a, b, ref_start,ref_len,ref_end;
    std::string r_tmp;
    bool a_is_rc;
    auto guide_read_len = guide_read.size();
    auto read_len = read.size();

    // decide if read should be fw or rc aligned to reference here by checking exact match of first and last strobe in the NAM
    bool fits = false;
    std::string ref_start_kmer;
    std::string ref_end_kmer;
    std::string read_start_kmer;
    std::string read_end_kmer;
    std::string read_rc_start_kmer;
    std::string read_rc_end_kmer;
//    std::cerr << "I'm here: n.ref_s:" << n.ref_s << " n.ref_e-k: " << n.ref_e-k <<  " ref_len: " << ref_seqs[n.ref_id].length() << " n.ref_e-k+k: " << n.ref_e-k+k << std::endl;

    ref_start_kmer = ref_seqs[n.ref_id].substr(n.ref_s, k);
    ref_end_kmer = ref_seqs[n.ref_id].substr(n.ref_e-k, k);

    if (!n.is_rc) {
//        std::cerr << "CHECK1: n.query_s:" << n.query_s << " n.query_e-k: " << n.query_e-k <<  " read_len: " << guide_read.length() << " n.query_e-k+k: " << n.query_e-k+k << std::endl;

        read_start_kmer = guide_read.substr(n.query_s, k);
        read_end_kmer = guide_read.substr(n.query_e-k, k);
        if ((ref_start_kmer == read_start_kmer) && (ref_end_kmer == read_end_kmer)) {
//            n.is_rc = false;
            fits = true;
        } else  {
            //  FALSE FORWARD TAKE CARE OF FALSE HITS HERE - it can be false forwards or false rc because of symmetrical hash values
            //    we need two extra checks for this - hopefully this will remove all the false hits we see (true hash collisions should be very few)

//              std::cerr << " CHECKING1!! " << std::endl;
//            std::cerr << read_start_kmer  << " " <<  ref_start_kmer << " " <<  read_end_kmer << " " << ref_end_kmer << std::endl;
            // false reverse hit, change coordinates in nam to forward
            if (!guide_rc_already_comp){
                guide_read_rc = reverse_complement(guide_read);
                guide_rc_already_comp = true;
            }

            int q_start_tmp = guide_read_len - n.query_e;
            int q_end_tmp = guide_read_len - n.query_s;
            read_start_kmer = guide_read_rc.substr(q_start_tmp, k);
            read_end_kmer = guide_read_rc.substr(q_end_tmp-k, k);
            if ((ref_start_kmer == read_start_kmer) && (ref_end_kmer == read_end_kmer)){
                fits = true;
                n.is_rc = true;
                n.query_s = q_start_tmp;
                n.query_e = q_end_tmp;
//                std::cerr << " RESCUE DETECTED FALSE RC FROM SYMM!! " << std::endl;
            }

        }
    } else {
        if (!guide_rc_already_comp){
            guide_read_rc = reverse_complement(guide_read);
            guide_rc_already_comp = true;
        }
//        std::cerr << "CHECK2: n.query_s:" << n.query_s << " n.query_e-k: " << n.query_e-k <<  " read_len: " << guide_read.length() << " n.query_e-k+k: " << n.query_e-k+k << std::endl;

        read_rc_start_kmer = guide_read_rc.substr(n.query_s, k);
        read_rc_end_kmer = guide_read_rc.substr(n.query_e-k, k);
        if ( (ref_start_kmer == read_rc_start_kmer) && (ref_end_kmer == read_rc_end_kmer) ) { // && (ref_segm.substr(n.query_e - k + (ref_diff - read_diff), k) == read_rc.substr(n.query_e - k, k)) ){
            n.is_rc = true;
            fits = true;
        } else{
            //  FALSE REVERSE TAKE CARE OF FALSE HITS HERE - it can be false forwards or false rc because of symmetrical hash values
            //    we need two extra checks for this - hopefully this will remove all the false hits we see (true hash collisions should be very few)

            int q_start_tmp = guide_read_len - n.query_e;
            int q_end_tmp = guide_read_len - n.query_s;
            read_start_kmer = guide_read.substr(q_start_tmp, k);
            read_end_kmer = guide_read.substr(q_end_tmp-k, k);
//            std::cerr << " CHECKING2!! " <<   n.query_s << " " <<   n.query_e << " " << std::endl;
//            std::cerr << read_start_kmer  << " " <<  ref_start_kmer << " " <<  read_end_kmer << " " << ref_end_kmer << std::endl;

            if ((ref_start_kmer == read_start_kmer) && (ref_end_kmer == read_end_kmer)){
                fits = true;
                n.is_rc = false;
                n.query_s = q_start_tmp;
                n.query_e = q_end_tmp;
//                std::cerr << " RESCUE DETECTED FALSE FW FROM SYMM!! " << std::endl;
            }
        }
    }

    if (n.is_rc){
        r_tmp = read;
        a = n.ref_s - n.query_s - (mu+5*sigma);
        b = n.ref_s - n.query_s + read_len/2; // at most half read overlap
        a_is_rc = false;
    }else{
        if (!rc_already_comp){
            read_rc = reverse_complement(read);
            rc_already_comp = true;
        }
        r_tmp = read_rc; // mate is rc since fr orientation
        a = n.ref_e + (read_len - n.query_e) - read_len/2; // at most half read overlap
        b = n.ref_e + (read_len - n.query_e) + (mu+5*sigma);
        a_is_rc = true;
    }

    ref_len = ref_len_map[n.ref_id];
    ref_start = std::max(0, std::min(a,ref_len));
    ref_end = std::min(ref_len, std::max(0, b));
//    if (ref_start == ref_len ){
//        std::cerr << "Rescue Bug1! ref start: " << ref_start << " ref end: " << ref_end << " ref len:  " << ref_len << std::endl;
//    }
//    if (ref_end == ref_len ){
//        std::cerr << "Rescue Bug2! ref start: " << ref_start << " ref end: " << ref_end << " ref len:  " << ref_len << std::endl;
//    }
//    if (ref_end <= ref_start){
//        std::cerr << "Rescue Bug3! ref start: " << ref_start << " ref end: " << ref_end << " ref len:  " << ref_len << std::endl;
//    }

    if (ref_end < ref_start + k){
        sam_aln.cigar = "*";
        sam_aln.ed = read_len;
        sam_aln.sw_score = 0;
        sam_aln.aln_score = 0;
        sam_aln.ref_start =  0;
        sam_aln.is_rc = n.is_rc;
        sam_aln.ref_id = n.ref_id;
        sam_aln.is_unaligned = true;
        sam_aln.not_proper = true;
//        std::cerr << "RESCUE: Caught Bug3! ref start: " << ref_start << " ref end: " << ref_end << " ref len:  " << ref_len << std::endl;
        return;
    }
    std::string ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_end - ref_start);
    aln_info info;

    // check that read shares at least some segment with ref otherwise abort
    int sub_size = 2*k/3;
    int step_size = k/3;
    std::string submer;
    bool found = false;
    for (size_t i = 0; i <= r_tmp.size()-k; i+=step_size) {
        submer = r_tmp.substr(i, sub_size);
        if (ref_segm.find( submer ) != std::string::npos) {
            found = true;
            break;
        }
    }
    if (!found){
        sam_aln.cigar = "*";
        sam_aln.ed = read_len;
        sam_aln.sw_score = 0;
        sam_aln.aln_score = 0;
        sam_aln.ref_start =  0;
        sam_aln.is_rc = n.is_rc;
        sam_aln.ref_id = n.ref_id;
        sam_aln.is_unaligned = true;
        sam_aln.not_proper = true;
//        std::cerr << "Avoided!" << std::endl;
        return;
//        std::cerr << "Aligning anyway at: " << ref_start << " to " << ref_end << "ref len:" << ref_len << " ref_id:" << n.ref_id << std::endl;
    }

    info = ssw_align(ref_segm, r_tmp, aln_params.match, aln_params.mismatch, aln_params.gap_open, aln_params.gap_extend);

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

    sam_aln.cigar = info.cigar;
    sam_aln.ed = info.ed;
    sam_aln.sw_score = info.sw_score;
    sam_aln.aln_score = sam_aln.sw_score;
    sam_aln.ref_start =  ref_start + info.ref_offset +1; // +1 because SAM is 1-based!
    sam_aln.is_rc = a_is_rc;
    sam_aln.ref_id = n.ref_id;
    sam_aln.is_unaligned = false;
    sam_aln.aln_length = info.length;
    tot_ksw_aligned ++;
    tot_rescued ++;
}


inline void align_PE(alignment_params &aln_params, Sam &sam, std::vector<nam> &all_nams1, std::vector<nam> &all_nams2,
                     const KSeq &record1, const KSeq &record2, int k, std::vector<unsigned int> &ref_len_map, const std::vector<std::string> &ref_seqs,
                     logging_variables &log_vars, float dropoff, i_dist_est &isize_est, int max_tries, size_t max_secondary) {
    auto mu = isize_est.mu;
    auto sigma = isize_est.sigma;
    std::string read1 = record1.seq;
    std::string read2 = record2.seq;
    std::string query_acc1 = record1.name;
    std::string query_acc2 = record2.name;
    double secondary_droppoff = 2 * aln_params.mismatch + aln_params.gap_open;
    std::string read1_rc;
    std::string read2_rc;
    bool rc_already_comp1 = false;
    bool rc_already_comp2 = false;
    int cnt = 0;
    double S = 0.0;
    float x;
    int cnt1 = 0;
    int cnt2 = 0;
    float score_dropoff1;
    float score_dropoff2;
    robin_hood::unordered_set<int> added_n1;
    robin_hood::unordered_set<int> added_n2;
//    int best_align_dist1 = ~0U >> 1;
//    int best_align_dist2 = ~0U >> 1;
    alignment sam_aln1;
    alignment sam_aln2;
    nam n_max1;
    nam n_max2;
    int mapq1, mapq2;

    if ((all_nams1.size() == 0) && (all_nams2.size() == 0)) { // None of the read pairs has any NAMs
        sam.add_unmapped_pair(record1, record2);
        return;
    } else if ((all_nams1.size() > 0) && (all_nams2.size() > 0)){
        n_max1 = all_nams1[0];
        n_max2 = all_nams2[0];
        score_dropoff1 = all_nams1.size() > 1 ? (float) all_nams1[1].n_hits / n_max1.n_hits : 0.0;
        score_dropoff2 = all_nams2.size() > 1 ? (float) all_nams2[1].n_hits / n_max2.n_hits : 0.0;
        score_dropoff1 = n_max1.n_hits > 2 ? score_dropoff1 : 1.0;
        score_dropoff2 = n_max2.n_hits > 2 ? score_dropoff2 : 1.0;
//        std::cerr << all_nams1.size() << " " << all_nams2.size() << " " << n_max1.n_hits << " " << n_max2.n_hits << " " << score_dropoff1 << " " << score_dropoff2 << " " << std::endl;

        int a = n_max1.ref_s - n_max1.query_s  > 0 ? n_max1.ref_s - n_max1.query_s : 0;
        int b = n_max2.ref_s - n_max2.query_s  > 0 ? n_max2.ref_s - n_max2.query_s : 0;
        bool r1_r2 = n_max2.is_rc && (a < b) && ((b-a) < 2000); // r1 ---> <---- r2
        bool r2_r1 = n_max1.is_rc && (b < a) && ((a-b) < 2000); // r2 ---> <---- r1
        if ( (score_dropoff1 < dropoff) && (score_dropoff2 < dropoff) && (n_max1.is_rc ^ n_max2.is_rc) && ( r1_r2 || r2_r1 ) ){ //( ((n_max1.ref_s - n_max2.ref_s) < mu + 4*sigma ) || ((n_max2.ref_s - n_max1.ref_s ) < mu + 4*sigma ) ) &&
//            std::cerr << "I'm here" << std::endl;
//            std::cerr << query_acc1 << std::endl;
            get_alignment(aln_params, n_max1, ref_len_map, ref_seqs, read1, read1_rc, sam_aln1, k, rc_already_comp1, log_vars.did_not_fit, log_vars.tot_ksw_aligned);
            log_vars.tot_all_tried ++;
//            std::cerr << query_acc2 << std::endl;
            get_alignment(aln_params, n_max2, ref_len_map, ref_seqs, read2, read2_rc, sam_aln2, k, rc_already_comp2, log_vars.did_not_fit, log_vars.tot_ksw_aligned);
            log_vars.tot_all_tried ++;
//            std::cerr<< "6" << std::endl;
            get_MAPQ(all_nams1, n_max1, mapq1);
            get_MAPQ(all_nams2, n_max2, mapq2);
//            std::cerr<< "7" << std::endl;
            sam.add_pair(sam_aln1, sam_aln2, record1, record2, read1_rc, read2_rc, mapq1, mapq2, mu, sigma, true);

            if ((isize_est.sample_size < 400) && ((sam_aln1.ed + sam_aln2.ed) < 3) && !sam_aln1.not_proper && !sam_aln2.not_proper ){
                int d = sam_aln1.ref_start > sam_aln2.ref_start ? sam_aln1.ref_start - sam_aln2.ref_start : sam_aln2.ref_start - sam_aln1.ref_start;
                if ( d < 2000){
//                    std::cerr<< "8 " << sample_size << std::endl;
                    float e;
                    e = d - isize_est.mu;
                    isize_est.mu = isize_est.mu + e/isize_est.sample_size; // (1.0/(sample_size +1.0)) * (sample_size*mu + d);
                    isize_est.SSE = isize_est.SSE + e*(d-isize_est.mu);
                    isize_est.V = isize_est.sample_size > 1 ? isize_est.SSE/(isize_est.sample_size -1.0) : isize_est.SSE; //d < 1000 ? ((sample_size +1.0)/sample_size) * ( (V*sample_size/(sample_size +1)) + ((mu-d)*(mu-d))/sample_size ) : V;
                    isize_est.sigma = std::sqrt( isize_est.V );
                    isize_est.sample_size = isize_est.sample_size + 1.0;
                    if (isize_est.mu < 0){
                        std::cerr<< "mu negative, mu: " << isize_est.mu << " sigma: " << isize_est.sigma << " SSE: " << isize_est.SSE << " sample size: " << isize_est.sample_size << std::endl;
                    }
                    if (isize_est.SSE < 0){
                        std::cerr<< "SSE negative, mu: " << isize_est.mu << " sigma: " << isize_est.sigma << " SSE: " << isize_est.SSE << " sample size: " << isize_est.sample_size << std::endl;
                    }
                }
            }

            return;
        }
        else{ // do full search of highest scoring pair
//            std::cerr << "Joint search" << std::endl;

            //////////////////////////// NEW ////////////////////////////////////
            /////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////
            // Get top hit counts for all locations. The joint hit count is the sum of hits of the two mates. Then align as long as score dropoff or cnt < 20

            std::vector<std::tuple<int,nam,nam>> joint_NAM_scores; // (score, aln1, aln2)
            get_best_scoring_NAM_locations(all_nams1, all_nams2, joint_NAM_scores, mu, sigma, added_n1, added_n2 );
            auto nam_max = joint_NAM_scores[0];
            auto max_score = std::get<0>(nam_max);

            robin_hood::unordered_map<int,alignment> is_aligned1;
            robin_hood::unordered_map<int,alignment> is_aligned2;
            alignment a1_indv_max;
//            a1_indv_max.sw_score = -10000;
            auto n1_max = all_nams1[0];
            get_alignment(aln_params, n1_max, ref_len_map, ref_seqs, read1, read1_rc, a1_indv_max, k, rc_already_comp1,
                          log_vars.did_not_fit, log_vars.tot_ksw_aligned);
            is_aligned1[n1_max.nam_id] = a1_indv_max;
            log_vars.tot_all_tried ++;
            alignment a2_indv_max;
//            a2_indv_max.sw_score = -10000;
            auto n2_max = all_nams2[0];
            get_alignment(aln_params, n2_max, ref_len_map, ref_seqs, read2, read2_rc, a2_indv_max, k, rc_already_comp2,
                          log_vars.did_not_fit, log_vars.tot_ksw_aligned);
            is_aligned2[n2_max.nam_id] = a2_indv_max;
            log_vars.tot_all_tried ++;

//            int a, b;
            std::string r_tmp;
//            int min_ed1, min_ed2 = 1000;
//            bool new_opt1, new_opt2 = false;
//            bool a1_is_rc, a2_is_rc;
//            int ref_start, ref_len, ref_end;
//            std::cerr << "LOOOOOOOOOOOOOOOOOOOL " << min_ed << std::endl;
            std::vector<std::tuple<double,alignment,alignment>> high_scores; // (score, aln1, aln2)
            for (auto &t : joint_NAM_scores) {
                auto score_ = std::get<0>(t);
                auto n1 = std::get<1>(t);
                auto n2 = std::get<2>(t);
                score_dropoff1 = (float) score_ / max_score;
//                std::cerr << "Min ed: " << min_ed << std::endl;
                if ( (cnt >= max_tries) || (score_dropoff1 < dropoff) ){ // only consider top 20 if there are more.
                    break;
                }

                //////// the actual testing of base pair alignment part start ////////
                //////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////
                alignment a1;
                if (n1.ref_s >= 0) {
                    if (is_aligned1.find(n1.nam_id) != is_aligned1.end() ){
//                    std::cerr << "Already aligned a1! " << std::endl;
                        a1 = is_aligned1[n1.nam_id];
                    } else {
//                    std::cerr << query_acc1 << std::endl;
                        get_alignment(aln_params, n1, ref_len_map, ref_seqs, read1, read1_rc, a1, k,
                                      rc_already_comp1,
                                      log_vars.did_not_fit, log_vars.tot_ksw_aligned);
                        is_aligned1[n1.nam_id] = a1;
                        log_vars.tot_all_tried++;
                    }
                } else { //rescue
//                    std::cerr << "RESCUE HERE1" << std::endl;
                    //////// Force SW alignment to rescue mate /////////
//                    std::cerr << query_acc2 << " RESCUE MATE 1" << a1.is_rc << " " n1.is_rc << std::endl;
                    rescue_mate(aln_params, n2, ref_len_map, ref_seqs, read2, read2_rc, rc_already_comp2, read1, read1_rc, a1, rc_already_comp1, log_vars.tot_ksw_aligned, mu, sigma, log_vars.tot_rescued, k);
//                    is_aligned1[n1.nam_id] = a1;
                    log_vars.tot_all_tried ++;
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
                        get_alignment(aln_params, n2, ref_len_map, ref_seqs, read2, read2_rc, a2, k,
                                      rc_already_comp2,
                                      log_vars.did_not_fit, log_vars.tot_ksw_aligned);
                        is_aligned2[n2.nam_id] = a2;
                        log_vars.tot_all_tried++;
                    }
                } else{
//                    std::cerr << "RESCUE HERE2" << std::endl;
                    //////// Force SW alignment to rescue mate /////////
//                    std::cerr << query_acc1 << " RESCUE MATE 2" << a1.is_rc << " " n1.is_rc << std::endl;
                    rescue_mate(aln_params, n1, ref_len_map, ref_seqs, read1, read1_rc, rc_already_comp1, read2, read2_rc, a2, rc_already_comp2, log_vars.tot_ksw_aligned, mu, sigma, log_vars.tot_rescued, k);
//                    is_aligned2[n2.nam_id] = a2;
                    log_vars.tot_all_tried ++;
                }
//                a2_indv_max = a2.sw_score >  a2_indv_max.sw_score ? a2 : a2_indv_max;
//                min_ed = a2.ed < min_ed ? a2.ed : min_ed;

                if (a2.sw_score >  a2_indv_max.sw_score){
                    a2_indv_max = a2;
//                    cnt = 0;
                }
                //////////////////////////////////////////////////////////////////

                bool r1_r2 = a2.is_rc && (a1.ref_start < a2.ref_start) && ((a2.ref_start - a1.ref_start) < mu+5*sigma); // r1 ---> <---- r2
                bool r2_r1 = a1.is_rc && (a2.ref_start < a1.ref_start) && ((a1.ref_start - a2.ref_start) < mu+5*sigma); // r2 ---> <---- r1

                if ( r1_r2 || r2_r1){
                    x = a1.ref_start > a2.ref_start ? (float) (a1.ref_start - a2.ref_start) : (float)(a2.ref_start - a1.ref_start);
                    S = (double)a1.sw_score + (double)a2.sw_score + log( normal_pdf(x, mu, sigma ) );  //* (1 - s2 / s1) * min_matches * log(s1);
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

            // Calculate joint MAPQ score
            int n_mappings = high_scores.size();
            if (n_mappings > 1) {
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
                if ( !same_pos || !same_ref){
                    get_joint_MAPQ_from_alingments(S1, S2, mapq1, mapq2);
                } else if (n_mappings > 2){ // individually highest alignment score was the same alignment as the joint highest score - calculate mapq relative to third best
                    auto third_aln_pair = high_scores[2];
                    auto S2 = std::get<0>(third_aln_pair);
//                    std::cerr << "FOR MAPQ " << S1 << " " << S2 << std::endl;
                    get_joint_MAPQ_from_alingments(S1, S2, mapq1, mapq2);

                } else { // there was no other alignment
                    mapq1 = 60;
                    mapq2 = 60;
                }
            } else{
                mapq1 = 60;
                mapq2 = 60;
            }

            if (max_secondary == 0) {
                auto best_aln_pair = high_scores[0];
                sam_aln1 = std::get<1>(best_aln_pair);
                sam_aln2 = std::get<2>(best_aln_pair);
//            get_MAPQ_aln(sam_aln1, sam_aln2);
                sam.add_pair(sam_aln1, sam_aln2, record1, record2, read1_rc, read2_rc,
                              mapq1, mapq2, mu, sigma, true);
            } else {
                int max_out = std::min(high_scores.size(), max_secondary);
//                std::cout << high_scores.size() << std::endl;
                // remove eventual duplicates - comes from, e.g., adding individual best alignments above (if identical to joint best alignment)
                auto best_aln_pair = high_scores[0];
                float s_max = std::get<0>(best_aln_pair);
                sam_aln1 = std::get<1>(best_aln_pair);
                sam_aln2 = std::get<2>(best_aln_pair);
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
//                    std::cout << i << " " << sam_aln1.ref_start << " " << sam_aln2.ref_start << " " << s_score << std::endl;
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

                    if ((s_max - s_score) < secondary_droppoff){
                        sam.add_pair(sam_aln1, sam_aln2, record1, record2, read1_rc, read2_rc,
                                      mapq1, mapq2, mu, sigma, is_primary);
                    } else{
                        break;
                    }

                    prev_start_m1 = sam_aln1.ref_start;
                    prev_start_m2 = sam_aln2.ref_start;
                    prev_ref_id_m1 = sam_aln1.ref_id;
                    prev_ref_id_m2 = sam_aln2.ref_id;
                }
            }

        }
    } else if (all_nams1.size() > 0 ) { // rescue read 2
//        std::cerr << "Rescue read 2 mode" << std::endl;
        n_max1 = all_nams1[0];
        std::vector<alignment> aln_scores1;
        std::vector<alignment> aln_scores2;
        read2_rc = reverse_complement(read2);
        for (auto &n : all_nams1) {
            score_dropoff1 = (float) n.n_hits / n_max1.n_hits;
            if ( (cnt1 >= max_tries) || score_dropoff1 < dropoff){ // only consider top 20 hits as minimap2 and break if alignment is exact match to reference or the match below droppoff cutoff.
                break;
            }
            //////// the actual testing of base pair alignment part start /////////
            alignment a1;
//            std::cerr << query_acc1 << " force rescue"  << std::endl;
            get_alignment(aln_params, n, ref_len_map, ref_seqs, read1, read1_rc, a1, k, rc_already_comp1, log_vars.did_not_fit, log_vars.tot_ksw_aligned);
            aln_scores1.push_back(a1);
            //////////////////////////////////////////////////////////////////

            //////// Force SW alignment to rescue mate /////////
            alignment a2;
//            std::cerr << query_acc2 << " force rescue" << std::endl;
            rescue_mate(aln_params, n, ref_len_map, ref_seqs, read1, read1_rc, rc_already_comp1, read2, read2_rc, a2, rc_already_comp2, log_vars.tot_ksw_aligned, mu, sigma, log_vars.tot_rescued, k);
            aln_scores2.push_back(a2);
            //////////////////////////////////////////////////////////////////

            cnt1 ++;
            log_vars.tot_all_tried ++;
        }
        std::sort(aln_scores1.begin(), aln_scores1.end(), score_sw);
        std::sort(aln_scores2.begin(), aln_scores2.end(), score_sw);

        // Calculate best combined score here
        std::vector<std::tuple<double,alignment,alignment>> high_scores; // (score, aln1, aln2)
        get_best_scoring_pair(aln_scores1, aln_scores2, high_scores, mu, sigma );

        // Calculate joint MAPQ score
        if (high_scores.size() > 1) {
            auto best_aln_pair = high_scores[0];
            auto S1 = std::get<0>(best_aln_pair);
            auto second_aln_pair = high_scores[1];
            auto S2 = std::get<0>(second_aln_pair);
            get_joint_MAPQ_from_alingments(S1, S2, mapq1, mapq2);

        } else{
            mapq1 = 60;
            mapq2 = 60;
        }

        // append both alignments to string here
        if (max_secondary == 0){
            auto best_aln_pair = high_scores[0];
            sam_aln1 = std::get<1>(best_aln_pair);
            sam_aln2 = std::get<2>(best_aln_pair);
//            get_MAPQ(all_nams1, n_max1, mapq1);
//            mapq2 = 0;
            sam.add_pair(sam_aln1, sam_aln2, record1, record2, read1_rc, read2_rc, mapq1, mapq2, mu, sigma, true);
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
                sam_aln1 = std::get<1>(aln_pair);
                sam_aln2 = std::get<2>(aln_pair);
                if ((s_max - s_score) < secondary_droppoff){
                    sam.add_pair(sam_aln1, sam_aln2, record1, record2, read1_rc, read2_rc,
                                  mapq1, mapq2, mu, sigma, is_primary);
                } else{
                    break;
                }
            }
        }

    } else if (all_nams2.size() > 0 ) { // rescue read 1
//        std::cerr << "Rescue read 1 mode" << std::endl;
        n_max2 = all_nams2[0];
        std::vector<alignment> aln_scores1;
        std::vector<alignment> aln_scores2;
        read1_rc = reverse_complement(read1);
        for (auto &n : all_nams2) {
            score_dropoff2 = (float) n.n_hits / n_max2.n_hits;
            if ( (cnt2 >= max_tries) || score_dropoff2 < dropoff){ // only consider top 20 hits as minimap2 and break if alignment is exact match to reference or the match below droppoff cutoff.
                break;
            }
//            std::cerr << "here1" << std::endl;

            //////// the actual testing of base pair alignment part start /////////
            alignment a2;
            get_alignment(aln_params, n, ref_len_map, ref_seqs, read2, read2_rc, a2, k, rc_already_comp2, log_vars.did_not_fit, log_vars.tot_ksw_aligned);
            aln_scores2.push_back(a2);
            //////////////////////////////////////////////////////////////////

//            std::cerr << "here2" << std::endl;

            //////// Force SW alignment to rescue mate /////////
            alignment a1;
            rescue_mate(aln_params, n, ref_len_map, ref_seqs, read2, read2_rc, rc_already_comp2, read1, read1_rc, a1, rc_already_comp1, log_vars.tot_ksw_aligned, mu, sigma, log_vars.tot_rescued, k);
            aln_scores1.push_back(a1);
            //////////////////////////////////////////////////////////////////

            cnt2 ++;
            log_vars.tot_all_tried ++;
        }
        std::sort(aln_scores1.begin(), aln_scores1.end(), score_sw);
        std::sort(aln_scores2.begin(), aln_scores2.end(), score_sw);

        // Calculate best combined score here
        std::vector<std::tuple<double,alignment,alignment>> high_scores; // (score, aln1, aln2)
        get_best_scoring_pair(aln_scores1, aln_scores2, high_scores, mu, sigma );

        // Calculate joint MAPQ score
        if (high_scores.size() > 1) {
            auto best_aln_pair = high_scores[0];
            auto S1 = std::get<0>(best_aln_pair);
            auto second_aln_pair = high_scores[1];
            auto S2 = std::get<0>(second_aln_pair);
            get_joint_MAPQ_from_alingments(S1, S2, mapq1, mapq2);

        } else{
            mapq1 = 60;
            mapq2 = 60;
        }

        // append both alignments to string here
        if(max_secondary == 0) {
            auto best_aln_pair = high_scores[0];
            sam_aln1 = std::get<1>(best_aln_pair);
            sam_aln2 = std::get<2>(best_aln_pair);
//            get_MAPQ(all_nams2, n_max2, mapq2);
//            mapq1 = 0;
            sam.add_pair(sam_aln1, sam_aln2, record1, record2, read1_rc, read2_rc,
                          mapq1, mapq2, mu, sigma, true);
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
                sam_aln1 = std::get<1>(aln_pair);
                sam_aln2 = std::get<2>(aln_pair);
                if ((s_max - s_score) < secondary_droppoff){
                    sam.add_pair(sam_aln1, sam_aln2, record1, record2, read1_rc, read2_rc,
                                  mapq1, mapq2, mu, sigma, is_primary);
                } else{
                    break;
                }
            }
        }

    }

}


inline void get_best_map_location(std::vector<std::tuple<int,nam,nam>> joint_NAM_scores, std::vector<nam> &nams1, std::vector<nam> &nams2, i_dist_est &isize_est, nam &best_nam1,  nam &best_nam2 ) {
    robin_hood::unordered_set<int> added_n1;
    robin_hood::unordered_set<int> added_n2;
    get_best_scoring_NAM_locations(nams1, nams2, joint_NAM_scores, isize_est.mu, isize_est.sigma, added_n1, added_n2 );
    nam n1_joint_max, n2_joint_max, n1_indiv_max, n2_indiv_max;
    float score_joint = 0;
    float score_indiv = 0;
    best_nam1.ref_s = -1; //Unmapped until proven mapped
    best_nam2.ref_s = -1; //Unmapped until proven mapped
    if (joint_NAM_scores.size() > 0) {
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
        if (nams1.size() > 0) {
            auto n1_indiv_max = nams1[0];
            score_indiv += n1_indiv_max.score - (n1_indiv_max.score/2.0); //Penalty for being mapped individually
            best_nam1 = n1_indiv_max;
        }
        if (nams2.size() > 0) {
            auto n2_indiv_max = nams2[0];
            score_indiv += n2_indiv_max.score - (n2_indiv_max.score/2.0); //Penalty for being mapped individually
            best_nam2 = n2_indiv_max;
        }
        if ( score_joint > score_indiv ){ // joint score is better than individual
//            std::cerr << "HERE " << score_joint << " " << score_indiv << std::endl;
            best_nam1 = n1_joint_max;
            best_nam2 = n2_joint_max;
        }

        if ((isize_est.sample_size < 400) && (score_joint > score_indiv) ){
            int d = n1_joint_max.ref_s > n2_joint_max.ref_s ? n1_joint_max.ref_s - n2_joint_max.ref_s : n2_joint_max.ref_s - n1_joint_max.ref_s;
//            std::cerr << "HERE " << d << " " << mu <<  " " << sigma << " "<< n1_joint_max.ref_s << " " <<  n2_joint_max.ref_s << " "<< n1_joint_max.score << " " <<  n2_joint_max.score << std::endl;
            if ( d < 2000){
                float e;
                e = d - isize_est.mu;
                isize_est.mu = isize_est.mu + e/isize_est.sample_size; // (1.0/(sample_size +1.0)) * (sample_size*mu + d);
                isize_est.SSE = isize_est.SSE + e*(d-isize_est.mu);
                isize_est.V = isize_est.sample_size > 1 ? isize_est.SSE/(isize_est.sample_size -1.0) : isize_est.SSE; //d < 1000 ? ((sample_size +1.0)/sample_size) * ( (V*sample_size/(sample_size +1)) + ((mu-d)*(mu-d))/sample_size ) : V;
                isize_est.sigma = std::sqrt( isize_est.V );
                isize_est.sample_size = isize_est.sample_size + 1.0;
            }
        }
    }

}


void align_PE_read(KSeq &record1, KSeq &record2, std::string &outstring, logging_variables &log_vars, i_dist_est &isize_est, alignment_params &aln_params,
        mapping_params &map_param, std::vector<unsigned int> &ref_lengths, std::vector<std::string> &ref_seqs, kmer_lookup &mers_index, mers_vector &flat_vector, ref_names &acc_map ){
    // Declare variables
    mers_vector_read query_mers1, query_mers2; // pos, chr_id, kmer hash value

    // generate mers here
    auto strobe_start = std::chrono::high_resolution_clock::now();
//                std::cerr << "Going in! " << std::endl;
    query_mers1 = seq_to_randstrobes2_read(map_param.n, map_param.k, map_param.w_min, map_param.w_max, record1.seq, 0, map_param.s, map_param.t_syncmer,
                                           map_param.q,
                                           map_param.max_dist);
//                std::cerr << "Lolz1 " << std::endl;
    query_mers2 = seq_to_randstrobes2_read(map_param.n, map_param.k, map_param.w_min, map_param.w_max, record2.seq, 0, map_param.s, map_param.t_syncmer,
                                           map_param.q,
                                           map_param.max_dist);
//                std::cerr << "Lolz2 " << std::endl;
    auto strobe_finish = std::chrono::high_resolution_clock::now();
    log_vars.tot_construct_strobemers += strobe_finish - strobe_start;
//                std::cerr << record1.name << " " << query_mers1.size() << std::endl;
//                std::cerr << record2.name << " " << query_mers2.size() << std::endl;

    // Find NAMs
    auto nam_start = std::chrono::high_resolution_clock::now();
    robin_hood::unordered_map< unsigned int, std::vector<hit>> hits_per_ref;
    hits_per_ref.reserve(100);
    std::vector<nam> nams1;
    std::vector<nam> nams2;
    std::pair<float, int> info1, info2;
    info1 = find_nams(nams1, hits_per_ref, query_mers1, flat_vector, mers_index, map_param.k, ref_seqs,
                      record1.seq, map_param.filter_cutoff);
    hits_per_ref.clear();
    info2 = find_nams(nams2, hits_per_ref, query_mers2, flat_vector, mers_index, map_param.k, ref_seqs,
                      record2.seq, map_param.filter_cutoff);
    hits_per_ref.clear();
    auto nam_finish = std::chrono::high_resolution_clock::now();
    log_vars.tot_find_nams += nam_finish - nam_start;


    if (map_param.R > 1) {
        auto rescue_start = std::chrono::high_resolution_clock::now();
        std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> hits_fw;
        std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> hits_rc;
        if ((nams1.size() == 0) || (info1.first < 0.7)) {
            hits_fw.reserve(5000);
            hits_rc.reserve(5000);
            log_vars.tried_rescue += 1;
            nams1.clear();
//                        std::cerr << "Rescue is_sam_out read 1: " << record1.name << info1.first <<  std::endl;
            find_nams_rescue(hits_fw, hits_rc, nams1, hits_per_ref, query_mers1,
                                     flat_vector,
                                     mers_index, map_param.k, ref_seqs,
                                     record1.seq, map_param.rescue_cutoff);
            hits_per_ref.clear();
            hits_fw.clear();
            hits_rc.clear();
//                    std::cerr << "Found: " << nams.size() <<  std::endl;
        }

        if ((nams2.size() == 0) || (info2.first < 0.7)) {
            hits_fw.reserve(5000);
            hits_rc.reserve(5000);
            log_vars.tried_rescue += 1;
            nams2.clear();
//                        std::cerr << "Rescue is_sam_out read 2: " << record2.name << info2.first <<  std::endl;
            find_nams_rescue(hits_fw, hits_rc, nams2, hits_per_ref, query_mers2,
                                     flat_vector,
                                     mers_index, map_param.k, ref_seqs,
                                     record2.seq, map_param.rescue_cutoff);
            hits_per_ref.clear();
            hits_fw.clear();
            hits_rc.clear();
//                    std::cerr << "Found: " << nams.size() <<  std::endl;
        }
        auto rescue_finish = std::chrono::high_resolution_clock::now();
        log_vars.tot_time_rescue += rescue_finish - rescue_start;
    }



    //Sort hits based on start choordinate on query sequence
    auto nam_sort_start = std::chrono::high_resolution_clock::now();
    std::sort(nams1.begin(), nams1.end(), score);
    std::sort(nams2.begin(), nams2.end(), score);
    auto nam_sort_finish = std::chrono::high_resolution_clock::now();
    log_vars.tot_sort_nams += nam_sort_finish - nam_sort_start;

//                std::cerr << record1.name << std::endl;
//                for (auto &n : nams1){
//                    std::cerr << "NAM ORG: " << n.ref_id << ": (" << n.score << ", " << n.n_hits << ", " << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e  << ")" << std::endl;
//                }
//                std::cerr << record2.name << std::endl;
//                for (auto &n : nams2){
//                    std::cerr << "NAM ORG: " << n.ref_id << ": (" << n.score << ", " << n.n_hits << ", " << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e  << ")" << std::endl;
//                }

    auto extend_start = std::chrono::high_resolution_clock::now();
    if (!map_param.is_sam_out) {
//                    output_hits_paf(output_streams[omp_get_thread_num()], nams1, record1.name, acc_map, k,
//                                    record1.seq.length(), ref_lengths);
        nam nam_read1;
        nam nam_read2;
        std::vector<std::tuple<int, nam, nam>> joint_NAM_scores;
        get_best_map_location(joint_NAM_scores, nams1, nams2, isize_est,
                              nam_read1,
                              nam_read2);
        output_hits_paf_PE(outstring, nam_read1, record1.name,
                           acc_map,
                           map_param.k,
                           record1.seq.length(), ref_lengths);
        output_hits_paf_PE(outstring, nam_read2, record2.name,
                           acc_map,
                           map_param.k,
                           record2.seq.length(), ref_lengths);
        joint_NAM_scores.clear();
    } else {
        Sam sam(outstring, acc_map);
        align_PE(aln_params, sam, nams1, nams2, record1,
                 record2,
                 map_param.k,
                 ref_lengths, ref_seqs, log_vars,
                 map_param.dropoff_threshold, isize_est, map_param.maxTries, map_param.max_secondary);
    }
    auto extend_finish = std::chrono::high_resolution_clock::now();
    log_vars.tot_extend += extend_finish - extend_start;
    nams1.clear();
    nams2.clear();

}



void align_SE_read(KSeq &record, std::string &outstring, logging_variables &log_vars, alignment_params &aln_params,
                   mapping_params &map_param, std::vector<unsigned int> &ref_lengths, std::vector<std::string> &ref_seqs, kmer_lookup &mers_index, mers_vector &flat_vector, ref_names &acc_map ){


        std::string seq, seq_rc;
        unsigned int q_id = 0;
        std::pair<float, int> info;
        mers_vector_read query_mers; // pos, chr_id, kmer hash value
        std::vector<nam> nams; // (r_id, r_pos_start, r_pos_end, q_pos_start, q_pos_end)
        robin_hood::unordered_map< unsigned int, std::vector<hit>> hits_per_ref;
        std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> hits_fw;
        std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> hits_rc;
        hits_per_ref.reserve(100);
        hits_fw.reserve(5000);
        hits_rc.reserve(5000);


        // generate mers here
        auto strobe_start = std::chrono::high_resolution_clock::now();
        query_mers = seq_to_randstrobes2_read(map_param.n, map_param.k, map_param.w_min, map_param.w_max, record.seq, q_id, map_param.s, map_param.t_syncmer, map_param.q, map_param.max_dist);
        auto strobe_finish = std::chrono::high_resolution_clock::now();
        log_vars.tot_construct_strobemers += strobe_finish - strobe_start;

        // Find NAMs
        auto nam_start = std::chrono::high_resolution_clock::now();
        info = find_nams(nams, hits_per_ref, query_mers, flat_vector, mers_index, map_param.k, ref_seqs, record.seq, map_param.filter_cutoff);
        hits_per_ref.clear();
        auto nam_finish = std::chrono::high_resolution_clock::now();
        log_vars.tot_find_nams += nam_finish - nam_start;

        if (map_param.R > 1) {
            auto rescue_start = std::chrono::high_resolution_clock::now();
            if ((nams.size() == 0) || (info.first < 0.7)) {
                log_vars.tried_rescue += 1;
                nams.clear();
                find_nams_rescue(hits_fw, hits_rc, nams, hits_per_ref, query_mers, flat_vector, mers_index, map_param.k, ref_seqs,
                                        record.seq, map_param.rescue_cutoff);
                hits_per_ref.clear();
                hits_fw.clear();
                hits_rc.clear();
            }
            auto rescue_finish = std::chrono::high_resolution_clock::now();
            log_vars.tot_time_rescue += rescue_finish - rescue_start;
        }

        //Sort hits on score
        auto nam_sort_start = std::chrono::high_resolution_clock::now();
        std::sort(nams.begin(), nams.end(), score);
        auto nam_sort_finish = std::chrono::high_resolution_clock::now();
        log_vars.tot_sort_nams += nam_sort_finish - nam_sort_start;

        auto extend_start = std::chrono::high_resolution_clock::now();
        if (!map_param.is_sam_out) {
            output_hits_paf(outstring, nams, record.name, acc_map, map_param.k,
                            record.seq.length(), ref_lengths);
        } else {
            Sam sam(outstring, acc_map);
            if (map_param.max_secondary > 0){
                // I created an entire new function here, duplicating a lot of the code as outputting secondary hits is has some overhead to the
                // original align_SE function (storing a vector of hits and sorting them)
                // Such overhead is not present in align_PE - which implements both options in the same function.

                align_SE_secondary_hits(aln_params, sam, nams, record, map_param.k,
                         ref_lengths, ref_seqs, log_vars, map_param.dropoff_threshold, map_param.maxTries, map_param.max_secondary);
            } else {
                align_SE(aln_params, sam, nams, record, map_param.k,
                         ref_lengths, ref_seqs, log_vars,  map_param.dropoff_threshold, map_param.maxTries);
            }
        }
        auto extend_finish = std::chrono::high_resolution_clock::now();
        log_vars.tot_extend += extend_finish - extend_start;
        q_id++;
        nams.clear();
}


