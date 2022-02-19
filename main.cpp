#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <assert.h>
#include <math.h>
#include <chrono>  // for high_resolution_clock
#include <omp.h>
#include <zlib.h>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <inttypes.h>
#include "source/kseq++.hpp"
using namespace klibpp;
#include "source/robin_hood.h"
#include "source/index.hpp"
//#include "gap_affine/affine_wavefront_align.h"
#include "source/ksw2.h"
#include "source/ssw_cpp.h"
//#include "source/parasail/parasail.h"


//develop
#include <chrono>
#include <thread>
#include <sstream>


typedef robin_hood::unordered_map<int, std::string > idx_to_acc;
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))


typedef struct { int match;
                int mismatch;
                int gap_open;
                int gap_extend;} alignment_params;


static uint64_t read_references(std::vector<std::string> &seqs, std::vector<unsigned int> &lengths, idx_to_acc &acc_map, std::string fn)
{
    uint64_t total_ref_seq_size = 0;
    std::ifstream file(fn);
    std::string line, seq;
    int ref_index = 0;
    while (getline(file, line)) {
        if (line[0] == '>') {
//            std::cerr << ref_index << " " << line << std::endl;
            if (seq.length() > 0){
//                seqs[ref_index -1] = seq;
                seqs.push_back(seq);
                lengths.push_back(seq.length());
                total_ref_seq_size += seq.length();
//                std::cerr << ref_index - 1 << " here " << seq << " " << seq.length() << " " << seq.size() << std::endl;
//                generate_kmers(h, k, seq, ref_index);
            }
//            acc_map[ref_index] = line.substr(1, line.length() -1); //line;
            acc_map[ref_index] = line.substr(1, line.find(' ') -1); // cutting at first space;
            ref_index++;
            seq = "";
        }
        else {
            seq += line;
        }
    }
    if (seq.length() > 0){
//        seqs[ref_index -1] = seq;
        seqs.push_back(seq);
        lengths.push_back(seq.length());
        total_ref_seq_size += seq.length();
//        std::cerr << ref_index -1 << " here2 " << seq << std::endl;
//        generate_kmers(h, k, seq, ref_index);
    }
    file.close();
    return total_ref_seq_size;
}



static inline void print_diagnostics(mers_vector &ref_mers, kmer_lookup &mers_index, std::string logfile_name, int k) {
    // Prins to csv file the statistics on the number of seeds of a particular length and what fraction of them them are unique in the index:
    // format:
    // seed_length, count, percentage_unique
    //

    int max_size = 100000;
    std::vector<int> log_count(max_size,0);  // stores count and each index represents the length
    std::vector<int> log_unique(max_size,0); // stores count unique and each index represents the length
    std::vector<int> log_repetitive(max_size,0); // stores count unique and each index represents the length


    int seed_length;
    for (auto &it : mers_index) {
        auto hash_refmer = it.first;
        auto ref_mer = it.second;

        auto offset = std::get<0>(ref_mer);
        auto count = std::get<1>(ref_mer);


        for (size_t j = offset; j < offset + count; ++j) {
            auto r = ref_mers[j];
            auto p = std::get<2>(r);
            int bit_alloc = 8;
            int r_id = (p >> bit_alloc);
            int mask=(1<<bit_alloc) - 1;
            int offset = (p & mask);
            seed_length =  offset + k;
            if (seed_length < max_size){

                log_count[seed_length] ++;
            } else {
               std::cerr << "Detected seed size over " << max_size << " bp (can happen, e.g., over centromere): " << seed_length << std::endl;
            }

        }

        if ( (count == 1) & (seed_length < max_size) ) {
            log_unique[seed_length] ++;
        }
        if ( (count >= 10) & (seed_length < max_size) ) {
            log_repetitive[seed_length] ++;
        }
    }

//    std::cerr << "Here" << std::endl;

    // printing
    std::ofstream log_file;
    log_file.open(logfile_name);
    for (int i=0 ; i < log_count.size(); ++i) {
        if (log_count[i] > 0) {
            log_file << i << ',' << log_count[i] << ',' << (float) log_unique[i] / (float) log_count[i] << ',' << (float) log_repetitive[i] / (float) log_count[i] << std::endl;
        }
    }
    log_file.close();

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

static inline std::pair<float,int> find_nams_rescue(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> hits_fw, std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> hits_rc, std::vector<nam> &final_nams, robin_hood::unordered_map< unsigned int, std::vector<hit>> &hits_per_ref, mers_vector_read &query_mers, mers_vector &ref_mers, kmer_lookup &mers_index, int k, std::vector<std::string> &ref_seqs, std::string &read, unsigned int filter_cutoff ){
    std::pair<float,int> info (0,0); // (nr_nonrepetitive_hits/total_hits, max_nam_n_hits)
    int nr_good_hits = 0, total_hits = 0;
    bool is_rc = true, no_rep_fw = true, no_rep_rc = true;
    std::pair<int, int> repeat_fw(0,0), repeat_rc(0,0);
    std::vector<std::pair<int, int>> repetitive_fw, repetitive_rc;
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
                if (count > filter_cutoff){
                    if (no_rep_rc){ //initialize
                        repeat_rc.first = query_s;
                        repeat_rc.second = query_e;
                        no_rep_rc = false;
                    }
                    else if (query_s >= repeat_rc.second){
                        repetitive_rc.push_back(repeat_rc);
                        repeat_rc.first = query_s;
                        repeat_rc.second = query_e;
                    } else{
                        repeat_rc.second = repeat_rc.second < query_e ? query_e : repeat_rc.second;
                    }
                } else{
                    nr_good_hits ++;
                }
            } else{
                std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool> s(count, offset, query_s, query_e, is_rc);
                hits_fw.push_back(s);
                if (count > filter_cutoff){
                    if (no_rep_fw){ //initialize
                        repeat_fw.first = query_s;
                        repeat_fw.second = query_e;
                        no_rep_fw = false;
                    }
                    else if (query_s >= repeat_fw.second ){
                        repetitive_fw.push_back(repeat_fw);
                        repeat_fw.first = query_s;
                        repeat_fw.second = query_e;
                    } else{
                        repeat_fw.second = repeat_fw.second < query_e ? query_e : repeat_fw.second;
                    }
                } else{
                    nr_good_hits ++;
                }
            }


        }
    }
    if (!no_rep_fw) {
        repetitive_fw.push_back(repeat_fw);
    }
    if (!no_rep_rc) {
        repetitive_rc.push_back(repeat_rc);
    }
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
//                ref_d = std::get<3>(r) + k - std::get<2>(r);
//                int diff = (h.query_e - h.query_s) - ref_d > 0 ? (h.query_e - h.query_s) -  ref_d : ref_d - (h.query_e - h.query_s);
//                if (diff <= min_diff ){
//                    min_diff = diff;
//                }
//            }

            for(size_t j = offset; j < offset+count; ++j)
            {
                auto r = ref_mers[j];
                h.ref_s = std::get<1>(r);
                auto p = std::get<2>(r);
                int bit_alloc = 8;
                int r_id = (p >> bit_alloc);
                int mask=(1<<bit_alloc) - 1;
                int offset = (p & mask);
                h.ref_e = h.ref_s + offset + k;
//                h.count = count;
//                hits_per_ref[std::get<1>(r)].push_back(h);

                int diff = (h.query_e - h.query_s) - (h.ref_e - h.ref_s) > 0 ? (h.query_e - h.query_s) - (h.ref_e - h.ref_s) : (h.ref_e - h.ref_s) - (h.query_e - h.query_s);
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
//                ref_d = std::get<3>(r) + k - std::get<2>(r);
//                int diff = (h.query_e - h.query_s) - ref_d > 0 ? (h.query_e - h.query_s) -  ref_d : ref_d - (h.query_e - h.query_s);
//                if (diff <= min_diff ){
//                    min_diff = diff;
//                }
//            }

            for(size_t j = offset; j < offset+count; ++j)
            {
                auto r = ref_mers[j];
                h.ref_s = std::get<1>(r);
                auto p = std::get<2>(r);
                int bit_alloc = 8;
                int r_id = (p >> bit_alloc);
                int mask=(1<<bit_alloc) - 1;
                int offset = (p & mask);
                h.ref_e = h.ref_s + offset + k;
//                h.count = count;
//                hits_per_ref[std::get<1>(r)].push_back(h);
                int diff = (h.query_e - h.query_s) - (h.ref_e - h.ref_s) > 0 ? (h.query_e - h.query_s) - (h.ref_e - h.ref_s) : (h.ref_e - h.ref_s) - (h.query_e - h.query_s);
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
    info.first = total_hits > 0 ? ((float) nr_good_hits) / ((float) total_hits) : 1.0;
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
            if (not is_added){
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
                        int n_max_span = MAX((n.query_e - n.query_s), (n.ref_e - n.ref_s));
                        int n_min_span = MIN((n.query_e - n.query_s), (n.ref_e - n.ref_s));
                        float n_score;
                        n_score = ( 2*n_min_span -  n_max_span) > 0 ? (float) (n.n_hits * ( 2*n_min_span -  n_max_span) ) : 1;   // this is really just n_hits * ( min_span - (offset_in_span) ) );
//                        n_score = n.n_hits * (n.query_e - n.query_s);
                        n.score = n_score;
                        final_nams.push_back(n);
                        max_nam_n_hits = n.n_hits > max_nam_n_hits ? n.n_hits : max_nam_n_hits;
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
            int n_max_span = MAX((n.query_e - n.query_s), (n.ref_e - n.ref_s));
            int n_min_span = MIN((n.query_e - n.query_s), (n.ref_e - n.ref_s));
            float n_score;
            n_score = ( 2*n_min_span -  n_max_span) > 0 ? (float) (n.n_hits * ( 2*n_min_span -  n_max_span) ) : 1;   // this is really just n_hits * ( min_span - (offset_in_span) ) );
//            n_score = n.n_hits * (n.query_e - n.query_s);
            n.score = n_score;
            final_nams.push_back(n);
            max_nam_n_hits = n.n_hits > max_nam_n_hits ? n.n_hits : max_nam_n_hits;
        }
    }

//    for (auto &n : final_nams){
//        std::cerr << "RESCUE NAM: " << n.ref_id << ": (" << n.score << ", " << n.n_hits << ", " << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e  << ")" << " " <<  n.is_rc << std::endl;
//    }
    info.second = max_nam_n_hits;
    return info;

}



static inline std::pair<float,int> find_nams(std::vector<nam> &final_nams, robin_hood::unordered_map< unsigned int, std::vector<hit>> &hits_per_ref, mers_vector_read &query_mers, mers_vector &ref_mers, kmer_lookup &mers_index, int k, std::vector<std::string> &ref_seqs, std::string &read, unsigned int filter_cutoff ){
//    std::cerr << "ENTER FIND NAMS " <<  std::endl;
//    robin_hood::unordered_map< unsigned int, std::vector<hit>> hits_per_ref; // [ref_id] -> vector( struct hit)
//    std::vector<std::vector<hit>> hits_per_ref(10);
//    int read_length = read.length();
//    std::cerr << " "  <<  std::endl;
    std::pair<float,int> info (0,0); // (nr_nonrepetitive_hits/total_hits, max_nam_n_hits)
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
//                unsigned int ref_id = std::get<0>(r);
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
//                    ref_d = std::get<3>(r) + k - std::get<2>(r);
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
//                    unsigned int ref_id = std::get<0>(r);
//                    unsigned int ref_s = std::get<1>(r);
//                    unsigned int ref_e = std::get<2>(r) + k; //ref_s + read_length/2;
                    h.ref_s = std::get<1>(r);
                    auto p = std::get<2>(r);
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
                    int diff = (h.query_e - h.query_s) - (h.ref_e - h.ref_s) > 0 ? (h.query_e - h.query_s) - (h.ref_e - h.ref_s) : (h.ref_e - h.ref_s) - (h.query_e - h.query_s);
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
            if (not is_added){
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
                        int n_max_span = MAX((n.query_e - n.query_s), (n.ref_e - n.ref_s));
                        int n_min_span = MIN((n.query_e - n.query_s), (n.ref_e - n.ref_s));
                        float n_score;
                        n_score = ( 2*n_min_span -  n_max_span) > 0 ? (float) (n.n_hits * ( 2*n_min_span -  n_max_span) ) : 1;   // this is really just n_hits * ( min_span - (offset_in_span) ) );
//                        n_score = n.n_hits * (n.query_e - n.query_s);
                        n.score = n_score;
                        final_nams.push_back(n);
                        max_nam_n_hits = n.n_hits > max_nam_n_hits ? n.n_hits : max_nam_n_hits;
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
            int n_max_span = MAX((n.query_e - n.query_s), (n.ref_e - n.ref_s));
            int n_min_span = MIN((n.query_e - n.query_s), (n.ref_e - n.ref_s));
            float n_score;
            n_score = ( 2*n_min_span -  n_max_span) > 0 ? (float) (n.n_hits * ( 2*n_min_span -  n_max_span) ) : 1;   // this is really just n_hits * ( min_span - (offset_in_span) ) );
//            n_score = n.n_hits * (n.query_e - n.query_s);
            n.score = n_score;
            final_nams.push_back(n);
            max_nam_n_hits = n.n_hits > max_nam_n_hits ? n.n_hits : max_nam_n_hits;
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

static inline bool score(const nam &a, const nam &b)
{
    return ( a.score > b.score );
}


static inline void output_hits_paf(std::string &paf_output, std::vector<nam> &all_nams, std::string query_acc, idx_to_acc &acc_map, int k, int read_len, std::vector<unsigned int> &ref_len_map) {
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

static inline void output_hits_paf_PE(std::string &paf_output, nam &n, std::string &query_acc, idx_to_acc &acc_map, int k, int read_len, std::vector<unsigned int> &ref_len_map) {
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

static inline std::string reverse_complement(std::string &read) {
    auto read_rev = read;
    std::reverse(read_rev.begin(), read_rev.end()); // reverse
//    std::cerr << read_rev << std::endl;
    for (size_t j = 0; j < read_rev.length(); ++j) { // complement
        if (read_rev[j] == 'A') read_rev[j] = 'T';
        else if (read_rev[j] == 'T') read_rev[j] = 'A';
        else if (read_rev[j] == 'C') read_rev[j] = 'G';
        else if (read_rev[j] == 'G') read_rev[j] = 'C';
    }
//    std::cerr << read_rev << std::endl;
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


inline aln_info ssw_align(std::string &ref, std::string &query, int read_len, int match_score, int mismatch_penalty, int gap_opening_penalty, int gap_extending_penalty) {

    aln_info aln;
    int32_t maskLen = strlen(query.c_str())/2;
    maskLen = maskLen < 15 ? 15 : maskLen;
    if (ref.length() > 2000){
        std::cerr << "ALIGNMENT TO REF LONGER THAN 2000bp - REPORT TO DEVELOPER. Happened for read: " <<  query << " ref len:" << ref.length() << std::endl;
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
    // Aligns the query to the ref
//    bool passed;
//    std::cerr << "I'm here!" << std::endl;
//    std::cerr << "read: " << query << std::endl;
//    std::cerr << "ref: "  << ref << std::endl;
//    passed =
      aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment_ssw, maskLen);
//    std::cerr << passed << std::endl;
//    if(!passed){
//        std::cerr << "Failed" << std::endl;
//        std::cerr << "read: " << query << std::endl;
//        std::cerr << "ref: "  << ref << std::endl;
//    }


//    cout << "===== SSW result =====" << endl;
//    cout << "Best Smith-Waterman score:\t" << alignment.sw_score << endl
//         << "Next-best Smith-Waterman score:\t" << alignment.sw_score_next_best << endl
//         << "Reference start:\t" << alignment.ref_begin << endl
//         << "Reference end:\t" << alignment.ref_end << endl
//         << "Query start:\t" << alignment.query_begin << endl
//         << "Query end:\t" << alignment.query_end << endl
//         << "Next-best reference end:\t" << alignment.ref_end_next_best << endl
//         << "Number of mismatches:\t" << alignment.mismatches << endl
//         << "Cigar: " << alignment.cigar_string << endl;

    aln.global_ed = alignment_ssw.global_ed;
    aln.ed = alignment_ssw.mismatches;
    aln.ref_offset = alignment_ssw.ref_begin;
    aln.cigar = alignment_ssw.cigar_string;
    aln.sw_score = alignment_ssw.sw_score;
    return aln;
}


inline aln_info ksw_align(const char *tseq, int tlen, const char *qseq, int qlen,
               int sc_mch, int sc_mis, int gapo, int gape, ksw_extz_t &ez) {
    int8_t a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
    int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a, b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
    const uint8_t *ts = reinterpret_cast<const uint8_t *>(tseq);
    const uint8_t *qs = reinterpret_cast<const uint8_t *>(qseq);
    memset(&ez, 0, sizeof(ksw_extz_t));
    ksw_extz2_sse(0, qlen, qs, tlen, ts, 5, mat, gapo, gape, -1, -1, 10000, KSW_EZ_EXTZ_ONLY, &ez);

    aln_info aln;
//    std::string cigar_mod;
//    cigar_mod.reserve(5*ez.n_cigar);
    unsigned int tstart_offset = 0;
    int eqx_len, switch_ind;
    std::stringstream cigar_string;
    int edit_distance = 0;
    int sw_score = 0;
    unsigned ref_pos = 0, read_pos = 0;
    for (int i = 0; i < ez.n_cigar; i++) {
        int count = ez.cigar[i] >> 4;
        char op = "MID"[ez.cigar[i] & 0xf];
//        std::cerr << "count: " << count << " op:" << op << std::endl;
        if ( (i==0) && op == 'D'){
            ref_pos += count;
            tstart_offset = ref_pos;
//            std::cerr << "First deletion " << i << " " << count << std::endl;
            continue;
        }
        if ( (i==ez.n_cigar-1) && op == 'D'){
            ref_pos += count;
//            std::cerr << "Last deletion " << i << " " << count << std::endl;
            continue;
        }
        cigar_string << count << op;
        switch (op) {
            case 'M': {
//                eqx_len = 0;
//                switch_ind = 0; // switch_ind 0 if prev was match, 1 if mismatch
//                char o = '=';
                for (int j = 0; j < count; j++, ref_pos++, read_pos++) {
                    if (tseq[ref_pos] != qseq[read_pos]) {
                        edit_distance++;
                        sw_score -= -b;
//                        if ((switch_ind == 0) && (j > 0)) { // prev was match
//                            cigar_string << eqx_len << '=';
//                            eqx_len = 0;
//                        }
//                        switch_ind = 1;
//                        o = 'X';
//                        eqx_len++;
                    } else{
                        sw_score += sc_mch;
//                        if (switch_ind == 1) { // prev was mismatch
//                            cigar_string << eqx_len << 'X';
//                            eqx_len = 0;
//                            o = '=';
//                            switch_ind = 0;
//                        }
//                        eqx_len++;
                    }
                }
//                cigar_string << eqx_len << o;
                break;
            }
            case 'D': {
                edit_distance += count;
                ref_pos += count;
                sw_score -= (gapo + (count - 1));
//                cigar_string << count << op;
                break;
            }
            case 'I': {
                edit_distance += count;
                read_pos += count;
                sw_score -= (gapo + (count - 1));
//                cigar_string << count << op;
                break;
            }
            default:assert(0);
        }
//        std::cerr << "ED " << edit_distance << std::endl;
    }
    aln.ed = edit_distance;
    aln.sw_score = sw_score;
    aln.ref_offset = tstart_offset;
    aln.cigar = cigar_string.str();
    free(ez.cigar); //free(ts); free(qs);
    return aln;
}

inline int HammingDistance(std::string &One, std::string &Two)
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


inline int HammingToCigarEQX2(std::string &One, std::string &Two, std::stringstream &cigar, int match, int mismatch, int &aln_score, int &soft_left, int &soft_right)
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
    bool prev_is_match = One[start_softclipp+1] == Two[start_softclipp+1] ? true : false;
    bool beginning = true;
    int hamming_mod = prev_is_match ? 0 : 1;
    bool curr_match;
    for(int i=start_softclipp+1; i < (end_softclipp+1); i++) {
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
//                needs_aln = counter > 2 ? true: false;
//            }
            counter = 0;
        }
        prev_is_match = curr_match;
        counter++;
    }

    // Print last
    if ( curr_match  ){
        cigar << counter << '=';
        aln_score += counter * match;
    } else{
        cigar << counter << 'X';
        hamming_mod += counter;
//        std::cout << "Added2: " << counter << " current: " << hamming_mod << std::endl;

//        needs_aln = counter > 2 ? true: false;
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


inline bool HammingToCigarEQX(std::string &One, std::string &Two, std::stringstream &cigar)
{
    if (One.length() != Two.length()){
        return true;
    }

    int counter = 1;
    bool prev_is_match = One[0] == Two[0] ? true : false;
    bool beginning = true;
    bool needs_aln = false; // prev_is_match == true ? false : true;
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
                needs_aln = counter > 2 ? true: false;
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
        needs_aln = counter > 2 ? true: false;
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


static inline void align_SE(alignment_params &aln_params, std::string &sam_string, std::vector<nam> &all_nams, std::string &query_acc, idx_to_acc &acc_map, int k, int read_len, std::vector<unsigned int> &ref_len_map, std::vector<std::string> &ref_seqs, std::string &read, std::string &qual, unsigned int &tot_ksw_aligned, unsigned int &tot_all_tried, float dropoff, unsigned int &did_not_fit, int max_tries ) {

    std::string read_rc;
    bool rc_already_comp = false;

    if (all_nams.size() == 0) {
        sam_string.append(query_acc);
        sam_string.append("\t4\t*\t0\t255\t*\t*\t0\t0\t");
        sam_string.append(read);
        sam_string.append("\t*\n");
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

        tot_all_tried ++;

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
            did_not_fit++;
            aln_did_not_fit = true;
        }

        // deal with any read hanging of ends of reference not to get 'std::out_of_range' what(): basic_string::substr
        int ref_tmp_start = n.ref_s - n.query_s;
        int ref_tmp_segm_size = read_len + diff;
        int ref_len = ref_len_map[n.ref_id];
        int ref_start = ref_tmp_start > 0 ? ref_tmp_start : 0;
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
                int diff_to_best = sw_score < best_align_sw_score ? best_align_sw_score - sw_score : sw_score - best_align_sw_score;
                min_mapq_diff = min_mapq_diff < diff_to_best ? min_mapq_diff : diff_to_best;
                std::stringstream cigar_string;
                int aln_score = 0;
                hamming_mod = HammingToCigarEQX2(r_tmp, ref_segm, cigar_string, aln_params.match, aln_params.mismatch, aln_score, soft_left, soft_right);
//                if (hamming_dist < best_align_dist){
                if (aln_score > best_align_sw_score){
                    min_mapq_diff = (sw_score - best_align_sw_score) > 0 ? (sw_score - best_align_sw_score)  : 0 ; // new distance to next best match

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
            int extra_ref_left = soft_left <= 50 ? soft_left : 50;
            int extra_ref_right = soft_right <= 50 ? soft_right: 50;
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
            info = ssw_align(ref_segm, r_tmp, read_len, aln_params.match, aln_params.mismatch, aln_params.gap_open, aln_params.gap_extend);
//            info.ed = info.global_ed; // read_len - info.sw_score;
            int diff_to_best = info.sw_score < best_align_sw_score ? best_align_sw_score - info.sw_score : info.sw_score - best_align_sw_score;
            min_mapq_diff = min_mapq_diff < diff_to_best ? min_mapq_diff : diff_to_best;
            tot_ksw_aligned ++;
//            if (info.global_ed <= best_align_dist){
            if (info.sw_score >= best_align_sw_score){
                min_mapq_diff = (info.sw_score - best_align_sw_score) > 0 ? (info.sw_score - best_align_sw_score)  : 0 ; // new distance to next best match
                best_align_dist = info.global_ed;
                sam_aln.cigar = info.cigar;
                sam_aln.ed = info.ed;
                sam_aln.ref_start =  ref_start + info.ref_offset +1; // +1 because SAM is 1-based!
                sam_aln.is_rc = is_rc;
                sam_aln.ref_id = n.ref_id;
                sam_aln.sw_score = info.sw_score;
                best_align_sw_score = info.sw_score;
                sam_aln.aln_score = info.sw_score;
            }
//            std::cout << "HERE 2 "  << sam_aln.ref_start << " global ed: " << info.global_ed  <<  ", hamming: " << hamming_dist << ", "<< read_len << ", " << (float) hamming_dist / (float) read_len << std::endl;

//            std::cout << "Aligned: " << n.score << ", "  << n.n_hits << ", " << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e  << ") ed:" << info.ed << ", best ed so far: " << best_align_dist  << std::endl;

        }
        cnt ++;
    }

    if (all_nams.size() > 0) {
        int o;
        std::string output_read;
        if (sam_aln.is_rc) {
            o = 16;
            output_read = read_rc;
        } else {
            o = 0;
            output_read = read;
        }
        sam_aln.mapq = min_mapq_diff <= 60 ? min_mapq_diff : 60;
        sam_string.append(query_acc);
        sam_string.append("\t");
        sam_string.append(std::to_string(o));
        sam_string.append("\t");
        sam_string.append(acc_map[sam_aln.ref_id]);
        sam_string.append("\t");
        sam_string.append(std::to_string(sam_aln.ref_start));
        sam_string.append("\t");
        sam_string.append(std::to_string(sam_aln.mapq));
        sam_string.append("\t");
        sam_string.append(sam_aln.cigar);
        sam_string.append("\t*\t0\t0\t");
        if (!sam_aln.is_unaligned) {
            sam_string.append(output_read);
            sam_string.append("\t");
            if (sam_aln.is_rc){
                auto qual_rev = qual;
                std::reverse(qual_rev.begin(), qual_rev.end()); // reverse
                sam_string.append(qual_rev);
            } else {
                sam_string.append(qual);
            }
            sam_string.append("\t");
            sam_string.append("NM:i:");
            sam_string.append(std::to_string(sam_aln.ed));
            sam_string.append("\t");
            sam_string.append("AS:i:");
            sam_string.append(std::to_string((int) sam_aln.aln_score));
        } else {
            sam_string.append(read);
            sam_string.append("\t");
            sam_string.append(qual);
        }
        sam_string.append("\n");

//        sam_string.append("\t*\tNM:i:");
//        sam_string.append(std::to_string(sam_aln.ed));
//        sam_string.append("\n");
    }
}


static inline void align_SE_secondary_hits(alignment_params &aln_params, std::string &sam_string, std::vector<nam> &all_nams, std::string &query_acc, idx_to_acc &acc_map, int k, int read_len, std::vector<unsigned int> &ref_len_map, std::vector<std::string> &ref_seqs, std::string &read, std::string &qual, unsigned int &tot_ksw_aligned, unsigned int &tot_all_tried, float dropoff, unsigned int &did_not_fit, int max_tries, int max_secondary ) {

    std::string read_rc;
    bool rc_already_comp = false;

    if (all_nams.size() == 0) {
        sam_string.append(query_acc);
        sam_string.append("\t4\t*\t0\t255\t*\t*\t0\t0\t");
        sam_string.append(read);
        sam_string.append("\t*\n");
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

        tot_all_tried ++;

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
            did_not_fit++;
            aln_did_not_fit = true;
        }

        // deal with any read hanging of ends of reference not to get 'std::out_of_range' what(): basic_string::substr
        int ref_tmp_start = n.ref_s - n.query_s;
        int ref_tmp_segm_size = read_len + diff;
        int ref_len = ref_len_map[n.ref_id];
        int ref_start = ref_tmp_start > 0 ? ref_tmp_start : 0;
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
                int diff_to_best = sw_score < best_align_sw_score ? best_align_sw_score - sw_score : sw_score - best_align_sw_score;
                min_mapq_diff = min_mapq_diff < diff_to_best ? min_mapq_diff : diff_to_best;
//                if (hamming_dist < best_align_dist) {
//                    min_mapq_diff = best_align_dist - hamming_dist; // new distance to next best match
//                }
                if (sw_score > best_align_sw_score){
                    min_mapq_diff = (sw_score - best_align_sw_score) > 0 ? (sw_score - best_align_sw_score)  : 0 ; // new distance to next best match
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
            int extra_ref_left = soft_left <= 50 ? soft_left : 50;
            int extra_ref_right = soft_right <= 50 ? soft_right: 50;
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
            info = ssw_align(ref_segm, r_tmp, read_len, aln_params.match, aln_params.mismatch, aln_params.gap_open, aln_params.gap_extend);
//            info.ed = info.global_ed; // read_len - info.sw_score;
            sw_score = info.sw_score;
            tot_ksw_aligned ++;
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
//            std::cout << "Aligned: " << n.score << ", "  << n.n_hits << ", " << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e  << ") ed:" << info.ed << ", best ed so far: " << best_align_dist  << std::endl;

        }
        if (sw_score > best_align_sw_score){
            best_align_sw_score = sw_score;
        }

        std::tuple<double, alignment> t (sam_aln.sw_score, sam_aln);
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
            int f = 0;
            if (i > 0) {
                f += 256;
                sam_aln.mapq = 255;
            } else {
                sam_aln.mapq = min_mapq_diff < 60 ? min_mapq_diff : 60;
            }


            std::string output_read;
            if (sam_aln.is_rc) {
                f += 16;
                output_read = read_rc;
            } else {
                output_read = read;
            }
            sam_string.append(query_acc);
            sam_string.append("\t");
            sam_string.append(std::to_string(f));
            sam_string.append("\t");
            sam_string.append(acc_map[sam_aln.ref_id]);
            sam_string.append("\t");
            sam_string.append(std::to_string(sam_aln.ref_start));
            sam_string.append("\t");
            sam_string.append(std::to_string(sam_aln.mapq));
            sam_string.append("\t");
            sam_string.append(sam_aln.cigar);
            sam_string.append("\t*\t0\t0\t");
            if (!sam_aln.is_unaligned) {
                sam_string.append(output_read);
                sam_string.append("\t");
                if (sam_aln.is_rc){
                    auto qual_rev = qual;
                    std::reverse(qual_rev.begin(), qual_rev.end()); // reverse
                    sam_string.append(qual_rev);
                } else {
                    sam_string.append(qual);
                }
                sam_string.append("\t");
                sam_string.append("NM:i:");
                sam_string.append(std::to_string(sam_aln.ed));
                sam_string.append("\t");
                sam_string.append("AS:i:");
                sam_string.append(std::to_string((int) sam_aln.aln_score));
            } else {
                sam_string.append(read);
                sam_string.append("\t");
                sam_string.append(qual);
            }
            sam_string.append("\n");
//            sam_string.append("\t*\tNM:i:");
//            sam_string.append(std::to_string(sam_aln.ed));
//            sam_string.append("\n");
        }
    }
}
//
//static inline void align_segment(alignment_params &aln_params, std::string &read_segm, std::string &ref_segm, int read_segm_len, int ref_segm_len, int ref_start,  int ext_left, int ext_right, bool aln_did_not_fit, bool is_rc, alignment &sam_aln_segm, unsigned int &tot_ksw_aligned) {
//    int hamming_dist = -1;
//    int soft_left = 50;
//    int soft_right = 50;
//    int hamming_mod;
//    int ref_segm_len_ham = ref_segm_len - ext_left - ext_right; // we send in the already extended ref segment to save time. This is not true in center alignment if merged match have diff length
////    std::cout << "ref_segm_len_ham" << ref_segm_len_ham << std::endl;
//    if ( (ref_segm_len_ham == read_segm_len) && (!aln_did_not_fit) ){
//        std::string ref_segm_ham = ref_segm.substr(ext_left, read_segm_len);
////        std::cout << "ref_segm_ham " << ref_segm_ham << std::endl;
//
//        hamming_dist = HammingDistance(read_segm, ref_segm_ham);
////        std::cout << "hamming_dist " << hamming_dist << std::endl;
//
//        if ( (hamming_dist >= 0) && (((float) hamming_dist / read_segm_len) < 0.05) ) { //Hamming distance worked fine, no need to ksw align
//            std::stringstream cigar_string;
//            int aln_score = 0;
//            hamming_mod = HammingToCigarEQX2(read_segm, ref_segm_ham, cigar_string, aln_params.match, aln_params.mismatch, aln_score, soft_left, soft_right);
//            sam_aln_segm.cigar = cigar_string.str();
//            sam_aln_segm.ed = hamming_mod;
//            sam_aln_segm.sw_score = aln_score; // aln_params.match*(read_len-hamming_dist) - aln_params.mismatch*hamming_dist;
//            sam_aln_segm.ref_start = ref_start + ext_left + soft_left+1; // +1 because SAM is 1-based!
//            sam_aln_segm.is_rc = is_rc;
//            sam_aln_segm.is_unaligned = false;
//            sam_aln_segm.aln_score = aln_score;
////            std::cout << "cigar " << sam_aln_segm.cigar << " sam_aln_segm.ref_start " << sam_aln_segm.ref_start << std::endl;
//            return;
//        }
//    }
//
//    aln_info info;
//    info = ssw_align(ref_segm, read_segm, read_segm_len, aln_params.match, aln_params.mismatch, aln_params.gap_open, aln_params.gap_extend);
//    tot_ksw_aligned ++;
//    sam_aln_segm.cigar = info.cigar;
//    sam_aln_segm.ed = info.ed;
////    std::cerr << r_tmp << " " << n.n_hits << " " << n.score << " " <<  diff << " " << sam_aln.ed << " "  <<  n.query_s << " "  << n.query_e << " "<<  n.ref_s << " "  << n.ref_e << " " << n.is_rc << " " << hamming_dist << " " << sam_aln.cigar << " " << info.sw_score << std::endl;
//    sam_aln_segm.sw_score = info.sw_score;
//    sam_aln_segm.ref_start =  ref_start + info.ref_offset +1; // +1 because SAM is 1-based!
//    sam_aln_segm.is_rc = is_rc;
//    sam_aln_segm.is_unaligned = false;
//    sam_aln_segm.aln_score = info.sw_score;
////    std::cerr << " ALIGN SCORE: " << sam_aln_segm.sw_score << " cigar: " << sam_aln_segm.cigar << std::endl;
//}
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
////    std::cerr << "n.ID " << n.nam_id  << " n.n_hits " << n.n_hits << " n.ref_s " <<  n.ref_s <<  " n.ref_e " << n.ref_e << " read " << read << std::endl;
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
//    int ext_left;
//    int ext_right;
//    int ref_tmp_segm_size;
//    int ref_len = ref_len_map[n.ref_id];
//    int ref_segm_size;
//    int ref_tmp_end;
//    int ref_tmp_start;
//    int ref_end;
//    int ref_start;
//    std::string ref_segm;
//    std::string read_segm;
//
//    if (true){ // full alignment
//        ref_tmp_start = n.ref_s - n.query_s > 0 ? n.ref_s - n.query_s : 0;
//        ext_left = ref_tmp_start < 50 ? ref_tmp_start : 50;
//        ref_start = ref_tmp_start - ext_left;
//
//        ref_tmp_segm_size = read_len + diff;
//        ext_right = ref_len - (n.ref_e +1) < 50 ? ref_len - (n.ref_e +1) : 50;
//
//        ref_segm_size = ref_tmp_segm_size + ext_left + ext_right;
//        ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_segm_size);
////        std::cerr << " ref_tmp_start " << ref_tmp_start << " ext left " << ext_left << " ext right " << ext_right << " ref_tmp_segm_size " << ref_tmp_segm_size << " ref_segm_size " << ref_segm_size << " ref_segm " << ref_segm << std::endl;
//        sam_aln.ref_id = n.ref_id;
//        align_segment(aln_params, r_tmp, ref_segm, read_len, ref_segm_size, ref_start, ext_left, ext_right, aln_did_not_fit, is_rc, sam_aln, tot_ksw_aligned);
//    } else{
//        // test full hamming based alignment first
//        ref_tmp_start = n.ref_s - n.query_s > 0 ? n.ref_s - n.query_s : 0;
//        int ref_start = ref_tmp_start > 0 ? ref_tmp_start : 0;
//        ref_tmp_segm_size = read_len + diff;
//        ref_segm_size = ref_tmp_segm_size < ref_len - ref_start ? ref_tmp_segm_size : ref_len - 1 - ref_start;
//
//        std::string ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_segm_size);
//        if ( (ref_segm_size == read_len) && fits ){
//            std::cout << "ref_segm_full_ham " << ref_segm << std::endl;
//
//            int hamming_dist = HammingDistance(r_tmp, ref_segm);
//            std::cout << "hamming_dist " << hamming_dist << std::endl;
//
//            if ( (hamming_dist >= 0) && (((float) hamming_dist / ref_segm_size) < 0.05) ) { //Hamming distance worked fine, no need to ksw align
//                std::stringstream cigar_string;
//                int aln_score = 0;
//                int soft_left = 0;
//                int soft_right = 0;
//                int hamming_mod = HammingToCigarEQX2(r_tmp, ref_segm, cigar_string, aln_params.match, aln_params.mismatch, aln_score, soft_left, soft_right);
//                sam_aln.cigar = cigar_string.str();
//                sam_aln.ed = hamming_mod;
//                sam_aln.sw_score = aln_score; // aln_params.match*(read_len-hamming_dist) - aln_params.mismatch*hamming_dist;
//                sam_aln.ref_start = ref_start + ext_left + soft_left+1; // +1 because SAM is 1-based!
//                sam_aln.is_rc = is_rc;
//                sam_aln.is_unaligned = false;
//                sam_aln.aln_score = aln_score;
//                std::cout << "FULL HAMMING , returning " << sam_aln.cigar << " sam_aln_segm.ref_start " << sam_aln.ref_start << std::endl;
//                return;
//            }
//        }
//
//        //// Didn't work with global Hamming - split into parts
//
//        // left alignment
//        ref_end = n.ref_s + k;
//        ext_left = ref_end < 50 ? ref_end : 50;
//        ref_start = ref_end - ext_left;
//
//        ref_tmp_segm_size = n.query_s + k;
//        ext_right = 0;
//
//        ref_segm_size = ref_tmp_segm_size + ext_left + ext_right;
//        ref_segm = ref_seqs[n.ref_id].substr(ref_end - ext_left, ref_segm_size);
//        alignment sam_aln_segm_left;
//        sam_aln_segm_left.ref_id = n.ref_id;
//        read_segm = r_tmp.substr(0, n.query_s+k);
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
//
//    }
//
//    }




static inline void get_alignment(alignment_params &aln_params, nam &n, std::vector<unsigned int> &ref_len_map, std::vector<std::string> &ref_seqs, std::string &read, std::string &read_rc, int read_len, alignment &sam_aln, int k, int cnt, bool &rc_already_comp, unsigned int &did_not_fit, unsigned int &tot_ksw_aligned){
    bool aln_did_not_fit = false;
    int ref_diff = n.ref_e - n.ref_s;
    int read_diff = n.query_e - n.query_s;
    int min_diff =  read_diff ^ ((ref_diff ^ read_diff) & -(ref_diff < read_diff));
    int max_diff = ref_diff ^ ((ref_diff ^ read_diff) & -(ref_diff < read_diff));
    int diff = max_diff - min_diff;
//    int max_allowed_mask = aln_params.gap_open/aln_params.match - 1 > 0 ? aln_params.gap_open/aln_params.match - 1 : 1;

    // deal with any read hanging of ends of reference not to get 'std::out_of_range' what(): basic_string::substr
    int ref_tmp_start = n.ref_s - n.query_s;
    int ref_tmp_segm_size = read_len + diff;
    int ref_len = ref_len_map[n.ref_id];
    int ref_start = ref_tmp_start > 0 ? ref_tmp_start : 0;
    int ref_segm_size = ref_tmp_segm_size < ref_len - ref_start ? ref_tmp_segm_size : ref_len - 1 - ref_start;

    std::string ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_segm_size);

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

//    std::cerr<< r_tmp << std::endl;
//    std::cerr<< ref_segm << std::endl;
//    std::cerr<< diff << std::endl;
    int soft_left = 50;
    int soft_right = 50;
    int hamming_mod;
//    bool needs_aln = false;
    if ( (ref_segm_size == read_len) && (!aln_did_not_fit) ){
        hamming_dist = HammingDistance(r_tmp, ref_segm);
//        std::cerr<< "Here " << hamming_dist << std::endl;
//        std::cerr<< aln_params.gap_open/aln_params.match  << std::endl;
        if ( (hamming_dist >= 0) && (((float) hamming_dist / read_len) < 0.05) ) { //Hamming distance worked fine, no need to ksw align
            std::stringstream cigar_string;
//            needs_aln = HammingToCigarEQX(r_tmp, ref_segm, cigar_string);
            int aln_score = 0;
            hamming_mod = HammingToCigarEQX2(r_tmp, ref_segm, cigar_string, aln_params.match, aln_params.mismatch, aln_score, soft_left, soft_right);

//            needs_aln = false;
            sam_aln.cigar = cigar_string.str();
//            sam_aln.cigar = std::to_string(read_len) + "M";
//            std::cerr<< "Here ham dist: " << hamming_dist << " ham mod: " << hamming_mod << " " << r_tmp.size() << " " << ref_segm.size()  << std::endl;
            sam_aln.ed = hamming_mod;
//            sam_aln.sw_score = aln_score;
            sam_aln.sw_score = aln_score; // aln_params.match*(read_len-hamming_dist) - aln_params.mismatch*hamming_dist;
            sam_aln.ref_start = ref_start + soft_left+1; // +1 because SAM is 1-based!
            sam_aln.is_rc = is_rc;
            sam_aln.ref_id = n.ref_id;
            sam_aln.is_unaligned = false;
            sam_aln.aln_score = aln_score;
            return;
//            if (hamming_mod == hamming_dist ){ // masked only what is justified by alingment parameters max_allowed_mask
//                return;
//            }
        }
        //TODO: Only do ksw of the ends outside the NAM to increase speed here
//        else{ // Segment(s) of read outside the NAM span is not fitting to reference, align the segments
//            std::cerr<< sam_aln.ed << " " << sam_aln.sw_score << " " <<   n.query_s << " " << n.query_e << std::endl;
//            std::cerr<< r_tmp << std::endl;
//            std::cerr<< ref_segm.substr(0,read_len) << std::endl;
//
//        }
    }

    // We didn't get away with hamming distance, do full ksw alignment
//    else {
//    std::cerr<< "3" << std::endl;

    int extra_ref_left = soft_left <= 50 ? soft_left : 50;
    int extra_ref_right = soft_right <= 50 ? soft_right: 50;
    int a = n.ref_s - n.query_s - extra_ref_left;
    ref_start = std::max(0, a);
    int b = n.ref_e + (read_len - n.query_e)+ extra_ref_right;
    int ref_end = std::min(ref_len, b);
    ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_end - ref_start);
//    ksw_extz_t ez;
//    const char *ref_ptr = ref_segm.c_str();
//    const char *read_ptr = r_tmp.c_str();
    aln_info info;
//    std::cerr<< "4" << std::endl;
//    info = ksw_align(ref_ptr, ref_segm.size(), read_ptr, r_tmp.size(), 1, 4, 6, 1, ez);
    info = ssw_align(ref_segm, r_tmp, read_len, aln_params.match, aln_params.mismatch, aln_params.gap_open, aln_params.gap_extend);
//    if (info.ed == 100000){
//        std::cerr<< "________________________________________" << std::endl;
//        std::cerr<< "NORMAL MODE" << std::endl;
//        std::cerr<< read << "   " << read_rc << std::endl;
//        std::cerr << r_tmp << " " << n.n_hits << " " << n.score << " " <<  diff << " " << sam_aln.ed << " "  <<  n.query_s << " "  << n.query_e << " "<<  n.ref_s << " "  << n.ref_e << " " << n.is_rc << " " << hamming_dist << " " << sam_aln.cigar << " " << info.sw_score << std::endl;
//        std::cerr << "a " << a << " b " << b << " ref_start " <<  ref_start << " ref_end " << ref_end << "  ref_end - ref_start "  <<  ref_end - ref_start << " extra_ref_left "  << extra_ref_left << " extra_ref_right "<<  extra_ref_right << "  n.is_flipped " <<  n.is_flipped << std::endl;
//        std::cerr<< "________________________________________" << std::endl;
//    }

//    std::cerr<< "5" << std::endl;
    sam_aln.cigar = info.cigar;
    sam_aln.ed = info.ed;
//    std::cerr << r_tmp << " " << n.n_hits << " " << n.score << " " <<  diff << " " << sam_aln.ed << " "  <<  n.query_s << " "  << n.query_e << " "<<  n.ref_s << " "  << n.ref_e << " " << n.is_rc << " " << hamming_dist << " " << sam_aln.cigar << " " << info.sw_score << std::endl;
    sam_aln.sw_score = info.sw_score;
    sam_aln.ref_start =  ref_start + info.ref_offset +1; // +1 because SAM is 1-based!
    sam_aln.is_rc = is_rc;
    sam_aln.ref_id = n.ref_id;
    sam_aln.is_unaligned = false;
    sam_aln.aln_score = info.sw_score;
    tot_ksw_aligned ++;
//    }
}

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

static inline void append_to_sam(std::string &sam_string, alignment &sam_aln1, alignment &sam_aln2, std::string &read1, std::string &read2, std::string &read1_rc, std::string &read2_rc, idx_to_acc &acc_map, std::string &query_acc1, std::string &query_acc2, int &mapq1, int &mapq2, float &mu, float &sigma, int read_len, bool is_primary, std::string &qual1, std::string &qual2){
    int f1 = 1;
    int f2 = 1; // template having multiple segments in sequencing
//    std::cerr <<  query_acc1 << " HERE 1: " << sam_aln1.is_unaligned << " " << query_acc2 << " HERE 2: " << sam_aln2.is_unaligned  << std::endl;
    if (!is_primary){
        f1 |= (1u << 8); // not primary
        f2 |= (1u << 8); // not primary
    }

    // Commented lines below because we do not longer mark a read as not proper just because of the non-matching hash
    // Proper or non proper reads are further below only decided based on the expected distance and relative orientation they align to
//    if (sam_aln1.ed < 5){ // Flag alignments previously deemed as 'not proper' (based on matching strobemer hash ) to proper because of small ed
//        sam_aln1.not_proper = false;
//    }
//    if (sam_aln2.ed < 5){ // Flag alignments previously deemed as 'not proper' (based on matching strobemer hash ) to proper because of small ed
//        sam_aln2.not_proper = false;
//    }

    int d, template_len1, template_len2;
    if (sam_aln1.ref_start < sam_aln2.ref_start){
        d = sam_aln2.ref_start - sam_aln1.ref_start;
        template_len1 = - d - read_len;
        template_len2 = d + read_len;
    }
    else{
        d = sam_aln1.ref_start - sam_aln2.ref_start;
        template_len1 = d + read_len;
        template_len2 = - d - read_len;
    }
    //    int d = sam_aln1.ref_start < sam_aln2.ref_start ? sam_aln2.ref_start - sam_aln1.ref_start : sam_aln1.ref_start - sam_aln2.ref_start;

    bool both_aligned = !sam_aln1.is_unaligned && !sam_aln2.is_unaligned;
    int a = sam_aln2.ref_start - sam_aln1.ref_start;
    bool r1_r2 = !sam_aln1.is_rc && sam_aln2.is_rc && a >= 0; // r1 ---> <---- r2
    bool r2_r1 = !sam_aln2.is_rc && sam_aln1.is_rc && a <= 0 ; // r2 ---> <---- r1
    bool rel_orientation_good = r1_r2 || r2_r1;
    bool insert_good = d <= mu + 6*sigma;
    if ( both_aligned && insert_good && rel_orientation_good ){
        sam_aln1.not_proper = false;
        sam_aln2.not_proper = false;
    } else {
        sam_aln1.not_proper = true;
        sam_aln2.not_proper = true;
    }

//    if ( d > (mu + 6*sigma) ){ // Flag alignments as 'not proper' because of too large deviation in insert size
//            sam_aln1.not_proper = true;
//            sam_aln2.not_proper = true;
//        }

    if ( (!sam_aln1.not_proper) && (!sam_aln2.not_proper)){ // if both segements in pair are properly aligned
        f1 |= (1u << 1);
        f2 |= (1u << 1);
    }

    std::string output_read1;
    output_read1 = read1;
    std::string output_read2;
    output_read2 = read2;
    f1 |= (1u << 6); // first segment in template
    f2 |= (1u << 7); // last segment in template
    if (sam_aln1.is_rc){ // set if alignnment1 is reverse complemented
        f1 |= (1u << 4);
        f2 |= (1u << 5);
        output_read1 = read1_rc;
    }
    if (sam_aln2.is_rc){ // set if alignnment2 is reverse complemented
        f1 |= (1u << 5);
        f2 |= (1u << 4);
        output_read2 = read2_rc;
    }

    std::string m1_chr;
    std::string m2_chr;
    if (sam_aln1.ref_id == sam_aln2.ref_id){
        m1_chr = "=";
        m2_chr = "=";
    } else{
        m1_chr = acc_map[sam_aln1.ref_id];
        m2_chr = acc_map[sam_aln2.ref_id];
    }

//    if ( (sam_aln1.is_unaligned) && (sam_aln2.is_unaligned) ){
//        f1 = 13;
//        f2 = 13;
//        m1_chr = "*";
//        m2_chr = "*";
//        sam_aln1.cigar = "*";
//        sam_aln2.cigar = "*";
//    } else if (sam_aln1.is_unaligned){
//        f1 = 5;
//        m1_chr = "*";
//        sam_aln1.cigar = "*";
//        f2 |= (1u << 3);
//        f2 -= 32;
//    } else if (sam_aln2.is_unaligned){
//        f2 = 5;
//        m2_chr = "*";
//        sam_aln2.cigar = "*";
//        f1 |= (1u << 3);
//        f1 -= 32;
//    }
    std::string ref1 = acc_map[sam_aln1.ref_id];
    std::string ref2 = acc_map[sam_aln2.ref_id];
    int ed1 = sam_aln1.ed;
    int ed2 = sam_aln2.ed;

    if (sam_aln1.is_unaligned && sam_aln2.is_unaligned){
        f1 |= (1u << 2);
        f1 |= (1u << 3);
        f2 |= (1u << 2);
        f2 |= (1u << 3);
        sam_aln1.ref_start = 0;
        sam_aln2.ref_start = 0;
        template_len1 = 0;
        template_len2 = 0;
        ref1 = "*";
        ref2 = "*";
        f1 |= (0u << 4);
        f1 |= (0u << 5);
        f2 |= (0u << 4);
        f2 |= (0u << 5);
        ed1 = 0;
        ed2 = 0;
        mapq1 = 255;
        mapq2 = 255;
    } else if (sam_aln1.is_unaligned){
        f1 |= (1u << 2);
        f1 |= (0u << 4);
        f2 |= (1u << 3);
        sam_aln1.ref_start = sam_aln2.ref_start;
        template_len1 = 0;
        template_len2 = 0;
        ed1 = 0;
        mapq1 = 255;
    } else if (sam_aln2.is_unaligned){
        f2 |= (1u << 2);
        f2 |= (0u << 4);
        f1 |= (1u << 3);
        sam_aln2.ref_start = sam_aln1.ref_start;
        template_len1 = 0;
        template_len2 = 0;
        ed2 = 0;
        mapq2 = 255;
    }

//    if ( (sam_aln1.ref_start == 0)){ // && !sam_aln1.is_unaligned  ){
//        std::cerr << "OMG1" << std::endl;
//        std::cerr << query_acc1 << std::endl;
//    }
//    if ( (sam_aln2.ref_start == 0)){ // && !sam_aln2.is_unaligned  ){
//        std::cerr << "OMG2" << std::endl;
//        std::cerr << query_acc2 << std::endl;
//    }

    sam_string.append(query_acc1);
    sam_string.append("\t");
    sam_string.append(std::to_string(f1));
    sam_string.append("\t");
    sam_string.append(ref1);
    sam_string.append("\t");
    sam_string.append(std::to_string(sam_aln1.ref_start));
    sam_string.append("\t");
    sam_string.append(std::to_string(mapq1));
    sam_string.append("\t");
    sam_string.append(sam_aln1.cigar);
    sam_string.append("\t");
    sam_string.append(m2_chr);
    sam_string.append("\t");
    sam_string.append(std::to_string(sam_aln2.ref_start));
    sam_string.append("\t");
    sam_string.append(std::to_string(template_len1));
    sam_string.append("\t");
//    sam_string.append("\t*\t0\t0\t");
    if (!sam_aln1.is_unaligned) {
        sam_string.append(output_read1);
        sam_string.append("\t");
        if (sam_aln1.is_rc){
            auto qual_rev = qual1;
            std::reverse(qual_rev.begin(), qual_rev.end()); // reverse
            sam_string.append(qual_rev);
        } else {
            sam_string.append(qual1);
        }
        sam_string.append("\t");
        sam_string.append("NM:i:");
        sam_string.append(std::to_string(ed1));
        sam_string.append("\t");
        sam_string.append("AS:i:");
        sam_string.append(std::to_string((int) sam_aln1.aln_score));
    } else {
        sam_string.append(read1);
        sam_string.append("\t");
        sam_string.append(qual1);
    }
    sam_string.append("\n");

    sam_string.append(query_acc2);
    sam_string.append("\t");
    sam_string.append(std::to_string(f2));
    sam_string.append("\t");
    sam_string.append(ref2);
    sam_string.append("\t");
    sam_string.append(std::to_string(sam_aln2.ref_start));
    sam_string.append("\t");
    sam_string.append(std::to_string(mapq2));
    sam_string.append("\t");
    sam_string.append(sam_aln2.cigar);
    sam_string.append("\t");
    sam_string.append(m1_chr);
    sam_string.append("\t");
    sam_string.append(std::to_string(sam_aln1.ref_start));
    sam_string.append("\t");
    sam_string.append(std::to_string(template_len2));
    sam_string.append("\t");
//    sam_string.append("\t*\t0\t0\t");
    if (!sam_aln2.is_unaligned) {
        sam_string.append(output_read2);
        sam_string.append("\t");
        if (sam_aln2.is_rc){
            auto qual2_rev = qual2;
            std::reverse(qual2_rev.begin(), qual2_rev.end()); // reverse
            sam_string.append(qual2_rev);
        } else {
            sam_string.append(qual2);
        }
        sam_string.append("\t");
        sam_string.append("NM:i:");
        sam_string.append(std::to_string(ed2));
        sam_string.append("\t");
        sam_string.append("AS:i:");
        sam_string.append(std::to_string((int) sam_aln2.aln_score));
    } else {
        sam_string.append(read2);
        sam_string.append("\t");
        sam_string.append(qual2);
    }
    sam_string.append("\n");
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
    float x;
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


static inline void rescue_mate(alignment_params &aln_params , nam &n, std::vector<unsigned int> &ref_len_map, std::vector<std::string> &ref_seqs, std::string &read, std::string &read_rc, int read_len, alignment &sam_aln, bool &rc_already_comp, unsigned int &tot_ksw_aligned, float &mu, float &sigma, unsigned int &tot_rescued, int k) {
    int a, b, ref_start,ref_len,ref_end;
    std::string r_tmp;
    bool a_is_rc;

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
    ref_end = std::min(ref_len, b);
    std::string ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_end - ref_start);
    aln_info info;

    // check that read shares at least some segment with ref otherwise abort
    int sub_size = 2*k/3;
    int step_size = k/3;
    std::string submer;
    bool found = false;
    for (int i = 0; i<=r_tmp.size()-k; i+=step_size) {
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
//        std::cerr << "LOOOOOOL!" << std::endl;
//        std::cerr << "Aligning anyway at: " << ref_start << " to " << ref_end << "ref len:" << ref_len << " ref_id:" << n.ref_id << std::endl;
//        std::cerr << "read: " << r_tmp << std::endl;
//        std::cerr << "ref: " << ref_segm << std::endl;
    }

//    std::cerr << "Aligning at: " << ref_start << " to " << ref_end << "ref len:" << ref_len << " ref_id:" << n.ref_id << std::endl;
//    std::cerr << "read: " << r_tmp << std::endl;
//    std::cerr << "ref: " << ref_segm << std::endl;
    info = ssw_align(ref_segm, r_tmp, read_len, aln_params.match, aln_params.mismatch, aln_params.gap_open, aln_params.gap_extend);

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
    tot_ksw_aligned ++;
    tot_rescued ++;
}

static inline void align_PE(alignment_params &aln_params, std::string &sam_string, std::vector<nam> &all_nams1, std::vector<nam> &all_nams2, KSeq &record1, KSeq &record2, idx_to_acc &acc_map, int k, std::vector<unsigned int> &ref_len_map, std::vector<std::string> &ref_seqs, unsigned int &tot_ksw_aligned, unsigned int &tot_all_tried, unsigned int &tot_rescued, float dropoff, unsigned int &did_not_fit, float &mu, float &sigma, float &sample_size, float &V, float &SSE, int max_tries, int max_secondary) {
    int read_len1 = record1.seq.length();
    int read_len2 = record2.seq.length();
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
//        std::cerr << "unaligned" << std::endl;
        sam_string.append(query_acc1);
        sam_string.append("\t77\t*\t0\t255\t*\t*\t0\t0\t");
        sam_string.append(read1);
        sam_string.append("\t*\n");
        sam_string.append(query_acc2);
        sam_string.append("\t141\t*\t0\t255\t*\t*\t0\t0\t");
        sam_string.append(read2);
        sam_string.append("\t*\n");
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
        if ( (score_dropoff1 < dropoff) && (score_dropoff2 < dropoff) && (n_max1.is_rc ^ n_max2.is_rc) && ( r1_r2 || r1_r2 ) ){ //( ((n_max1.ref_s - n_max2.ref_s) < mu + 4*sigma ) || ((n_max2.ref_s - n_max1.ref_s ) < mu + 4*sigma ) ) &&
//            std::cerr << "I'm here" << std::endl;
//            std::cerr << query_acc1 << std::endl;
            get_alignment(aln_params, n_max1, ref_len_map, ref_seqs, read1, read1_rc, read_len1, sam_aln1, k, cnt1, rc_already_comp1, did_not_fit, tot_ksw_aligned);
            tot_all_tried ++;
//            std::cerr << query_acc2 << std::endl;
            get_alignment(aln_params, n_max2, ref_len_map, ref_seqs, read2, read2_rc, read_len2, sam_aln2, k, cnt2, rc_already_comp2, did_not_fit, tot_ksw_aligned);
            tot_all_tried ++;
//            std::cerr<< "6" << std::endl;
            get_MAPQ(all_nams1, n_max1, mapq1);
            get_MAPQ(all_nams2, n_max2, mapq2);
//            std::cerr<< "7" << std::endl;
            append_to_sam(sam_string,sam_aln1, sam_aln2, read1, read2, read1_rc, read2_rc, acc_map, query_acc1, query_acc2, mapq1, mapq2, mu, sigma, read_len1, true, record1.qual, record2.qual);

            if ((sample_size < 400) && ((sam_aln1.ed + sam_aln2.ed) < 3) && !sam_aln1.not_proper && !sam_aln2.not_proper ){
                int d = sam_aln1.ref_start > sam_aln2.ref_start ? sam_aln1.ref_start - sam_aln2.ref_start : sam_aln2.ref_start - sam_aln1.ref_start;
                if ( d < 2000){
//                    std::cerr<< "8 " << sample_size << std::endl;
                    float e;
                    e = d - mu;
                    mu = mu + e/sample_size; // (1.0/(sample_size +1.0)) * (sample_size*mu + d);
                    SSE = SSE + e*(d-mu);
                    V = sample_size > 1 ? SSE/(sample_size -1.0) : SSE; //d < 1000 ? ((sample_size +1.0)/sample_size) * ( (V*sample_size/(sample_size +1)) + ((mu-d)*(mu-d))/sample_size ) : V;
                    sigma = std::sqrt( V );
                    sample_size = sample_size + 1.0;
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
            get_alignment(aln_params, n1_max, ref_len_map, ref_seqs, read1, read1_rc, read_len1, a1_indv_max, k, cnt1, rc_already_comp1,
                          did_not_fit, tot_ksw_aligned);
            is_aligned1[n1_max.nam_id] = a1_indv_max;
            tot_all_tried ++;
            alignment a2_indv_max;
//            a2_indv_max.sw_score = -10000;
            auto n2_max = all_nams2[0];
            get_alignment(aln_params, n2_max, ref_len_map, ref_seqs, read2, read2_rc, read_len2, a2_indv_max, k, cnt2, rc_already_comp2,
                          did_not_fit, tot_ksw_aligned);
            is_aligned2[n2_max.nam_id] = a2_indv_max;
            tot_all_tried ++;

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
                        get_alignment(aln_params, n1, ref_len_map, ref_seqs, read1, read1_rc, read_len1, a1, k, cnt1,
                                      rc_already_comp1,
                                      did_not_fit, tot_ksw_aligned);
                        is_aligned1[n1.nam_id] = a1;
                        tot_all_tried++;
                    }
                } else { //rescue
//                    std::cerr << "RESCUE HERE1" << std::endl;
                    //////// Force SW alignment to rescue mate /////////
//                    std::cerr << query_acc2 << " RESCUE MATE 1" << a1.is_rc << " " n1.is_rc << std::endl;
                    rescue_mate(aln_params, n2, ref_len_map, ref_seqs, read1, read1_rc, read_len1, a1, rc_already_comp1, tot_ksw_aligned, mu, sigma, tot_rescued, k);
//                    is_aligned1[n1.nam_id] = a1;
                    tot_all_tried ++;
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
                        get_alignment(aln_params, n2, ref_len_map, ref_seqs, read2, read2_rc, read_len2, a2, k, cnt2,
                                      rc_already_comp2,
                                      did_not_fit, tot_ksw_aligned);
                        is_aligned2[n2.nam_id] = a2;
                        tot_all_tried++;
                    }
                } else{
//                    std::cerr << "RESCUE HERE2" << std::endl;
                    //////// Force SW alignment to rescue mate /////////
//                    std::cerr << query_acc1 << " RESCUE MATE 2" << a1.is_rc << " " n1.is_rc << std::endl;
                    rescue_mate(aln_params, n1, ref_len_map, ref_seqs, read2, read2_rc, read_len2, a2, rc_already_comp2, tot_ksw_aligned, mu, sigma, tot_rescued, k);
//                    is_aligned2[n2.nam_id] = a2;
                    tot_all_tried ++;
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
                append_to_sam(sam_string, sam_aln1, sam_aln2, read1, read2, read1_rc, read2_rc, acc_map, query_acc1,
                              query_acc2, mapq1, mapq2, mu, sigma, read_len1, true, record1.qual, record2.qual);
            } else {
                int max_out = high_scores.size() < max_secondary ? high_scores.size() : max_secondary;
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
                        append_to_sam(sam_string, sam_aln1, sam_aln2, read1, read2, read1_rc, read2_rc, acc_map, query_acc1,
                                      query_acc2, mapq1, mapq2, mu, sigma, read_len1, is_primary, record1.qual, record2.qual);
                    } else{
                        break;
                    }

                    prev_start_m1 = sam_aln1.ref_start;
                    prev_start_m2 = sam_aln2.ref_start;
                    prev_ref_id_m1 = sam_aln1.ref_id;
                    prev_ref_id_m2 = sam_aln2.ref_id;
                }
            }
            //////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////


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
            get_alignment(aln_params, n, ref_len_map, ref_seqs, read1, read1_rc, read_len1, a1, k, cnt1, rc_already_comp1, did_not_fit, tot_ksw_aligned);
            aln_scores1.push_back(a1);
            //////////////////////////////////////////////////////////////////

            //////// Force SW alignment to rescue mate /////////
            alignment a2;
//            std::cerr << query_acc2 << " force rescue" << std::endl;
            rescue_mate(aln_params, n, ref_len_map, ref_seqs, read2, read2_rc, read_len2, a2, rc_already_comp2, tot_ksw_aligned, mu, sigma,tot_rescued, k);
            aln_scores2.push_back(a2);
            //////////////////////////////////////////////////////////////////

            cnt1 ++;
            tot_all_tried ++;
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
            append_to_sam(sam_string,sam_aln1, sam_aln2, read1, read2, read1_rc, read2_rc, acc_map, query_acc1, query_acc2, mapq1, mapq2, mu, sigma, read_len1, true, record1.qual, record2.qual);
        } else {
            int max_out = high_scores.size() < max_secondary ? high_scores.size() : max_secondary;
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
                    append_to_sam(sam_string, sam_aln1, sam_aln2, read1, read2, read1_rc, read2_rc, acc_map, query_acc1,
                                  query_acc2, mapq1, mapq2, mu, sigma, read_len1, is_primary, record1.qual, record2.qual);
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
            //////// the actual testing of base pair alignment part start /////////
            alignment a2;
            get_alignment(aln_params, n, ref_len_map, ref_seqs, read2, read2_rc, read_len2, a2, k, cnt2, rc_already_comp2, did_not_fit, tot_ksw_aligned);
            aln_scores2.push_back(a2);
            //////////////////////////////////////////////////////////////////

            //////// Force SW alignment to rescue mate /////////
            alignment a1;
            rescue_mate(aln_params, n, ref_len_map, ref_seqs, read1, read1_rc, read_len1, a1, rc_already_comp1, tot_ksw_aligned, mu, sigma, tot_rescued, k);
            aln_scores1.push_back(a1);
            //////////////////////////////////////////////////////////////////

            cnt2 ++;
            tot_all_tried ++;
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
            append_to_sam(sam_string, sam_aln1, sam_aln2, read1, read2, read1_rc, read2_rc, acc_map, query_acc1,
                          query_acc2, mapq1, mapq2, mu, sigma, read_len1, true, record1.qual, record2.qual);
        } else {
            int max_out = high_scores.size() < max_secondary ? high_scores.size() : max_secondary;
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
                    append_to_sam(sam_string, sam_aln1, sam_aln2, read1, read2, read1_rc, read2_rc, acc_map, query_acc1,
                                  query_acc2, mapq1, mapq2, mu, sigma, read_len1, is_primary, record1.qual, record2.qual);
                } else{
                    break;
                }
            }
        }

    }

}


static inline void get_best_map_location(std::vector<std::tuple<int,nam,nam>> joint_NAM_scores, std::vector<nam> &nams1, std::vector<nam> &nams2, float &mu, float &sigma, float &sample_size, float &V, float &SSE, nam &best_nam1,  nam &best_nam2 ) {
    robin_hood::unordered_set<int> added_n1;
    robin_hood::unordered_set<int> added_n2;
    get_best_scoring_NAM_locations(nams1, nams2, joint_NAM_scores, mu, sigma, added_n1, added_n2 );
    nam n1_joint_max, n2_joint_max, n1_indiv_max, n2_indiv_max;
    float score_joint = 0;
    float score_indiv = 0;
    best_nam1.ref_s = -1; //Unmapped until proven mapped
    best_nam2.ref_s = -1; //Unmapped until proven mapped
    if (joint_NAM_scores.size() > 0) {
        // get best joint score
        for (auto &t : joint_NAM_scores) { // already sorted by descending score
            auto score_ = std::get<0>(t);
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

        if ((sample_size < 400) && (score_joint > score_indiv) ){
            int d = n1_joint_max.ref_s > n2_joint_max.ref_s ? n1_joint_max.ref_s - n2_joint_max.ref_s : n2_joint_max.ref_s - n1_joint_max.ref_s;
//            std::cerr << "HERE " << d << " " << mu <<  " " << sigma << " "<< n1_joint_max.ref_s << " " <<  n2_joint_max.ref_s << " "<< n1_joint_max.score << " " <<  n2_joint_max.score << std::endl;
            if ( d < 2000){
                float e;
                e = d - mu;
                mu = mu + e/sample_size; // (1.0/(sample_size +1.0)) * (sample_size*mu + d);
                SSE = SSE + e*(d-mu);
                V = sample_size > 1 ? SSE/(sample_size -1.0) : SSE; //d < 1000 ? ((sample_size +1.0)/sample_size) * ( (V*sample_size/(sample_size +1)) + ((mu-d)*(mu-d))/sample_size ) : V;
                sigma = std::sqrt( V );
                sample_size = sample_size + 1.0;
            }
        }
    }

}



void print_usage() {
    std::cerr << "\n";
    std::cerr << "StrobeAlign VERSION 0.5 \n";
    std::cerr << "\n";
    std::cerr << "StrobeAlign [options] <ref.fa> <reads1.fast[a/q.gz]> [reads2.fast[a/q.gz]]\n";
    std::cerr << "options:\n";
    std::cerr << "\n";
    std::cerr << "Resources:\n";
    std::cerr << "\t-t INT number of threads [3]\n";

    std::cerr << "\n";
    std::cerr << "Input/output:\n";
    std::cerr << "\t-o STR redirect output to file [stdout]\n";
    std::cerr << "\t-x Only map reads, no base level alignment (produces paf file)\n";
    std::cerr << "\t-N INT retain at most INT secondary alignments (is upper bounded by -M, and depends on -S) [0]  \n";
    std::cerr << "\t-L STR Print statistics of indexing to logfie [log.csv] \n";

    std::cerr << "\n";
    std::cerr << "Seeding:\n";
//    std::cerr << "\t-n INT number of strobes [2]\n";
    std::cerr << "\t-r INT Approximate read length. Sets suitable parameters for -k, -l, -u and -q. [150] \n";
    std::cerr << "\t-m INT Maximum seed length. Defaults to r - 50. For reasonable values on -l and -u, the seed length distribution is usually determined by\n\t   parameters l and u. Then, this parameter is only active in regions where syncmers are very sparse.\n";
    std::cerr << "\t-k INT strobe length, has to be below 32. [20]\n";
    std::cerr << "\t-l INT Lower syncmer offset from k/(k-s+1). Start sample second syncmer k/(k-s+1) + l syncmers downstream [0]\n";
    std::cerr << "\t-u INT Upper syncmer offset from k/(k-s+1). End sample second syncmer k/(k-s+1) + u syncmers downstream [7]\n";
    std::cerr << "\t-c INT Bitcount length between 2 and 63. [8]\n";
    std::cerr << "\t-s INT Submer size used for creating syncmers [k-4]. Only even numbers on k-s allowed.\n\t   A value of s=k-4 roughly represents w=10 as minimizer window [k-4]. It is recommended not to change this parameter\n\t   unless you have a good understanding of syncmenrs as it will drastically change the memory usage and results with non default values. \n";


    std::cerr << "\n";
    std::cerr << "Alignment:\n";
//    std::cerr << "\t-n INT number of strobes [2]\n";
    std::cerr << "\t-A INT matching score [2]\n";
    std::cerr << "\t-B INT mismatch penalty [8]\n";
    std::cerr << "\t-O INT gap open penalty [12]\n";
    std::cerr << "\t-E INT gap extension penalty [1]\n";


    std::cerr << "\n";
    std::cerr << "Search parameters:\n";
    std::cerr << "\t-f FLOAT top fraction of repetitive strobemers to filter out from sampling [0.0002]\n";
    std::cerr << "\t-S FLOAT Try candidate sites with mapping score at least S of maximum mapping score [0.5]\n";
    std::cerr << "\t-M INT Maximum number of mapping sites to try [20]\n";
    std::cerr << "\t-R INT Rescue level. Perform additional search for reads with many repetitive seeds filtered out.\n\t   This search includes seeds of R*repetitive_seed_size_filter (default: R=2). Higher R than default makes StrobeAlign\n\t   significantly slower but more accurate. R <= 1 deactivates rescue and is the fastest. \n";
}


int main (int argc, char **argv)
{

    if (argc < 3) {
        print_usage();
        return 0;
    }

    // Default parameters
    std::string choice = "randstrobes";
    int A = 2;
    int B = 8;
    int O = 12;
    int E = 1;

    int max_secondary = 0;
    bool is_sam_out = true; // true = align, false=map, default is_sam_out is align
    int n_threads = 3;
    int n = 2;
    int k = 20;
    int s = k - 4;
    float f = 0.0002;
    int R = 2;
    float dropoff_threshold = 0.5;
    int maxTries = 20;
    int l = 0;
    int u = 7;
    int c = 8;
    int r = 150;
    int max_seed_len;
    int max_dist = r - 50 < 255 ? r-50 : 255;
    bool r_set = false;
    bool max_seed_len_set = false;
    bool write_to_stdout = true;
    std::string output_file_name;
    std::string logfile_name = "log.csv";
    bool s_set = false;
    bool index_log = false;
    int opn = 1;
    while (opn < argc) {
        bool flag = false;
        if (argv[opn][0] == '-') {
//            if (argv[opn][1] == 'n') {
//                n = std::stoi(argv[opn + 1]);
//                opn += 2;
//                flag = true;
//            } else
            if (argv[opn][1] == 't') {
                n_threads = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'k') {
                k = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'o') {
                output_file_name = argv[opn + 1];
                opn += 2;
                flag = true;
                write_to_stdout = false;
            } else if (argv[opn][1] == 'L') {
                logfile_name = argv[opn + 1];
                index_log = true;
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 's') {
                s = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
                s_set = true;
            } else if (argv[opn][1] == 'f') {
                f = std::stof(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'x') {
                is_sam_out = false;
                opn += 1;
                flag = true;
            } else if (argv[opn][1] == 'R') {
                R = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'l') {
                l = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'u') {
                u = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'c') {
                c = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'r') {
                r = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
                r_set = true;
            } else if (argv[opn][1] == 'm') {
                max_seed_len = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
                max_seed_len_set = true;
            } else if (argv[opn][1] == 'S') {
                dropoff_threshold = std::stof(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'M') {
                maxTries = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'A') {
                A = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'B') {
                B = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'O') {
                O = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'E') {
                E = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'N') {
                max_secondary = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            }

            else {
                print_usage();
            }
        }
        if (!flag)
            break;
    }

    if (r_set) {
        if (r <= 125) { // based on params for 100
            k = 20;
            l = -2;
            u = 2;
        } else if ((r > 125) && (r <= 175)) { // based on params for 150
            k = 20;
            l = 1;
            u = 7;
        } else if ((r > 175) && (r <= 275)) { // based on params for 200 and 250
            k = 20;
            l = 4;
            u = 13;
        } else { // based on params for 300
            k = 22;
            l = 2;
            u = 12;
        }
    }

    if (!max_seed_len_set){
        max_dist = r - 70 > k ? r - 70 : k;
    } else {
        max_dist = max_seed_len - k; //convert to distance in start positions
    }

    if ( (!s_set ) ){
        s = k - 4; // Update default s to k - 4 if user has not set s parameter
    }
    uint64_t q;
    if ( (c < 64) && (c > 0)){
        q = pow (2, c) - 1;
    } else{
        std::cerr << "Warning wrong value for parameter c, setting c=8" << std::endl;
        q = pow (2, 8) - 1;
    }
    omp_set_num_threads(n_threads); // set number of threads in "parallel" blocks
    int w_min = k/(k-s+1) + l > 1 ? k/(k-s+1) + l : 1;
    int w_max = k/(k-s+1) + u;
    int t_syncmer = (k-s)/2 + 1;

    alignment_params aln_params;
    aln_params.match = A;
    aln_params.mismatch = B;
    aln_params.gap_open = O;
    aln_params.gap_extend = E;

    std::cerr << "Using" << std::endl;
    std::cerr << "k: " << k << std::endl;
    std::cerr << "s: " << s << std::endl;
    std::cerr << "w_min: " << w_min << std::endl;
    std::cerr << "w_max: " << w_max << std::endl;
    std::cerr << "maximum seed length: " << max_dist +k << std::endl;
    std::cerr << "threads: " << n_threads << std::endl;
    std::cerr << "R: " << R << std::endl;
    std::cerr << "[w_min, w_max] under thinning w roughly corresponds to sampling from downstream read coordinates (expected values): [" << (k-s+1)*w_min << ", " << (k-s+1)*w_max << "]" << std::endl;
    std::cerr << "A: " << A << std::endl;
    std::cerr << "B: " << B << std::endl;
    std::cerr << "O: " << O << std::endl;
    std::cerr << "E: " << E << std::endl;

//    assert(k <= (w/2)*w_min && "k should be smaller than (w/2)*w_min to avoid creating short strobemers");
    assert(k > 7 && "You should really not use too small strobe size!");
    assert(k <= 32 && "k have to be smaller than 32!");
    assert( ( s <= k ) && " s have to be smaller or equal to k!");
    assert( ( (k-s) % 2 == 0) && " k - s have to be an even number to create canonical syncmers. Set s to e.g., k-2, k-4, k-6, k-8.");
    assert(max_dist <= 255 && " -m (maximum seed length have to be smaller than 255 + k in v0.4 and up. If you need longer seeds, use v0.3");

//    assert(n == 2 && "Currently only n=2 is implemented");
    // File name to reference
    std::string ref_filename = argv[opn];
//    opn++;
//    const char *reads_filename = argv[opn];
//    opn++;
//    const char *reads_filename_PE2 = argv[opn];
    opn++;
    const char *reads_filename1;
    const char *reads_filename2;
    bool is_SE = false;
    if (opn == argc - 1) {
        reads_filename1 = argv[opn];
        is_SE = true;
    } else if (opn == argc - 2) {
        reads_filename1 = argv[opn];
        opn++;
        reads_filename2 = argv[opn];
    } else {
        print_usage();
        return 0;
    }


    //////////// CREATE INDEX OF REF SEQUENCES /////////////////

    // Record index creation start time

    auto start = std::chrono::high_resolution_clock::now();
    auto start_read_refs = std::chrono::high_resolution_clock::now();
    std::vector<std::string> ref_seqs;
    std::vector<unsigned int> ref_lengths;
    uint64_t total_ref_seq_size;
    idx_to_acc acc_map;
    total_ref_seq_size = read_references(ref_seqs, ref_lengths, acc_map, ref_filename);
    auto finish_read_refs = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_read_refs = finish_read_refs - start_read_refs;
    std::cerr << "Time reading references: " << elapsed_read_refs.count() << " s\n" <<  std::endl;

    auto start_flat_vector = std::chrono::high_resolution_clock::now();

    mers_vector flat_vector;
    int approx_vec_size = total_ref_seq_size / (k-s+1);
    std::cerr << "ref vector approximate size: " << approx_vec_size << std::endl;
    flat_vector.reserve(approx_vec_size);
    unsigned int mer_cnt = 0;
    for(size_t i = 0; i < ref_seqs.size(); ++i)
    {
        mers_vector randstrobes2; // pos, chr_id, kmer hash value
//        std::cerr << i << " " << i_mod << std::endl;
        randstrobes2 = seq_to_randstrobes2(n, k, w_min, w_max, ref_seqs[i], i, s, t_syncmer, q, max_dist);
        for (auto &t : randstrobes2)
        {
            flat_vector.push_back(t);
        }
    }
    std::cerr << "Ref vector actual size: " << flat_vector.size() << std::endl;
    flat_vector.shrink_to_fit();

    auto finish_generating_seeds = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_generating_seeds = finish_generating_seeds - start_flat_vector;
    std::cerr << "Time generating seeds: " << elapsed_generating_seeds.count() << " s\n" <<  std::endl;


//    create vector of vectors here nr_threads
//    std::vector<std::vector<std::tuple<uint64_t, unsigned int, unsigned int, unsigned int>>> vector_per_ref_chr(n_threads);
//    for(size_t i = 0; i < ref_seqs.size(); ++i)
//    {
//        mers_vector randstrobes2; // pos, chr_id, kmer hash value
//        std::cerr << "Started thread: " << omp_get_thread_num() << " chr size: " << ref_lengths[i] << " acc map:" << acc_map[i] << std::endl;
//        randstrobes2 = seq_to_randstrobes2(n, k, w_min, w_max, ref_seqs[i], i, s, t);
//        for (auto &t : randstrobes2)
//        {
//            vector_per_ref_chr[omp_get_thread_num()].push_back(t);
//        }
//        std::cerr << "Completed thread: " << omp_get_thread_num() << " chr size: " << ref_lengths[i] << " acc map:" << acc_map[i] << std::endl;
//    }

    uint64_t unique_mers = 0;
    auto start_sorting = std::chrono::high_resolution_clock::now();
//    std::cerr << "Reserving flat vector size: " << approx_vec_size << std::endl;
//    all_mers_vector_tmp.reserve(approx_vec_size); // reserve size corresponding to sum of lengths of all sequences divided by expected sampling
    process_flat_vector(flat_vector, unique_mers);
    auto finish_sorting = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_sorting_seeds = finish_sorting - start_sorting;
    std::cerr << "Time sorting seeds: " << elapsed_sorting_seeds.count() << " s\n" <<  std::endl;
    std::cerr << "Unique strobemers: " << unique_mers  <<  std::endl;

    auto finish_flat_vector = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_flat_vector = finish_flat_vector - start_flat_vector;
    std::cerr << "Total time generating flat vector: " << elapsed_flat_vector.count() << " s\n" <<  std::endl;

    auto start_hash_index = std::chrono::high_resolution_clock::now();
    kmer_lookup mers_index; // k-mer -> (offset in flat_vector, occurence count )
    mers_index.reserve(unique_mers);
    unsigned int filter_cutoff;
    filter_cutoff = index_vector(flat_vector, mers_index, f); // construct index over flat array
    auto finish_hash_index = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_hash_index = finish_hash_index - start_hash_index;
    std::cerr << "Total time generating hash table index: " << elapsed_hash_index.count() << " s\n" <<  std::endl;

//    mers_vector_reduced all_mers_vector;
//    all_mers_vector = remove_kmer_hash_from_flat_vector(flat_vector);
    /* destroy vector */
//    flat_vector.clear();


////////////////////////////////////////////////////////////////////////


//    std::cerr << "Wrote index to disc" << std::endl;

//    std::chrono::milliseconds timespan(10000); // or whatever
//    std::this_thread::sleep_for(timespan);
    //////////////////////////////////////////////////////////////////////////

    // Record index creation end time
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cerr << "Total time indexing: " << elapsed.count() << " s\n" <<  std::endl;

    if (index_log){
        std::cerr << "Printing log stats" << std::endl;
        print_diagnostics(flat_vector, mers_index, logfile_name, k);
        std::cerr << "Finished printing log stats" << std::endl;

    }

//    std::chrono::milliseconds timespan2(1000000); // or whatever
//    std::this_thread::sleep_for(timespan2);

    ///////////////////////////// MAP ///////////////////////////////////////

    // Record matching time
    auto start_aln_part = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> tot_read_file;
    std::chrono::duration<double> tot_construct_strobemers;
    std::chrono::duration<double> tot_find_nams;
    std::chrono::duration<double> tot_time_rescue;
    std::chrono::duration<double> tot_find_nams_alt;
    std::chrono::duration<double> tot_sort_nams;
    std::chrono::duration<double> tot_extend;
    std::chrono::duration<double> tot_rc;
    std::chrono::duration<double> tot_write_file;

    unsigned int tot_ksw_aligned = 0;
    unsigned int tot_rescued = 0;
    unsigned int tot_all_tried = 0;
    unsigned int did_not_fit = 0;
    unsigned int tried_rescue = 0;
//    std::ifstream query_file(reads_filename);

    int rescue_cutoff = R < 100 ? R*filter_cutoff : 1000;
    std::cerr << "Using rescue cutoff: " << rescue_cutoff <<  std::endl;

    std::streambuf *buf;
    std::ofstream of;

    if(!write_to_stdout) {
        of.open(output_file_name);
        buf = of.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }

    std::ostream out(buf);
//    std::ofstream out;
//    out.open(output_file_name);

    std::stringstream sam_output;
    std::stringstream paf_output;

    if (is_sam_out) {
        for (auto &it : acc_map) {
            out << "@SQ\tSN:" << it.second << "\tLN:" << ref_lengths[it.first] << "\n";
        }
        out << "@PG\tID:strobealign\tPN:strobealign\tVN:0.5\tCL:strobealign\n";
    }

    if(is_SE) {
        std::cerr << "Running SE mode" <<  std::endl;
        //    KSeq record;
        gzFile fp = gzopen(reads_filename1, "r");
        auto ks = make_ikstream(fp, gzread);
        int n_q_chunk_size = 1000000;
        KSeq record;
        std::string seq, seq_rc;
        unsigned int q_id = 0;
        std::pair<float, int> info;
        mers_vector_read query_mers; // pos, chr_id, kmer hash value
        std::vector<nam> nams; // (r_id, r_pos_start, r_pos_end, q_pos_start, q_pos_end)
//        std::vector<alignment> alignments;
//        alignments.reserve(maxTries);
        robin_hood::unordered_map< unsigned int, std::vector<hit>> hits_per_ref;
        std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> hits_fw;
        std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> hits_rc;
        hits_per_ref.reserve(100);
        hits_fw.reserve(5000);
        hits_rc.reserve(5000);
//    mers_vector_read query_mers_rc; // pos, chr_id, kmer hash value

//    std::vector<std::stringstream> output_streams(n_threads);
        std::vector<std::string> output_streams(n_threads);
        for (int i = 0; i < n_threads; ++i) {
            output_streams[i].reserve((n_q_chunk_size / n_threads + 1) *
                                      450); // Reserve sufficient space for appending multiple SAM records (400 is an upper setimate on the number of characters for each sam record of a 200-300bp read)
        }

        while (ks) {

            auto read_start = std::chrono::high_resolution_clock::now();
            auto records = ks.read(n_q_chunk_size);  // read a chunk of 500000 records
            auto read_finish = std::chrono::high_resolution_clock::now();
            tot_read_file += read_finish - read_start;
//        std::chrono::duration<double> elapsed_read = read_finish - read_start;
//        std::cerr << "Total time reading from file: " << elapsed_read.count() << " s\n" <<  std::endl;

            int n_it = records.size();
            std::cerr << "Mapping chunk of " << n_it << " query sequences... " << std::endl;
            #pragma omp parallel for num_threads(n_threads) shared(aln_params, output_streams, out, q_id, tot_all_tried, did_not_fit, tot_ksw_aligned, tried_rescue) private(sam_output, paf_output, record, seq_rc, query_mers, nams, hits_per_ref, info)
            for (int i = 0; i < n_it; ++i) {
                auto record = records[i];
                // generate mers here
                auto strobe_start = std::chrono::high_resolution_clock::now();
                query_mers = seq_to_randstrobes2_read(n, k, w_min, w_max, record.seq, q_id, s, t_syncmer, q, max_dist);
                auto strobe_finish = std::chrono::high_resolution_clock::now();
                tot_construct_strobemers += strobe_finish - strobe_start;

//            // Find NAMs alternative function
//            auto nam_alt_start = std::chrono::high_resolution_clock::now();
//            std::vector<nam> nams_alt; // (r_id, r_pos_start, r_pos_end, q_pos_start, q_pos_end)
//            nams_alt = find_nams_alt(query_mers, all_mers_vector, mers_index, k, ref_seqs, record.seq, hit_upper_window_lim,
//                                       filter_cutoff);
//            auto nam_alt_finish = std::chrono::high_resolution_clock::now();
//            tot_find_nams_alt += nam_alt_finish - nam_alt_start;

                // Find NAMs
//                std::cerr << "mapping " << record.name << std::endl;
                auto nam_start = std::chrono::high_resolution_clock::now();
                info = find_nams(nams, hits_per_ref, query_mers, flat_vector, mers_index, k, ref_seqs, record.seq, filter_cutoff);
                hits_per_ref.clear();
                auto nam_finish = std::chrono::high_resolution_clock::now();
                tot_find_nams += nam_finish - nam_start;

                if (R > 1) {
                    auto rescue_start = std::chrono::high_resolution_clock::now();
                    if ((nams.size() == 0) || (info.first < 0.7)) {
                        tried_rescue += 1;
                        nams.clear();
//                    std::cerr << "Rescue is_sam_out: " << record.name <<  std::endl;
                        info = find_nams_rescue(hits_fw, hits_rc, nams, hits_per_ref, query_mers, flat_vector, mers_index, k, ref_seqs,
                                                record.seq, rescue_cutoff);
                        hits_per_ref.clear();
                        hits_fw.clear();
                        hits_rc.clear();
//                    std::cerr << "Found: " << nams.size() <<  std::endl;
                    }
                    auto rescue_finish = std::chrono::high_resolution_clock::now();
                    tot_time_rescue += rescue_finish - rescue_start;
                }


                //Sort hits on score
                auto nam_sort_start = std::chrono::high_resolution_clock::now();
                std::sort(nams.begin(), nams.end(), score);
                auto nam_sort_finish = std::chrono::high_resolution_clock::now();
                tot_sort_nams += nam_sort_finish - nam_sort_start;

                auto extend_start = std::chrono::high_resolution_clock::now();
                if (!is_sam_out) {
                    output_hits_paf(output_streams[omp_get_thread_num()], nams, record.name, acc_map, k,
                                    record.seq.length(), ref_lengths);
//                output_streams[omp_get_thread_num()].append(paf_output.str()); // << paf_output.str();
                } else {
//                auto rc_start = std::chrono::high_resolution_clock::now();
//                auto rc_finish = std::chrono::high_resolution_clock::now();
//                tot_rc += rc_finish - rc_start;
//                    align_SE_multimap_deprecated(aln_params, output_streams[omp_get_thread_num()], nams, record.name, acc_map, k, record.seq.length(),
//                             ref_lengths, ref_seqs, record.seq,
//                             tot_ksw_aligned, tot_all_tried, dropoff_threshold, did_not_fit, maxTries, max_secondary);
                    if (max_secondary > 0){
                        // I created an entire new function here, dupliocating a lot of the code as outputting secondary hits is has some overhead to the
                        // original align_SE function (storing a vector of hits and sorting them)
                        // Such overhead is not present in align_PE - which implements both options in the same function.
                        align_SE_secondary_hits(aln_params,output_streams[omp_get_thread_num()], nams, record.name, acc_map, k, record.seq.length(),
                                 ref_lengths, ref_seqs, record.seq, record.qual,
                                 tot_ksw_aligned, tot_all_tried, dropoff_threshold, did_not_fit, maxTries, max_secondary);
                    } else {
                        align_SE(aln_params,output_streams[omp_get_thread_num()], nams, record.name, acc_map, k, record.seq.length(),
                                 ref_lengths, ref_seqs, record.seq, record.qual,
                                 tot_ksw_aligned, tot_all_tried, dropoff_threshold, did_not_fit, maxTries);
                    }

//                output_streams[omp_get_thread_num()] << sam_output.str();
                }
                auto extend_finish = std::chrono::high_resolution_clock::now();
                tot_extend += extend_finish - extend_start;
                q_id++;
                nams.clear();
            }
            // Output results
            auto write_start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < n_threads; ++i) {
                out << output_streams[i];
                output_streams[i].clear();
            }
            auto write_finish = std::chrono::high_resolution_clock::now();
            tot_write_file += write_finish - write_start;
        }
        gzclose(fp);
    }
    else{
        std::cerr << "Running PE mode" <<  std::endl;
        //    KSeq record;
//        int max_nam_n_hits1,max_nam_n_hits2;
        std::pair<float, int> info1, info2;
//        std::vector<int> isizes;
        gzFile fp1 = gzopen(reads_filename1, "r");
        auto ks1 = make_ikstream(fp1, gzread);
        gzFile fp2 = gzopen(reads_filename2, "r");
        auto ks2 = make_ikstream(fp2, gzread);
        int n_q_chunk_size = 1000000;
        KSeq record1;
        KSeq record2;
        std::string seq1, seq2, seq_rc1, seq_rc2;
        unsigned int q_id = 0;
        mers_vector_read query_mers1, query_mers2; // pos, chr_id, kmer hash value
        std::vector<nam> nams1;
        std::vector<nam> nams2;
        robin_hood::unordered_map< unsigned int, std::vector<hit>> hits_per_ref;
        std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> hits_fw;
        std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> hits_rc;
        hits_per_ref.reserve(100);
        hits_fw.reserve(5000);
        hits_rc.reserve(5000);
        std::vector<std::tuple<int,nam,nam>> joint_NAM_scores; // (score, aln1, aln2)
        std::vector<std::string> output_streams(n_threads);
        for (int i = 0; i < n_threads; ++i) {
            output_streams[i].reserve((n_q_chunk_size / n_threads + 1) *
                                      450); // Reserve sufficient space for appending multiple SAM records (400 is an upper setimate on the number of characters for each sam record of a 200-300bp read)
        }

        while (ks1) {

            auto read_start = std::chrono::high_resolution_clock::now();
            auto records1 = ks1.read(n_q_chunk_size);  // read a chunk of 1000000 records
            auto records2 = ks2.read(n_q_chunk_size);  // read a chunk of 1000000 records
            auto read_finish = std::chrono::high_resolution_clock::now();
            tot_read_file += read_finish - read_start;
            float sample_size = 1;
            float mu = 300;
            float sigma = 100;
            float V = 10000;
            float SSE = 10000;
            int n_it = records1.size();
            std::cerr << "Mapping chunk of " << n_it << " query sequences... " << std::endl;
            #pragma omp parallel for num_threads(n_threads) shared(aln_params, output_streams, out, q_id, tot_all_tried, did_not_fit, tot_ksw_aligned, tried_rescue, sample_size, mu, sigma, V, SSE) private(sam_output, paf_output, record1, record2, seq_rc1, seq_rc2, query_mers1, query_mers2, nams1, nams2, hits_per_ref, joint_NAM_scores, info1, info2)
            for (int i = 0; i < n_it; ++i) {
                auto record1 = records1[i];
                auto record2 = records2[i];
                // generate mers here
                auto strobe_start = std::chrono::high_resolution_clock::now();
//                std::cerr << "Going in! " << std::endl;
                query_mers1 = seq_to_randstrobes2_read(n, k, w_min, w_max, record1.seq, q_id, s, t_syncmer, q, max_dist);
//                std::cerr << "Lolz1 " << std::endl;
                query_mers2 = seq_to_randstrobes2_read(n, k, w_min, w_max, record2.seq, q_id, s, t_syncmer, q, max_dist);
//                std::cerr << "Lolz2 " << std::endl;
                auto strobe_finish = std::chrono::high_resolution_clock::now();
                tot_construct_strobemers += strobe_finish - strobe_start;
//                std::cerr << record1.name << " " << query_mers1.size() << std::endl;
//                std::cerr << record2.name << " " << query_mers2.size() << std::endl;

                // Find NAMs
                auto nam_start = std::chrono::high_resolution_clock::now();
                info1 = find_nams(nams1, hits_per_ref, query_mers1, flat_vector, mers_index, k, ref_seqs, record1.seq, filter_cutoff);
                hits_per_ref.clear();
                info2 = find_nams(nams2, hits_per_ref, query_mers2, flat_vector, mers_index, k, ref_seqs, record2.seq, filter_cutoff);
                hits_per_ref.clear();
                auto nam_finish = std::chrono::high_resolution_clock::now();
                tot_find_nams += nam_finish - nam_start;

//                if ( ((info1.second + info2.second) <= 10) || (info1.first < 0.7) || (info2.first < 0.7) ){
//                    tried_rescue +=2;
////                    nams1.clear();
////                    nams2.clear();
////                    std::cerr << "Rescue is_sam_out joint!: " << std::endl;
//                    info1 = find_nams(nams1, hits_per_ref, query_mers1, flat_vector, mers_index, k, ref_seqs, record1.seq, hit_upper_window_lim, 50000);
//                    hits_per_ref.clear();
//                    info2 = find_nams(nams2, hits_per_ref, query_mers2, flat_vector, mers_index, k, ref_seqs, record2.seq, hit_upper_window_lim, 50000);
//                    hits_per_ref.clear();
//                }

                if (R > 1) {
                    auto rescue_start = std::chrono::high_resolution_clock::now();
                    if ((nams1.size() == 0) || (info1.first < 0.7)) {
                        tried_rescue += 1;
                        nams1.clear();
//                        std::cerr << "Rescue is_sam_out read 1: " << record1.name << info1.first <<  std::endl;
                        info1 = find_nams_rescue(hits_fw, hits_rc, nams1, hits_per_ref, query_mers1, flat_vector, mers_index, k, ref_seqs,
                                                 record1.seq, rescue_cutoff);
                        hits_per_ref.clear();
                        hits_fw.clear();
                        hits_rc.clear();
//                    std::cerr << "Found: " << nams.size() <<  std::endl;
                    }

                    if ((nams2.size() == 0) || (info2.first < 0.7)) {
                        tried_rescue += 1;
                        nams2.clear();
//                        std::cerr << "Rescue is_sam_out read 2: " << record2.name << info2.first <<  std::endl;
                        info2 = find_nams_rescue(hits_fw, hits_rc, nams2, hits_per_ref, query_mers2, flat_vector, mers_index, k, ref_seqs,
                                                 record2.seq, rescue_cutoff);
                        hits_per_ref.clear();
                        hits_fw.clear();
                        hits_rc.clear();
//                    std::cerr << "Found: " << nams.size() <<  std::endl;
                    }
                    auto rescue_finish = std::chrono::high_resolution_clock::now();
                    tot_time_rescue += rescue_finish - rescue_start;
                }



                //Sort hits based on start choordinate on query sequence
                auto nam_sort_start = std::chrono::high_resolution_clock::now();
                std::sort(nams1.begin(), nams1.end(), score);
                std::sort(nams2.begin(), nams2.end(), score);
                auto nam_sort_finish = std::chrono::high_resolution_clock::now();
                tot_sort_nams += nam_sort_finish - nam_sort_start;

//                std::cerr << record1.name << std::endl;
//                for (auto &n : nams1){
//                    std::cerr << "NAM ORG: " << n.ref_id << ": (" << n.score << ", " << n.n_hits << ", " << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e  << ")" << std::endl;
//                }
//                std::cerr << record2.name << std::endl;
//                for (auto &n : nams2){
//                    std::cerr << "NAM ORG: " << n.ref_id << ": (" << n.score << ", " << n.n_hits << ", " << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e  << ")" << std::endl;
//                }

                auto extend_start = std::chrono::high_resolution_clock::now();
                if (!is_sam_out) {
//                    output_hits_paf(output_streams[omp_get_thread_num()], nams1, record1.name, acc_map, k,
//                                    record1.seq.length(), ref_lengths);
                    nam nam_read1;
                    nam nam_read2;
                    std::vector<std::tuple<int,nam,nam>> joint_NAM_scores;
                    get_best_map_location(joint_NAM_scores, nams1, nams2, mu, sigma, sample_size, V, SSE, nam_read1, nam_read2);
                    output_hits_paf_PE(output_streams[omp_get_thread_num()], nam_read1,  record1.name, acc_map, k, record1.seq.length(), ref_lengths);
                    output_hits_paf_PE(output_streams[omp_get_thread_num()], nam_read2,  record2.name, acc_map, k, record2.seq.length(), ref_lengths);
                    joint_NAM_scores.clear();
                } else {
                    align_PE(aln_params, output_streams[omp_get_thread_num()], nams1, nams2, record1, record2, acc_map, k,
                             ref_lengths, ref_seqs, tot_ksw_aligned, tot_all_tried, tot_rescued, dropoff_threshold, did_not_fit, mu, sigma, sample_size, V, SSE, maxTries, max_secondary);
                }
                auto extend_finish = std::chrono::high_resolution_clock::now();
                tot_extend += extend_finish - extend_start;
                q_id++;
                nams1.clear();
                nams2.clear();
            }

            std::cerr << "Estimated diff in start coordinates b/t mates, (mean: " << mu << ", stddev: " << sigma << ") " << std::endl;
//            std::cerr << "Len: " << isizes.size() << std::endl;

            // Output results
            auto write_start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < n_threads; ++i) {
                out << output_streams[i];
                output_streams[i].clear();
            }
            auto write_finish = std::chrono::high_resolution_clock::now();
            tot_write_file += write_finish - write_start;
        }
        gzclose(fp1);
        gzclose(fp2);

    }


//    out.close();

    std::cerr << "Total mapping sites tried: " << tot_all_tried << std::endl;
    std::cerr << "Total calls to ssw: " << tot_ksw_aligned << std::endl;
    std::cerr << "Calls to ksw (rescue mode): " << tot_rescued << std::endl;
    std::cerr << "Did not fit strobe start site: " << did_not_fit  << std::endl;
    std::cerr << "Tried rescue: " << tried_rescue  << std::endl;
    // Record mapping end time
    auto finish_aln_part = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> tot_aln_part = finish_aln_part - start_aln_part;
    std::cerr << "Total time mapping: " << tot_aln_part.count() << " s." <<  std::endl;
    std::cerr << "Total time reading read-file(s): " << tot_read_file.count() << " s." <<  std::endl;
    std::cerr << "Total time creating strobemers: " << tot_construct_strobemers.count()/n_threads << " s." <<  std::endl;
    std::cerr << "Total time finding NAMs (non-rescue mode): " << tot_find_nams.count()/n_threads  << " s." <<  std::endl;
    std::cerr << "Total time finding NAMs (rescue mode): " << tot_time_rescue.count()/n_threads  << " s." <<  std::endl;
//    std::cerr << "Total time finding NAMs ALTERNATIVE (candidate sites): " << tot_find_nams_alt.count()/n_threads  << " s." <<  std::endl;
    std::cerr << "Total time sorting NAMs (candidate sites): " << tot_sort_nams.count()/n_threads  << " s." <<  std::endl;
    std::cerr << "Total time reverse compl seq: " << tot_rc.count()/n_threads  << " s." <<  std::endl;
    std::cerr << "Total time base level alignment (ssw): " << tot_extend.count()/n_threads  << " s." <<  std::endl;
    std::cerr << "Total time writing alignment to files: " << tot_write_file.count() << " s." <<  std::endl;

    //////////////////////////////////////////////////////////////////////////


    /////////////////////// FIND AND OUTPUT NAMs ///////////////////////////////





}

