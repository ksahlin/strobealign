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

#include "source/kseq++.hpp"
using namespace klibpp;
#include "source/robin_hood.h"
#include "source/index.hpp"
//#include "gap_affine/affine_wavefront_align.h"
#include "source/ksw2.h"

//develop
#include <chrono>
#include <thread>
#include <sstream>


typedef robin_hood::unordered_map< unsigned int, std::string > idx_to_acc;

//typedef robin_hood::unordered_map< uint64_t, std::tuple<uint64_t, unsigned int >> vector_index;




static uint64_t read_references(std::vector<std::string> &seqs, std::vector<unsigned int> &lengths, idx_to_acc &acc_map, std::string fn)
{
    uint64_t total_ref_seq_size = 0;
    std::ifstream file(fn);
    std::string line, seq;
    int ref_index = 0;
    while (getline(file, line)) {
        if (line[0] == '>') {
//            std::cout << ref_index << " " << line << std::endl;
            if (seq.length() > 0){
//                seqs[ref_index -1] = seq;
                seqs.push_back(seq);
                lengths.push_back(seq.length());
                total_ref_seq_size += seq.length();
//                std::cout << ref_index - 1 << " here " << seq << " " << seq.length() << " " << seq.size() << std::endl;
//                generate_kmers(h, k, seq, ref_index);
            }
            acc_map[ref_index] = line.substr(1, line.length() -1); //line;
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
//        std::cout << ref_index -1 << " here2 " << seq << std::endl;
//        generate_kmers(h, k, seq, ref_index);
    }
    file.close();
    return total_ref_seq_size;
}



static inline void print_diagnostics_new4(mers_vector_reduced &mers_vector, kmer_lookup mers_index ) {
    uint64_t tot_flat_vector_size = 0;
    for (size_t i = 0; i < mers_vector.size(); ++i)
    {
        // access using []
        auto t = mers_vector[i];
//        std::cout << "(" << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << ", " << std::get<3>(t) << ", " << std::get<4>(t) << "), ";
        tot_flat_vector_size += sizeof(t);
    }
    std::cout << "Total size of flat mers-vector : " << tot_flat_vector_size/1000000  << " Mb." << std::endl;
    std::cout << "Total entries in flat mers-vector : " << mers_vector.size()  << std::endl;
//    uint64_t tot_hashtable_index_size = 0;
//    for (auto &it : mers_index)
//    {
////        std::cout << it.first << ": (" << std::get<0>(it.second) << ", " << std::get<1>(it.second) << "), " ;
//        tot_hashtable_index_size += sizeof(it.first);
//        tot_hashtable_index_size += sizeof(it.second);
//    }
//    std::cout << "Total size of hash table index : " << tot_hashtable_index_size/1000000  << " Mb." << std::endl;

    std::cout << "Total entries in hash table : " << mers_index.size()  << std::endl;
    std::cout << "Total size of hash table entries (ignoring internal nodes) : " << (mers_index.size() * sizeof(kmer_lookup::value_type))/1000000 << " Mb." << "\n";
    // https://stackoverflow.com/questions/720507/how-can-i-estimate-memory-usage-of-stdmap/720520
    std::cout << "Total size of hash table (applying approximation that about 2/3 of the memory is from internal nodes and 1/3 from the leafs, i.e., hash entries) : " << 3*(mers_index.size() * sizeof(kmer_lookup::value_type))/1000000 << " Mb." << "\n";
    std::cout << "" << std::endl;
    std::cout << "Total index size: " <<  tot_flat_vector_size/1000000 +  3*(mers_index.size() * sizeof(kmer_lookup::value_type))/1000000 << " Mb." << std::endl;

}


//static inline bool sort_hits(const hit &a, const hit &b)
//{
//    // first sort on ref ID, then on reference starts, then on query starts
//    return (a.ref_id < b.ref_id) ||
//           ( (a.ref_id == b.ref_id) && (a.ref_s < b.ref_s) ) ||
//           ((a.ref_id == b.ref_id) && (a.ref_s < b.ref_s) && (a.query_s == b.query_s )) ;
//}
//static inline std::vector<nam> find_nams_alt(mers_vector_read &query_mers, mers_vector_reduced &ref_mers, kmer_lookup &mers_index, int k, std::vector<std::string> &ref_seqs, std::string &read, unsigned int hit_upper_window_lim, unsigned int filter_cutoff ) {
////    std::cout << "ENTER FIND NAMS " <<  std::endl;
//    std::vector<hit> all_hits;
//    for (auto &q : query_mers)
////    for (size_t i = 0; i < query_mers.size(); ++i)
//    {
////        std::cout << "Q " << h.query_s << " " << h.query_e << " read length:" << read_length << std::endl;
//        uint64_t mer_hashv = std::get<0>(q);
//        if (mers_index.find(mer_hashv) != mers_index.end()) { //  In  index
//            hit h;
//            h.query_s = std::get<2>(q);
//            h.query_e = std::get<3>(q) + k; // h.query_s + read_length/2;
//            h.is_rc = std::get<4>(q);
//            std::tuple<uint64_t, unsigned int> mer;
//            mer = mers_index[mer_hashv];
//            uint64_t offset = std::get<0>(mer);
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
//        o.query_last_hit_pos = h.query_s;
//        o.ref_last_hit_pos = h.ref_s;
//        o.n_hits = 1;
//        o.is_rc = h.is_rc;
//    }
//
//    hit h;
//    for(size_t i = 1; i < all_hits.size(); ++i) // all but first element
//    {
//        h = all_hits[i];
//        if ( (o.ref_id == h.ref_id) && ( o.is_rc == h.is_rc) && ( o.query_last_hit_pos < h.query_s) && (h.query_s <= o.query_e ) && ( o.ref_last_hit_pos < h.ref_s) && (h.ref_s <= o.ref_e)){
//            if ( (h.query_e > o.query_e) && (h.ref_e > o.ref_e) ) {
//                o.query_e = h.query_e;
//                o.ref_e = h.ref_e;
////                o.previous_query_start = h.query_s;
////                o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
//                o.query_last_hit_pos = h.query_s; // log the last strobemer hit in case of outputting paf
//                o.ref_last_hit_pos = h.ref_s; // log the last strobemer hit in case of outputting paf
//                o.n_hits ++;
//            }
//            else if ((h.query_e <= o.query_e) && (h.ref_e <= o.ref_e)) {
////                o.previous_query_start = h.query_s;
////                o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
//                o.query_last_hit_pos = h.query_s; // log the last strobemer hit in case of outputting paf
//                o.ref_last_hit_pos = h.ref_s; // log the last strobemer hit in case of outputting paf
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
//            o.query_last_hit_pos = h.query_s;
//            o.ref_last_hit_pos = h.ref_s;
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
////        std::cout << "NAM ALT: " << n.ref_id << ": (" << n.n_hits << ", " << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e  << ")" << std::endl;
////    }
////    std::cout << " " << std::endl;
//
//    return final_nams;
//}


static inline std::vector<nam> find_nams(mers_vector_read &query_mers, mers_vector_reduced &ref_mers, kmer_lookup &mers_index, int k, std::vector<std::string> &ref_seqs, std::string &read, unsigned int hit_upper_window_lim, unsigned int filter_cutoff ){
//    std::cout << "ENTER FIND NAMS " <<  std::endl;
    robin_hood::unordered_map< unsigned int, std::vector<hit>> hits_per_ref; // [ref_id] -> vector( struct hit)
//    int read_length = read.length();
//    std::cout << " "  <<  std::endl;
    for (auto &q : query_mers)
//    for (size_t i = 0; i < query_mers.size(); ++i)
    {
//        std::cout << "Q " << h.query_s << " " << h.query_e << " read length:" << read_length << std::endl;
        uint64_t mer_hashv = std::get<0>(q);
        if (mers_index.find(mer_hashv) != mers_index.end()){ //  In  index
            hit h;
            h.query_s = std::get<2>(q);
            h.query_e = std::get<3>(q) + k; // h.query_s + read_length/2;
            h.is_rc = std::get<4>(q);

            std::tuple<uint64_t, unsigned int> mer;
            mer = mers_index[mer_hashv];
            uint64_t offset = std::get<0>(mer);
            unsigned int count = std::get<1>(mer);
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
                for(size_t j = offset; j < offset+count; ++j)
                {
                    auto r = ref_mers[j];
                    unsigned int ref_id = std::get<0>(r);
                    unsigned int ref_s = std::get<1>(r);
                    unsigned int ref_e = std::get<2>(r) + k; //ref_s + read_length/2;

                    h.ref_s = ref_s;
                    h.ref_e = ref_e;
                    hits_per_ref[ref_id].push_back(h);
//                    h.hit_count = count;
//                    std::cout << "Found: " <<  h.query_s << " " << h.query_e << " ref: " <<  h.ref_s << " " << h.ref_e << " " << h.is_rc << std::endl;
//                    hit_count_all ++;

                }
            }

        }
    }

//    std::cout << "NUMBER OF HITS GENERATED: " << hit_count_all << std::endl;
    std::vector<nam> open_nams;
    std::vector<nam> final_nams; // [ref_id] -> vector(struct nam)

    for (auto &it : hits_per_ref)
    {
        unsigned int ref_id = it.first;
        std::vector<hit> hits = it.second;
        open_nams = std::vector<nam> (); // Initialize vector
        uint64_t prev_q_start = 0;
        for (auto &h : hits){
            bool is_added = false;
//            std::cout << "HIT " << h.is_rc << " " << h.query_s <<  ", " << h.query_e << ", " << h.ref_s <<  ", " << h.ref_e << std::endl;

            for (auto & o : open_nams) {

                // Extend NAM
                if ( ( o.is_rc == h.is_rc) && ( o.query_last_hit_pos < h.query_s) && (h.query_s <= o.query_e ) && ( o.ref_last_hit_pos < h.ref_s) && (h.ref_s <= o.ref_e) ){
                    if ( (h.query_e > o.query_e) && (h.ref_e > o.ref_e) ) {
                        o.query_e = h.query_e;
                        o.ref_e = h.ref_e;
//                        o.previous_query_start = h.query_s;
//                        o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
                        o.query_last_hit_pos = h.query_s; // log the last strobemer hit in case of outputting paf
                        o.ref_last_hit_pos = h.ref_s; // log the last strobemer hit in case of outputting paf
                        o.n_hits ++;
//                        o.score += (float)1/ (float)h.hit_count;
                        is_added = true;
                        break;
                    }
                    else if ((h.query_e <= o.query_e) && (h.ref_e <= o.ref_e)) {
//                        o.previous_query_start = h.query_s;
//                        o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
                        o.query_last_hit_pos = h.query_s; // log the last strobemer hit in case of outputting paf
                        o.ref_last_hit_pos = h.ref_s; // log the last strobemer hit in case of outputting paf
                        o.n_hits ++;
//                        o.score += (float)1/ (float)h.hit_count;
                        is_added = true;
                        break;
                    }

                }


            }

            // Add the hit to open matches
            if (not is_added){
                nam n;
                n.query_s = h.query_s;
                n.query_e = h.query_e;
                n.ref_s = h.ref_s;
                n.ref_e = h.ref_e;
                n.ref_id = ref_id;
//                n.previous_query_start = h.query_s;
//                n.previous_ref_start = h.ref_s;
                n.query_last_hit_pos = h.query_s;
                n.ref_last_hit_pos = h.ref_s;
                n.n_hits = 1;
                n.is_rc = h.is_rc;
//                n.score += (float)1/ (float)h.hit_count;
                open_nams.push_back(n);
            }


            // Only filter if we have advanced at least k nucleotides
            if (h.query_s > prev_q_start + k) {

                // Output all NAMs from open_matches to final_nams that the current hit have passed
                for (auto &n : open_nams) {
                    if (n.query_e < h.query_s) {
                        final_nams.push_back(n);
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
            final_nams.push_back(n);
        }
    }

//    for (auto &n : final_nams){
////        std::cout << n.ref_id << ": (" << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e << ")" << std::endl;
//        std::cout << "NAM ORG: " << n.ref_id << ": (" << n.n_hits << ", " << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e  << ")" << std::endl;
////        std::cout << n.ref_id << ": (" << n.n_hits << ", " << n.ref_e - n.ref_s  <<  ", " << n.query_e - n.query_s << ")" << std::endl;
////        std::cout << " " << std::endl;
//    }

//    std::cout << "DONE" << std::endl;

    return final_nams;
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

static inline bool score(const nam &a, const nam &b)
{
    return ( (a.n_hits * (a.query_e - a.query_s)) > (b.n_hits * (b.query_e - b.query_s)) );
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
//    paf_output << query_acc << "\t" << read_len <<  "\t" << n.query_s << "\t" << n.query_last_hit_pos + k << "\t" << o  <<  "\t" << acc_map[n.ref_id] << "\t" << ref_len_map[n.ref_id] << "\t" << n.ref_s << "\t" << n.ref_last_hit_pos + k << "\t" << n.n_hits << "\t" << n.ref_last_hit_pos + k - n.ref_s << "\t" << "-" << "\n";
    paf_output.append(query_acc);
    paf_output.append("\t");
    paf_output.append(std::to_string(read_len));
    paf_output.append("\t");
    paf_output.append(std::to_string(n.query_s));
    paf_output.append("\t");
    paf_output.append(std::to_string(n.query_last_hit_pos + k));
    paf_output.append("\t");
    paf_output.append(o);
    paf_output.append("\t");
    paf_output.append(acc_map[n.ref_id]);
    paf_output.append("\t");
    paf_output.append(std::to_string( ref_len_map[n.ref_id]));
    paf_output.append("\t");
    paf_output.append(std::to_string(n.ref_s));
    paf_output.append("\t");
    paf_output.append(std::to_string( n.ref_last_hit_pos + k));
    paf_output.append("\t");
    paf_output.append(std::to_string( n.n_hits));
    paf_output.append("\t");
    paf_output.append(std::to_string( n.ref_last_hit_pos + k - n.ref_s));
    paf_output.append("\t-\n");
}

static inline std::string reverse_complement(std::string &read) {
    auto read_rev = read;
    std::reverse(read_rev.begin(), read_rev.end()); // reverse
//    std::cout << read_rev << std::endl;
    for (size_t j = 0; j < read_rev.length(); ++j) { // complement
        if (read_rev[j] == 'A') read_rev[j] = 'T';
        else if (read_rev[j] == 'T') read_rev[j] = 'A';
        else if (read_rev[j] == 'C') read_rev[j] = 'G';
        else if (read_rev[j] == 'G') read_rev[j] = 'C';
    }
    return read_rev;
}


//static inline std::string reverse_complement(std::string &read)
//{
//    auto read_rev = read;
//    reverse(read_rev.begin(), read_rev.end()); // reverse
//    for (std::size_t i = 0; i < read_rev.length(); ++i){
//        switch (read_rev[i]){
//            case 'A':
//                read_rev[i] = 'T';
//                break;
//            case 'C':
//                read_rev[i] = 'G';
//                break;
//            case 'G':
//                read_rev[i] = 'C';
//                break;
//            case 'T':
//                read_rev[i] = 'A';
//                break;
//        }
//    }
//    return read_rev;
//}


inline aln_info ksw_align(const char *tseq, int tlen, const char *qseq, int qlen,
               int sc_mch, int sc_mis, int gapo, int gape, ksw_extz_t &ez) {
    int8_t a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
    int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a, b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
    const uint8_t *ts = reinterpret_cast<const uint8_t *>(tseq);
    const uint8_t *qs = reinterpret_cast<const uint8_t *>(qseq);
    memset(&ez, 0, sizeof(ksw_extz_t));
    ksw_extz2_sse(0, qlen, qs, tlen, ts, 5, mat, gapo, gape, -1, -1, 0, 0, &ez);

    aln_info aln;
//    std::string cigar_mod;
//    cigar_mod.reserve(5*ez.n_cigar);
    unsigned int tstart_offset = 0;

    std::stringstream cigar_string;
    int edit_distance = 0;
    unsigned ref_pos = 0, read_pos = 0;
    for (int i = 0; i < ez.n_cigar; i++) {
        int count = ez.cigar[i] >> 4;
        char op = "MID"[ez.cigar[i] & 0xf];
        if ( (i==0) && op == 'D'){
            ref_pos += count;
            tstart_offset = ref_pos;
//            std::cout << "First deletion " << i << " " << count << std::endl;
            continue;
        }
        if ( (i==ez.n_cigar-1) && op == 'D'){
            ref_pos += count;
//            std::cout << "Last deletion " << i << " " << count << std::endl;
            continue;
        }
        cigar_string << count << op;
        switch (op) {
            case 'M':
                for (int j = 0; j < count; j++, ref_pos++, read_pos++) {
                    if (tseq[ref_pos] != qseq[read_pos])
                        edit_distance++;
                }
                break;
            case 'D':edit_distance += count;
                ref_pos += count;
                break;
            case 'I':edit_distance += count;
                read_pos += count;
                break;
            default:assert(0);
        }
//        std::cout << "ED " << edit_distance << std::endl;
    }
    aln.ed = edit_distance;
    aln.ref_offset = tstart_offset;
    aln.cigar = cigar_string.str();
    free(ez.cigar); //free(ts); free(qs);
    return aln;
}

inline int HammingDistance(std::string One, std::string Two)
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


static inline void align_SE(std::string &sam_string, std::vector<nam> &all_nams, std::string &query_acc, idx_to_acc &acc_map, int k, int read_len, std::vector<unsigned int> &ref_len_map, std::vector<std::string> &ref_seqs, std::string &read, unsigned int &tot_ksw_aligned, unsigned int &tot_all_tried, float dropoff, unsigned int &did_not_fit ) {

//    std::stringstream sam_string;
//    std::cout << "" << std::endl;
//    std::cout << query_acc << std::endl;
    std::string read_rc;
    bool rc_already_comp = false;

    if (all_nams.size() == 0) {
//        std::cout << "LEWL!!! "  << std::endl;
        sam_string.append(query_acc);
        sam_string.append("\t4\t*\t0\t255\t*\t*\t0\t0\t");
        sam_string.append(read);
        sam_string.append("\t*\n");
        return;
//        sam_string << query_acc << "\t" << 4 << "\t" << "*" << "\t" << 0
//                    << "\t" << 255 << "\t" << "*" << "\t" << "*" << "\t"
//                    << 0 << "\t" << 0 << "\t" << read << "\t" << "*" << "\n";
//        return sam_string;
    }

//        std::cout << "HERE!!! " << all_nams.size()  << std::endl;
    // Output results

    int cnt = 0;
    float score_dropoff;
//    float hits_max = (float) all_nams[0].n_hits;
    nam n_max = all_nams[0];
    float s1 = (float) (n_max.n_hits * (n_max.query_e - n_max.query_s));
    int mapq = 60; // MAPQ = 40(1−s2/s1) ·min{1,|M|/10} · log s1
    if (all_nams.size() > 1) {
        nam n_second = all_nams[1];
        float s2 = (float) (n_second.n_hits * (n_second.query_e - n_second.query_s));
//        ref_start = ref_tmp_start > 0 ? ref_tmp_start : 0;
        float min_matches;
        min_matches  = (float)n_max.n_hits/10 > 1 ? (float)n_max.n_hits/10 : 1;
        mapq = 40*(1 - s2/s1)*min_matches*log(s1) < 60 ? 40*(1 - s2/s1)*min_matches*log(s1) : 60 ;
    }

    int best_align_dist = ~0U >> 1;
    int best_align_index = 0; // assume by default it is the nam with most hits and most similar span length
    bool aln_did_not_fit;
//    std::cout << "best_align_dist: " << best_align_dist << std::endl;
    final_alignment sam_aln;
    // Only output single best hit based on: Firstly: number of randstrobe-hits. Secondly the concordance the span of the hits between ref and query (more simmilar ranked higher)
//    std::cout << "" << std::endl;


//    int n_it =  all_nams.size();
//    for(int i = 0; i < n_it; ++i){
//        auto n = all_nams[i];
    for (auto &n : all_nams) {
        aln_did_not_fit = false;
        score_dropoff = (float) n.n_hits / n_max.n_hits;

        if ( (cnt >= 20) || best_align_dist == 0 || score_dropoff < dropoff){ // only consider top 20 hits as minimap2 and break if alignment is exact match to reference or the match below droppoff cutoff.
            break;
        }

        tot_all_tried ++;

        unsigned int ref_diff = n.ref_e - n.ref_s;
        unsigned int read_diff = n.query_e - n.query_s;
        unsigned int min_diff =  read_diff ^ ((ref_diff ^ read_diff) & -(ref_diff < read_diff));
        unsigned int max_diff = ref_diff ^ ((ref_diff ^ read_diff) & -(ref_diff < read_diff));
        unsigned int diff = max_diff - min_diff;

        // deal with any read hanging of ends of reference not to get 'std::out_of_range' what(): basic_string::substr
        int ref_tmp_start = n.ref_s - n.query_s;
        int ref_tmp_segm_size = read_len + diff;
        int ref_len = ref_len_map[n.ref_id];
        int ref_start = ref_tmp_start > 0 ? ref_tmp_start : 0;
        int ref_segm_size = ref_tmp_segm_size < ref_len - ref_start ? ref_tmp_segm_size : ref_len - 1 - ref_start;

        std::string ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_segm_size);

        // decide if read should be fw or rc aligned to reference here by checking exact match of first and last strobe in the NAM

        if ( (ref_segm.substr(n.query_s, k) == read.substr(n.query_s, k) ) ) { //&& (ref_segm.substr(n.query_e - k + (ref_diff - read_diff), k) == read.substr(n.query_e - k, k)) ){
            n.is_rc = false;
        }
        else {
            if (!rc_already_comp){
                read_rc = reverse_complement(read);
                rc_already_comp = true;
            }

            if ((ref_segm.substr(n.query_s, k) == read_rc.substr(n.query_s, k))) { // && (ref_segm.substr(n.query_e - k + (ref_diff - read_diff), k) == read_rc.substr(n.query_e - k, k)) ){
                n.is_rc = true;
            } else {
                did_not_fit++;
                aln_did_not_fit = true;
            }
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

        if (ref_segm.length() >= read_len){
            hamming_dist = HammingDistance(r_tmp, ref_segm.substr(0,read_len));
        }

        if ( (hamming_dist) >=0 && (diff == 0)) { //Substitutions only
            if (hamming_dist < best_align_dist){
                best_align_index = cnt;
                best_align_dist = hamming_dist;
                sam_aln.cigar = std::to_string(read_len) + "M";
                sam_aln.ed = hamming_dist;
                sam_aln.ref_start = n.ref_s - n.query_s +1; // +1 because SAM is 1-based!
                sam_aln.is_rc = is_rc;
                sam_aln.ref_id = n.ref_id;
            }
        } else if ( (best_align_dist > 1) || did_not_fit ){
            int a = n.ref_s - n.query_s;
            int ref_start = std::max(0, a);
            int b = n.ref_e + (read_len - n.query_e);
            int ref_len = ref_seqs[n.ref_id].size();
            int ref_end = std::min(ref_len, b);
            std::string ref_segm = ref_seqs[n.ref_id].substr(ref_start, ref_end - ref_start);
            ksw_extz_t ez;
            const char *ref_ptr = ref_segm.c_str();
            const char *read_ptr = r_tmp.c_str();
            aln_info info;
            info = ksw_align(ref_ptr, ref_segm.size(), read_ptr, r_tmp.size(), 1, 4, 6, 1, ez);

            tot_ksw_aligned ++;
            if (info.ed < best_align_dist){
                best_align_index = cnt;
                best_align_dist = info.ed;
                sam_aln.cigar = info.cigar;
                sam_aln.ed = info.ed;
                sam_aln.ref_start =  a + info.ref_offset +1; // +1 because SAM is 1-based!
                sam_aln.is_rc = is_rc;
                sam_aln.ref_id = n.ref_id;
            }

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
        //TODO: Best way to calc Alignment score?
//        std::stringstream ss;
//        ss << query_acc << "\t" << o << "\t" << acc_map[sam_aln.ref_id] << "\t" << sam_aln.ref_start
//                    << "\t" << mapq << "\t" << sam_aln.cigar << "\t" << "*" << "\t"
//                    << 0 << "\t" << 0 << "\t" << output_read << "\t" << "*" << "\tNM:i:" << sam_aln.ed << "\n";
        sam_string.append(query_acc);
        sam_string.append("\t");
        sam_string.append(std::to_string(o));
        sam_string.append("\t");
        sam_string.append(acc_map[sam_aln.ref_id]);
        sam_string.append("\t");
        sam_string.append(std::to_string(sam_aln.ref_start));
        sam_string.append("\t");
        sam_string.append(std::to_string(mapq));
        sam_string.append("\t");
        sam_string.append(sam_aln.cigar);
        sam_string.append("\t*\t0\t0\t");
        sam_string.append(output_read);
        sam_string.append("\t*\tNM:i:");
        sam_string.append(std::to_string(sam_aln.ed));
        sam_string.append("\n");
//        return sam_string;

    }
//    return sam_string;
}


void print_usage() {
    std::cerr << "\n";
    std::cerr << "StrobeAlign VERSION 0.0.3\n";
    std::cerr << "\n";
    std::cerr << "StrobeAlign [options] <ref.fa> <reads1.fast[a/q.gz]> [reads2.fast[a/q.gz]]\n";
    std::cerr << "options:\n";
    std::cerr << "\t-t INT number of threads [3]\n";
    std::cerr << "\t-n INT number of strobes [2]\n";
    std::cerr << "\t-k INT strobe length [22]\n";
    std::cerr << "\t-o name of output SAM-file to print results to [mapped.sam]\n";
    std::cerr << "\t-x Only map reads, no base level alignment (produces paf file)\n";
    std::cerr << "\t-s INT syncmer thinning parameter to sample strobes. A value of s=k-4 roughly represents w=10 as minimizer window [k-4]. \n";
    std::cerr << "\t-f FLOAT top fraction of repetitive syncmers to filter out from sampling [0.0002]\n";
}


int main (int argc, char **argv)
{

    if (argc < 3) {
        print_usage();
        return 0;
    }

    // Default parameters
    std::string choice = "randstrobes";
    bool mode = true; // true = align, false=map, default mode is align

    int n_threads = 3;
    int n = 2;
    int k = 22;
    int s = k - 4;
    float f = 0.0002;
    std::string output_file_name = "mapped.sam";
    bool s_set = false;

    int opn = 1;
    while (opn < argc) {
        bool flag = false;
        if (argv[opn][0] == '-') {
            if (argv[opn][1] == 'n') {
                n = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 't') {
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
                mode = false;
                opn += 1;
                flag = true;
            }

            else {
                print_usage();
            }
        }
        if (!flag)
            break;
    }

    if ( (!s_set ) ){
        s = k - 4; // Update default s to k - 4 if user has not set s parameter
    }
    omp_set_num_threads(n_threads); // set number of threads in "parallel" blocks
    int w_min = k/(k-s+1)+2;
    int w_max = k/(k-s+1) + 10;
    int hit_upper_window_lim = (k-s+1)*w_max;
    float dropoff = 0.5;
    int t = (k-s)/2 + 1;
    std::cout << "Using" << std::endl;
    std::cout << "n: " << n << std::endl;
    std::cout << "k: " << k << std::endl;
    std::cout << "s: " << s << std::endl;
    std::cout << "t: " << t << std::endl;
    std::cout << "w_min: " << w_min << std::endl;
    std::cout << "w_max: " << w_max << std::endl;
    std::cout << "[w_min, w_max] under thinning w roughly corresponds to sampling from downstream read coordinates (under random minimizer sampling): [" << (k-s+1)*w_min << ", " << hit_upper_window_lim << "]" << std::endl;

//    assert(k <= (w/2)*w_min && "k should be smaller than (w/2)*w_min to avoid creating short strobemers");
    assert(k > 7 && "You should really not use too small strobe size!");
    assert(k <= 32 && "k have to be smaller than 32!");
    assert( ( s <= k ) && " s have to be smaller or equal to k!");
    assert( ( (k-s) % 2 == 0) && " k - s have to be an even number to create canonical syncmers. Set s to e.g., k-2, k-4, k-6, k-8.");
    assert(n == 2 && "Currently only n=2 is implemented");
    // File name to reference
    std::string ref_filename = argv[opn];
//    opn++;
//    const char *reads_filename = argv[opn];
//    opn++;
//    const char *reads_filename_PE2 = argv[opn];
    opn++;
    const char *reads_filename;
    const char *reads_filename_PE2;
    bool is_SE = false;
    if (opn == argc - 1) {
        reads_filename = argv[opn];
        is_SE = true;
    } else if (opn == argc - 2) {
        reads_filename = argv[opn];
        opn++;
        reads_filename_PE2 = argv[opn];
    } else {
        print_usage();
        return 0;
    }

    ///////////////////// INPUT /////////////////////////
//    std::string filename  = "test_ploy2.txt";
//    std::string reads_filename  = "test_ploy2.txt";

//    std::string filename  = "example_repeats.txt";
//    std::string reads_filename  = "example_repeats.txt";

//    std::string filename  = "example3.txt";
//    std::string reads_filename  = "example3.txt";

//    std::string filename  = "ecoli_repeats.txt";
//    std::string reads_filename  = "ecoli_repeats.txt";

//    std::string filename  = "ecoli_bug.txt";
//    std::string reads_filename  = "ecoli_bug.txt";

//    std::string filename  = "ecoli_randmer_bug.txt";
//    std::string reads_filename  = "ecoli_randmer_bug.txt";

//    std::string filename  = "/Users/kxs624/Documents/workspace/strobemers/cmake-build-debug/ecoli.fa";
//    std::string reads_filename  = "ecoli.fa";
//    std::string reads_filename  = "/Users/kxs624/Documents/workspace/strobemers/cmake-build-debug/SRR8187994_1_250k_subset.fasta";
//    std::string reads_filename  = "/Users/kxs624/Documents/workspace/strobemers/cmake-build-debug/SRR8187994_1_50_subset.fasta";
//
//    std::string filename  = "hg38_chr21.fa";
//    std::string reads_filename  = "hg38_chr21.fa";


//    std::string filename  = "/Users/kxs624/Documents/data/genomes/human/chm13_chr21.fa";
//    std::string reads_filename  = "/Users/kxs624/Documents/workspace/StrobeAlign/data/chm13_chr21_reads.fa";

//    std::string filename  = "/Users/kxs624/Documents/data/genomes/human/chm13_chr1.fa";
//    std::string reads_filename  = "/Users/kxs624/Documents/workspace/StrobeAlign/data/chm13_chr1_100k_reads.fa";


//    std::string filename  = "/Users/kxs624/Documents/data/genomes/human/hg38_chr21.fa";
//    std::string reads_filename  = "/Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr21_10k_reads_high_error.fa";
//    std::string filename  = "/Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr21_bug_ref.fa";
//    std::string reads_filename  = "/Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr21_1M_reads.fa";
//    std::string reads_filename  = "/Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr21_reads.fa";
//    std::string reads_filename  = "/Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr21_erroneous.fa";
//    std::string reads_filename  = "/Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr21_100k_reads.fa";
//    std::string reads_filename  = "/Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr21_bug2_reads.fa";


//    std::string filename  = "/Users/kxs624/Documents/data/genomes/human/hg38_chr1.fa";
//    std::string reads_filename  = "/Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr1_1M_reads.fa";
//    std::string output_filename  = "/Users/kxs624/Documents/workspace/StrobeAlign/data/tmp_out.sam";

//    std::string filename  = "hg21_bug.txt";
//    std::string reads_filename  = "hg21_bug.txt";

//    std::string filename  = "/Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr21_10k_unaligned_ref.fa";
//    std::string reads_filename  = "/Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr21_10k_unaligned_reads.fa";
//    std::string choice = "kmers";


    //////////////////////////////////////////////////////


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
    std::cout << "Time reading references: " << elapsed_read_refs.count() << " s\n" <<  std::endl;

    auto start_flat_vector = std::chrono::high_resolution_clock::now();
    pos_index tmp_index; // hash table holding all reference mers

    if (choice == "kmers" ){

        for(size_t i = 0; i < ref_seqs.size(); ++i)
        {
            mers_vector kmers; //  kmer hash value, pos, chr_id
//            kmers = seq_to_kmers(k, ref_seqs[i], i);
            tmp_index[i] = kmers;
        }
    }
    else if (choice == "randstrobes" ){
        if (n == 2 ){
            for(size_t i = 0; i < ref_seqs.size(); ++i)
            {
                mers_vector randstrobes2; // pos, chr_id, kmer hash value
                randstrobes2 = seq_to_randstrobes2(n, k, w_min, w_max, ref_seqs[i], i, s, t);
                tmp_index[i] = randstrobes2;
            }

        }
        else if (n == 3){
            for(size_t i = 0; i < ref_seqs.size(); ++i)
            {
                mers_vector randstrobes3; // pos, chr_id, kmer hash value
//                randstrobes3 = seq_to_randstrobes3(n, k, w_min, w_max, ref_seqs[i], i, w);
                tmp_index[i] = randstrobes3;
            }
        }
    }
    else {
        std::cout << choice << "not implemented : " << std::endl;
    }


    mers_vector all_mers_vector_tmp;
    uint64_t unique_mers = 0;
//    uint64_t approx_vec_size = total_ref_seq_size / (k-s+1);
//    std::cout << "Reserving flat vector size: " << approx_vec_size << std::endl;
//    all_mers_vector_tmp.reserve(approx_vec_size); // reserve size corresponding to sum of lengths of all sequences divided by expected sampling
    all_mers_vector_tmp = construct_flat_vector(tmp_index, unique_mers);
    std::cout << "Unique strobemers: " << unique_mers  <<  std::endl;


    auto finish_flat_vector = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_flat_vector = finish_flat_vector - start_flat_vector;
    std::cout << "Total time generating flat vector: " << elapsed_flat_vector.count() << " s\n" <<  std::endl;

    auto start_hash_index = std::chrono::high_resolution_clock::now();
    kmer_lookup mers_index; // k-mer -> (offset in flat_vector, occurence count )
    mers_index.reserve(unique_mers);
    unsigned int filter_cutoff;
    filter_cutoff = index_vector(all_mers_vector_tmp, mers_index, f); // construct index over flat array
    tmp_index.clear();
    mers_vector_reduced all_mers_vector;

    auto finish_hash_index = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_hash_index = finish_hash_index - start_hash_index;
    std::cout << "Total time generating hash table index: " << elapsed_hash_index.count() << " s\n" <<  std::endl;

    all_mers_vector = remove_kmer_hash_from_flat_vector(all_mers_vector_tmp);
    print_diagnostics_new4(all_mers_vector, mers_index);


////////////////////////////////////////////////////////////////////////


//    std::cout << "Wrote index to disc" << std::endl;

//    std::chrono::milliseconds timespan(10000); // or whatever
//    std::this_thread::sleep_for(timespan);
    //////////////////////////////////////////////////////////////////////////

    // Record index creation end time
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Total time indexing: " << elapsed.count() << " s\n" <<  std::endl;

//    std::chrono::milliseconds timespan2(1000000); // or whatever
//    std::this_thread::sleep_for(timespan2);

    ///////////////////////////// MAP ///////////////////////////////////////

    // Record matching time
    auto start_aln_part = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> tot_read_file;
    std::chrono::duration<double> tot_construct_strobemers;
    std::chrono::duration<double> tot_find_nams;
    std::chrono::duration<double> tot_find_nams_alt;
    std::chrono::duration<double> tot_sort_nams;
    std::chrono::duration<double> tot_extend;
    std::chrono::duration<double> tot_rc;
    std::chrono::duration<double> tot_write_file;

    unsigned int tot_ksw_aligned = 0;
    unsigned int tot_all_tried = 0;
    unsigned int did_not_fit = 0;

//    std::ifstream query_file(reads_filename);

    //    KSeq record;
    gzFile fp = gzopen(reads_filename, "r");
    auto ks = make_ikstream(fp, gzread);
    int n_q_chunk_size = 1000000;
    KSeq record;
    std::ofstream output_file;
    std::stringstream sam_output;
    std::stringstream paf_output;
    output_file.open(output_file_name);

    if (mode) {
        for (auto &it : acc_map) {
            output_file << "@SQ\tSN:" << it.second << "\tLN:" << ref_lengths[it.first] << "\n";
        }
        output_file << "@PG\tID:strobealign\tPN:strobealign\tVN:0.0.1\tCL:strobealign\n";
    }

    if(is_SE) {
        std::cout << "Woho SE mode! Let's go!" <<  std::endl;
        std::string line, seq, seq_rc, prev_acc;
        unsigned int q_id = 0;
        mers_vector_read query_mers; // pos, chr_id, kmer hash value
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
//        std::cout << "Total time reading from file: " << elapsed_read.count() << " s\n" <<  std::endl;

            int n_it = records.size();
            std::cout << "Mapping chunk of " << n_it << " query sequences... " << std::endl;
#pragma omp parallel for num_threads(n_threads) shared(output_streams, output_file, q_id, tot_all_tried, did_not_fit, tot_ksw_aligned) private(sam_output, paf_output, record, seq_rc, query_mers)
            for (int i = 0; i < n_it; ++i) {
                auto record = records[i];
                // generate mers here
                auto strobe_start = std::chrono::high_resolution_clock::now();
                query_mers = seq_to_randstrobes2_read(n, k, w_min, w_max, record.seq, q_id, s, t);
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
                auto nam_start = std::chrono::high_resolution_clock::now();
                std::vector<nam> nams; // (r_id, r_pos_start, r_pos_end, q_pos_start, q_pos_end)
                nams = find_nams(query_mers, all_mers_vector, mers_index, k, ref_seqs, record.seq, hit_upper_window_lim,
                                 filter_cutoff);
                auto nam_finish = std::chrono::high_resolution_clock::now();
                tot_find_nams += nam_finish - nam_start;

                //Sort hits based on start choordinate on query sequence
                auto nam_sort_start = std::chrono::high_resolution_clock::now();
                std::sort(nams.begin(), nams.end(), score);
                auto nam_sort_finish = std::chrono::high_resolution_clock::now();
                tot_sort_nams += nam_sort_finish - nam_sort_start;

                auto extend_start = std::chrono::high_resolution_clock::now();
                if (!mode) {
                    output_hits_paf(output_streams[omp_get_thread_num()], nams, record.name, acc_map, k,
                                    record.seq.length(), ref_lengths);
//                output_streams[omp_get_thread_num()].append(paf_output.str()); // << paf_output.str();
                } else {
//                auto rc_start = std::chrono::high_resolution_clock::now();
//                auto rc_finish = std::chrono::high_resolution_clock::now();
//                tot_rc += rc_finish - rc_start;
                    align_SE(output_streams[omp_get_thread_num()], nams, record.name, acc_map, k, record.seq.length(),
                             ref_lengths, ref_seqs, record.seq,
                             tot_ksw_aligned, tot_all_tried, dropoff, did_not_fit);
//                output_streams[omp_get_thread_num()] << sam_output.str();
                }
                auto extend_finish = std::chrono::high_resolution_clock::now();
                tot_extend += extend_finish - extend_start;
                q_id++;
            }
            // Output results
            auto write_start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < n_threads; ++i) {
                output_file << output_streams[i];
                output_streams[i].clear();
            }
            auto write_finish = std::chrono::high_resolution_clock::now();
            tot_write_file += write_finish - write_start;
        }
    }
    else{
        std::cout << "Woho PE mode! Let's go!" <<  std::endl;
    }


//    query_file.close();
    output_file.close();
    gzclose(fp);
    std::cout << "Total mapping sites tried: " << tot_all_tried << std::endl;
    std::cout << "Total calls to ksw: " << tot_ksw_aligned << std::endl;
    std::cout << "Did not fit strobe start site: " << did_not_fit  << std::endl;
    // Record mapping end time
    auto finish_aln_part = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> tot_aln_part = finish_aln_part - start_aln_part;
    std::cout << "Total time mapping: " << tot_aln_part.count() << " s." <<  std::endl;
    std::cout << "Total time reading read-file(s): " << tot_read_file.count() << " s." <<  std::endl;
    std::cout << "Total time creating strobemers: " << tot_construct_strobemers.count()/n_threads << " s." <<  std::endl;
    std::cout << "Total time finding NAMs (candidate sites): " << tot_find_nams.count()/n_threads  << " s." <<  std::endl;
//    std::cout << "Total time finding NAMs ALTERNATIVE (candidate sites): " << tot_find_nams_alt.count()/n_threads  << " s." <<  std::endl;
    std::cout << "Total time sorting NAMs (candidate sites): " << tot_sort_nams.count()/n_threads  << " s." <<  std::endl;
    std::cout << "Total time reverse compl seq: " << tot_rc.count()/n_threads  << " s." <<  std::endl;
    std::cout << "Total time extending alignment: " << tot_extend.count()/n_threads  << " s." <<  std::endl;
    std::cout << "Total time writing alignment to files: " << tot_write_file.count() << " s." <<  std::endl;

    //////////////////////////////////////////////////////////////////////////


    /////////////////////// FIND AND OUTPUT NAMs ///////////////////////////////





}

