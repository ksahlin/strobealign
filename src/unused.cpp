#include <string>
#include <vector>
#include <deque>
#include <sstream>
#include <cmath>
#include <bitset>
#include "xxhash.h"
#include "index.hpp"
#include "sam.hpp"
#include "aln.hpp"
#include "nam.hpp"

int main() {
}

struct hit {
    int query_s;
    int query_e;
    int ref_s;
    int ref_e;
    bool is_rc = false;
};


// defined in index.hpp
void update_window(std::deque <uint64_t> &q, std::deque <unsigned int> &q_pos, uint64_t &q_min_val, int &q_min_pos, uint64_t new_strobe_hashval, int i, bool &new_minimizer);
unsigned char seq_nt4_table[256];


/**********************************************************
 *
 * hash kmer into uint64
 *
 * *******************************************************/
// copy from http://www.cse.yorku.ca/~oz/hash.html:

uint64_t hash(const std::string& kmer)
{
    unsigned long hash = 5381;
    int c;
    for (std::string::size_type i=0; i< kmer.length(); i++) {
        c = kmer[i];
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
    }
    return hash;
}

/**********************************************************
 *
 * hash kmer into uint64
 *
 * *******************************************************/
// copy from minimap2:sketch.c :
static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}


static inline uint64_t sahlin_dna_hash(uint64_t key, uint64_t mask)
{
    key = (key << 3)|(key >> (64 - 3)); // rotate left with 11
    key = ~key; //flip
    key = (key << 13)|(key >> (64 - 13)); // rotate left with 13
    return key;
}

static inline uint64_t kmer_to_uint64(std::string &kmer, uint64_t kmask)
{
    uint64_t bkmer = 0;

    for (char i : kmer) {
        int c = seq_nt4_table[(uint8_t)i];
        bkmer = (bkmer << 2 | c) & kmask;

    }
    return bkmer;
}


//mers_vector seq_to_kmers(int k, std::string &seq, unsigned int ref_index)
//{
//    mers_vector kmers;
//    int l;
//    int i;
//    uint64_t mask=(1ULL<<2*k) - 1;
//    uint64_t x = 0;
//    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
//    int cnt = 0;
//    for (int i = l = 0; i <= seq.length()-1; i++) {
//        int c = seq_nt4_table[(uint8_t)seq[i]];
//        if (c < 4) { // not an "N" base
//            x = (x << 2 | c) & mask;                  // forward strand
//            if (++l >= k) { // we find a k-mer
//                uint64_t hash_k = x;
//                std::tuple<uint64_t, unsigned int, unsigned int, unsigned int> s (hash_k, ref_index, i-k+1, i-k+1);
//                kmers.push_back(s);
//                cnt ++;
//                if ((cnt % 1000000) == 0 ){
//                    std::cerr << cnt << " kmers created." << std::endl;
//                }
//            }
//        }
//        else {
//            l = 0, x = 0; // if there is an "N", restart
//        }
//
//    }
//    return  kmers;
//}


//mers_vector_reduced remove_kmer_hash_from_flat_vector(mers_vector &flat_vector){
//    mers_vector_reduced flat_vector_reduced;
//    for ( auto &t : flat_vector ) {
//        std::tuple<unsigned int, unsigned int,unsigned int> s( std::get<1>(t), std::get<2>(t), std::get<3>(t) );
//        flat_vector_reduced.push_back(s);
//    }
//    return flat_vector_reduced;
//}

//// initialize queue and current minimum and position
//static inline void initialize_window(std::vector<uint64_t> &string_hashes, std::deque <uint64_t> &q, uint64_t &q_min_val, int &q_min_pos, int w_min, int w_max, int k){
//    robin_hood::hash<std::string> robin_hash;
//    for (int i = w_min; i < w_max; i++) {
//        uint64_t hv = string_hashes[i];
//        q.push_back(hv);
//        if (hv < q_min_val) {
//            q_min_val = hv;
//            q_min_pos = i;
//        }
//    }
//}


//static inline void make_string_to_hashvalues_open_syncmers(std::string &seq, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, uint64_t kmask, int k, uint64_t smask, int s, int t) {
//    // initialize the deque
//    std::deque <uint64_t> qs;
//    std::deque <unsigned int> qs_pos;
//    int seq_length = seq.length();
//    int qs_size = 0;
//    uint64_t qs_min_val = UINT64_MAX;
//    int qs_min_pos = -1;
//
//
//    robin_hood::hash<uint64_t> robin_hash;
////    std::vector<std::tuple<uint64_t, unsigned int, unsigned int> > kmers;
//    unsigned int hash_count = 0;
//    int l;
////    int i;
//    uint64_t xk = 0;
//    uint64_t xs = 0;
//    for (int i = l = 0; i < seq_length; i++) {
//        int c = seq_nt4_table[(uint8_t) seq[i]];
//        if (c < 4) { // not an "N" base
//            xk = (xk << 2 | c) & kmask;
//            xs = (xs << 2 | c) & smask;
//            if (++l >= s) { // we find an s-mer
//                uint64_t hash_s = robin_hash(xs);
//                // que not initialized yet
//                if (qs_size < k - s + 1) {
//                    qs.push_back(hash_s);
//                    qs_pos.push_back(i - s + 1);
//                    qs_size++;
//                }
//                else{
//                    bool new_minimizer = false;
//                    update_window(qs, qs_pos, qs_min_val, qs_min_pos, hash_s, i - s + 1, new_minimizer );
//                    if (qs_min_pos == qs_pos[t-1]) { // occurs at t:th position in k-mer
//                        uint64_t hash_k = robin_hash(xk);
//                        string_hashes.push_back(hash_k);
//                        pos_to_seq_choord.push_back(i - k + 1);
//                        hash_count++;
//                    }
//                }
//            }
//        } else {
//            l = 0, xk = 0, xs = 0; // if there is an "N", restart
//            qs_size = 0;
//            qs_min_val = UINT64_MAX;
//            qs_min_pos = -1;
//            qs.clear();
//            qs_pos.clear();
//        }
//    }
////    std::cerr << hash_count << " values produced from string of length " <<   seq_length << std::endl;
////    for(auto t: pos_to_seq_choord){
////        std::cerr << t << " ";
////    }
////    std::cerr << " " << std::endl;
//}

static inline void make_string_to_hashvalues_closed_syncmers_canonical(std::string &seq, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, uint64_t kmask, int k, uint64_t smask, int s, int t, int startpos) {
    // initialize the deque
    std::deque <uint64_t> qs;
    std::deque <unsigned int> qs_pos;
    int seq_length = seq.length();
    int qs_size = 0;
    uint64_t qs_min_val = UINT64_MAX;
    int qs_min_pos = -1;


//    robin_hood::hash<uint64_t> robin_hash;
    uint64_t mask = (1ULL<<2*k) - 1;
//    std::cerr << mask << std::endl;

//    std::vector<std::tuple<uint64_t, unsigned int, unsigned int> > kmers;
    unsigned int hash_count = 0;
    int l;
    uint64_t xk[2];
    xk[0] = xk[1] = 0;
    uint64_t xs[2];
    xs[0] = xs[1] = 0;
    uint64_t kshift = (k - 1) * 2;
    uint64_t sshift = (s - 1) * 2;
    for (int i = l = 0; i < seq_length; i++) {
        int c = seq_nt4_table[(uint8_t) seq[i]];
        if (c < 4) { // not an "N" base
            xk[0] = (xk[0] << 2 | c) & kmask;                  // forward strand
            xk[1] = xk[1] >> 2 | (uint64_t)(3 - c) << kshift;  // reverse strand
            xs[0] = (xs[0] << 2 | c) & smask;                  // forward strand
            xs[1] = xs[1] >> 2 | (uint64_t)(3 - c) << sshift;  // reverse strand
            if (++l >= s) { // we find an s-mer
                uint64_t ys = std::min(xs[0], xs[1]);
//                uint64_t hash_s = robin_hash(ys);
                uint64_t hash_s = ys;
//                uint64_t hash_s = hash64(ys, mask);
//                uint64_t hash_s = XXH64(&ys, 8,0);
                // que not initialized yet
                if (qs_size < k - s ) {
                    qs.push_back(hash_s);
                    qs_pos.push_back(i - s + 1);
                    qs_size++;
                }
                else if (qs_size == k - s ) { // We are here adding the last s-mer and have filled queue up, need to decide for this k-mer (the first encountered) if we are adding it/
                    qs.push_back(hash_s);
                    qs_pos.push_back(i - s + 1);
                    qs_size++;
//                    std::cerr << qs_size << " "<< i - k + 1 << std::endl;
                    for (int j = 0; j < qs_size; j++) {
//                        std::cerr << qs_pos[j] << " " << qs[j] << " " << qs_min_val << std::endl;
                        if (qs[j] < qs_min_val) {
                            qs_min_val = qs[j];
                            qs_min_pos = qs_pos[j];
                        }
                    }
//                    if (qs_min_pos == qs_pos[t-1]) { // occurs at t:th position in k-mer
                    if ( (qs_min_pos == qs_pos[k - s]) || (qs_min_pos == qs_pos[0]) ) { // occurs at first or last position in k-mer
                        uint64_t yk = std::min(xk[0], xk[1]);
//                        uint64_t hash_k = robin_hash(yk);
//                        uint64_t hash_k = yk;
//                        uint64_t hash_k =  hash64(yk, mask);
                        uint64_t hash_k = XXH64(&yk, 8,0);
//                        uint64_t hash_k =  sahlin_dna_hash(yk, mask);
                        string_hashes.push_back(hash_k);
                        pos_to_seq_choord.push_back(startpos + i - k + 1);
                        hash_count++;
//                        std::cerr << i - s + 1 << " " << i - k + 1 << " " << (xk[0] < xk[1]) << std::endl;
//                        std::cerr <<  "Sampled gap " << std::endl;
                    }
                }
                else{
                    bool new_minimizer = false;
                    update_window(qs, qs_pos, qs_min_val, qs_min_pos, hash_s, i - s + 1, new_minimizer );
//                    if (qs_min_pos == qs_pos[t-1]) { // occurs at t:th position in k-mer
                    if ( (qs_min_pos == qs_pos[k - s]) || (qs_min_pos == qs_pos[0]) ) { // occurs at first or last position in k-mer


                        uint64_t yk = xk[0] < xk[1]? xk[0] : xk[1];
//                        uint64_t hash_k = robin_hash(yk);
//                        uint64_t hash_k = yk;
//                        uint64_t hash_k = hash64(yk, mask);
                        uint64_t hash_k = XXH64(&yk, 8, 0);
//                        uint64_t hash_k =  sahlin_dna_hash(yk, mask);
                        string_hashes.push_back(hash_k);
                        pos_to_seq_choord.push_back(startpos + i - k + 1);
//                        std::cerr << i - k + 1 << std::endl;
                        hash_count++;
//                        std::cerr << i - s + 1 << " " << i - k + 1 << " " << (xk[0] < xk[1]) << std::endl;
//                        std::cerr <<  "Sampled gap " << std::endl;
                    }
                }
            }
        } else {
            qs_min_val = UINT64_MAX;
            qs_min_pos = -1;
            l = 0, xs[0] = xs[1] = 0,  xk[0] = xk[1] = 0; // if there is an "N", restart
            qs_size = 0;
            qs.clear();
            qs_pos.clear();
        }
    }
//    std::cerr << hash_count << " values produced from string of length " <<   seq_length << std::endl;
//    for(auto t: pos_to_seq_choord){
//        std::cerr << t << " ";
//    }
//    std::cerr << " " << std::endl;
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


// from index.cpp
//
//static inline void make_string_to_hashvalues_random_minimizers(std::string &seq, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, int k, uint64_t kmask, int w) {
//    // initialize the deque
//    std::deque <uint64_t> q;
//    std::deque <unsigned int> q_pos;
//    int seq_length = seq.length();
//    int q_size = 0;
//    uint64_t q_min_val = UINT64_MAX;
//    int q_min_pos = -1;
//
//
//    robin_hood::hash<uint64_t> robin_hash;
////    std::vector<std::tuple<uint64_t, unsigned int, unsigned int> > kmers;
//    unsigned int hash_count = 0;
//    int l;
//    int i;
//    uint64_t x = 0;
//    for (int i = l = 0; i < seq_length; i++) {
//        int c = seq_nt4_table[(uint8_t) seq[i]];
//        if (c < 4) { // not an "N" base
//            x = (x << 2 | c) & kmask;                  // forward strand
//            if (++l >= k) { // we find a k-mer
//                uint64_t hash_k = robin_hash(x);
//
//                // que not initialized yet
//                if (q_size < w){
//                    q.push_back(hash_k);
//                    q_pos.push_back(i-k+1);
//                    q_size ++;
//                }
//                // que just filled up
//                else if (q_size == w - 1){
//                    q.push_back(hash_k);
//                    q_pos.push_back(i-k+1);
//                    q_min_val = UINT64_MAX;
//                    q_min_pos = -1;
//                    for (int j = 0; j <= w-1; j++) {
//                        if (q[j] < q_min_val) {
//                            q_min_val = q[j];
//                            q_min_pos = q_pos[j];
//                        }
//                    }
//                    string_hashes.push_back(q_min_val);
//                    pos_to_seq_choord.push_back(q_min_pos);
//                    hash_count ++;
//                }
//                // sliding the queue
//                else{
//                    bool new_minimizer = false;
//                    update_window(q, q_pos, q_min_val, q_min_pos, hash_k, i - k + 1, new_minimizer );
//                    if (new_minimizer) {
//                        string_hashes.push_back(q_min_val);
//                        pos_to_seq_choord.push_back(q_min_pos);
//                        hash_count++;
//                    }
//                }
//
//            }
//        } else {
//            l = 0, x = 0; // if there is an "N", restart
//        }
//    }
////    std::cerr << hash_count << " values produced from string of length " <<   seq_length << std::endl;
////    for(auto t: pos_to_seq_choord){
////        std::cerr << t << " ";
////    }
////    std::cerr << " " << std::endl;
//}

static inline void get_next_strobe(const std::vector<uint64_t> &string_hashes, uint64_t strobe_hashval, unsigned int &strobe_pos_next, uint64_t &strobe_hashval_next, unsigned int w_start, unsigned int w_end, uint64_t q){
    uint64_t min_val = UINT64_MAX;
//    int max_val = INT_MIN;
//    int min_val = INT_MAX;
//    int res;
//    std::bitset<64> b1,b2;
//        uint64_t rot_strobe_hashval = (strobe_hashval << q)|(strobe_hashval >> (64 - q));
//    uint64_t shift_strobe_hashval = strobe_hashval >> 5;
    std::bitset<64> b;
//    int a,b;
//    int p = pow (2, 4) - 1;
//    a = strobe_hashval & p;
//    int c = b1.count();
//    int d;

//    unsigned int min_pos;
//    min_pos = -1;
    for (auto i = w_start; i <= w_end; i++) {
//         Method 2
//        uint64_t res = (strobe_hashval + string_hashes[i]) & q;

//         Method 1
//        uint64_t res = (strobe_hashval + string_hashes[i]) % q;

        // Method 3 - seems to give the best tradeoff in speed and accuracy at this point
//        b = (strobe_hashval ^ string_hashes[i]);
//        uint64_t res = b.count();

        // Method 3' skew sample more for prob exact matching
        b = (strobe_hashval ^ string_hashes[i])  & q;
        uint64_t res = b.count();

        // Method by Lidon Gao (other strobemers library) and Giulio Ermanno Pibiri @giulio_pibiri
//        uint64_t res = (strobe_hashval ^ string_hashes[i]) ;

        // Method 6 Sahlin introduce skew (Method 3 and 3' are symmetrical for comp value of (s1,s2) and (s2,s1)
        // Methods 6 introduce asymmetry to reduce prob that we pick (s1,s2) and (s2,s1) as strobes to minimize fw and rc collisions
//        b = (shift_strobe_hashval ^ string_hashes[i])  & q;
//        uint64_t res = b.count();

        // Method 7 minimize collisions while still keeping small values space:
//        b = string_hashes[i] & p;
//        int res = a - b;

        if (res < min_val){
            min_val = res;
            strobe_pos_next = i;
//            std::cerr << strobe_pos_next << " " << min_val << std::endl;
            strobe_hashval_next = string_hashes[i];
        }
    }
//    std::cerr << "Offset: " <<  strobe_pos_next - w_start << " val: " << min_val <<  ", P exact:" <<  1.0 - pow ( (float) (8-min_val)/9, strobe_pos_next - w_start) << std::endl;

}


//mers_vector seq_to_randstrobes3(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index, int w)
//{
//    mers_vector randstrobes3;
//
//    if (seq.length() < 2*w_max) {
//        return randstrobes3;
//    }
//
//    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
//    uint64_t kmask=(1ULL<<2*k) - 1;
//    uint64_t q = pow (2, 16) - 1;
//    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
//    std::vector<uint64_t> string_hashes;
////    std::vector<uint64_t> pos_to_seq_choord;
//    std::vector<unsigned int> pos_to_seq_choord;
////    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
//    make_string_to_hashvalues_random_minimizers(seq, string_hashes, pos_to_seq_choord, k, kmask, w);
//    unsigned int seq_length = string_hashes.size();
//
////    std::cerr << seq << std::endl;
//
//    // create the randstrobes
//    for (unsigned int i = 0; i <= seq_length; i++) {
//
////        if ((i % 10) == 0 ){
////            std::cerr << i << " randstrobes created." << std::endl;
////        }
//        uint64_t strobe_hash;
//        strobe_hash = string_hashes[i];
//
//        unsigned int strobe_pos_next1;
//        uint64_t strobe_hashval_next1;
//        unsigned int strobe_pos_next2;
//        uint64_t strobe_hashval_next2;
//
//        if (i + 2*w_max < seq_length){
//            unsigned int w1_start = i+w_min;
//            unsigned int w1_end = i+w_max;
//            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next1, strobe_hashval_next1, w1_start, w1_end, q);
//
//            unsigned int w2_start = i+w_max + w_min;
//            unsigned int w2_end = i+2*w_max;
////            uint64_t conditional_next = strobe_hash ^ strobe_hashval_next1;
//            get_next_strobe(string_hashes, strobe_hashval_next1, strobe_pos_next2, strobe_hashval_next2, w2_start, w2_end, q);
//        }
//
//        else if ((i + 2*w_min + 1 < seq_length) && (seq_length <= i + 2*w_max) ){
//
//            int overshot;
//            overshot = i + 2*w_max - seq_length;
//            unsigned int w1_start = i+w_min;
//            unsigned int w1_end = i+w_max - overshot/2;
//            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next1, strobe_hashval_next1, w1_start, w1_end, q);
//
//            unsigned int w2_start = i+w_max - overshot/2 + w_min;
//            unsigned int w2_end = i+2*w_max - overshot;
////            uint64_t conditional_next = strobe_hash ^ strobe_hashval_next1;
//            get_next_strobe(string_hashes, strobe_hashval_next1, strobe_pos_next2, strobe_hashval_next2, w2_start, w2_end, q);
//        }
//        else{
//            return randstrobes3;
//        }
//
//        uint64_t hash_randstrobe3 = (strobe_hash/3) + (strobe_hashval_next1/4) + (strobe_hashval_next2/5);
//
//        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
////        unsigned int seq_pos_strobe2 =  pos_to_seq_choord[strobe_pos_next1]; //seq_pos_strobe1 + (strobe_pos_next1 - i); //
//        unsigned int seq_pos_strobe3 =  pos_to_seq_choord[strobe_pos_next2]; //seq_pos_strobe1 + (strobe_pos_next2 - i); //
////        std::cerr << i << " " << strobe_pos_next1 << " " << strobe_pos_next2 << " " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << seq_pos_strobe3 << " " << pos_to_seq_choord.size() << std::endl;
//
////         TODO: Take care of corner case (tmep if statement below. Some values in end of string produce a cororidnate of 0 for the last strobe. Probably an off-by-one error in the calculation of the strobe coord in the last strobe window
////        if (strobe_pos_next2 ==  seq_length){
//////            std::cerr << "OMGGGGGGG " << i << " " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << seq_pos_strobe3 << std::endl;
////            seq_pos_strobe3 = seq_length-1;
////        }
//        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int> s (hash_randstrobe3, ref_index, seq_pos_strobe1, seq_pos_strobe3);
//        randstrobes3.push_back(s);
//
//
////        auto strobe1 = seq.substr(i, k);
////        std::cerr << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next1 - (i+k), ' ') << std::string(k, 'X') << std::string(strobe_pos_next2 - strobe_pos_next1 - k, ' ') << std::string(k, 'X') << std::endl;
////        std::cerr << i << " " << strobe_pos_next1 << " " << strobe_pos_next2 << " " << seq_length << std::endl;
//    }
//    return randstrobes3;
//}


// from aln.cpp

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


// from aln.cpp
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


// from aln.cpp
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


// from aln.cpp


static inline bool sort_lowest_ed_scores_single(const std::tuple<int, alignment> &a,
                                                const std::tuple<int, alignment> &b)
{
    return (std::get<0>(a) < std::get<0>(b));
}


// from aln.cpp

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

// from aln.cpp

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


// index.cpp (StrobemerIndex::populate)
//    create vector of vectors here nr_threads
//    std::vector<std::vector<std::tuple<uint64_t, unsigned int, unsigned int, unsigned int>>> vector_per_ref_chr(opt.n_threads);
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


struct Hit {
    unsigned int count;
    unsigned int offset;
    unsigned int query_s;
    unsigned int query_e;
    bool is_rc;

    bool operator< (const Hit& rhs) const {
        return std::tie(count, offset, query_s, query_e, is_rc)
            < std::tie(rhs.count, rhs.offset, rhs.query_s, rhs.query_e, rhs.is_rc);
    }
};

static inline bool sort_hits(const hit &a, const hit &b)
{
    // first sort on query starts, then on reference starts
    return (a.query_s < b.query_s) || ( (a.query_s == b.query_s) && (a.ref_s < b.ref_s) );
}


static inline void find_nams_rescue(
    std::vector<Nam> &final_nams,
    robin_hood::unordered_map<unsigned int, std::vector<hit>> &hits_per_ref,
    const QueryRandstrobeVector &query_mers,
    const RefRandstrobeVector &ref_mers,
    RandstrobeMap &mers_index,
    int k,
    const std::vector<std::string> &ref_seqs,
    const std::string &read,
    unsigned int filter_cutoff
) {
    std::vector<Hit> hits_fw;
    std::vector<Hit> hits_rc;
    hits_fw.reserve(5000);
    hits_rc.reserve(5000);

//    std::pair<float,int> info (0,0); // (nr_nonrepetitive_hits/total_hits, max_nam_n_hits)
    int nr_good_hits = 0, total_hits = 0;
    bool is_rc = true, no_rep_fw = true, no_rep_rc = true;
//    std::pair<int, int> repeat_fw(0,0), repeat_rc(0,0);
//    std::vector<std::pair<int, int>> repetitive_fw, repetitive_rc;
    for (auto &q : query_mers)
    {
        auto mer_hashv = q.hash;
        if (mers_index.find(mer_hashv) != mers_index.end()){ //  In  index
            total_hits ++;
            auto ref_hit = mers_index[mer_hashv];
            auto offset = ref_hit.offset();
            auto count = ref_hit.count();
            auto query_s = q.start;
            auto query_e = q.end;
            is_rc = q.is_reverse;
            if (is_rc){
                Hit s{count, offset, query_s, query_e, is_rc};
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
                Hit s{count, offset, query_s, query_e, is_rc};
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
        auto count = q.count;
        auto offset = q.offset;
        h.query_s = q.query_s;
        h.query_e = q.query_e; // h.query_s + read_length/2;
        h.is_rc = q.is_rc;

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
                h.ref_s = r.position;
                h.ref_e = h.ref_s + r.strobe2_offset() + k;
//                h.count = count;
//                hits_per_ref[std::get<0>(r)].push_back(h);

                int diff = std::abs((h.query_e - h.query_s) - (h.ref_e - h.ref_s));
                if (diff <= min_diff ){
                    hits_per_ref[r.reference_index()].push_back(h);
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
        auto count = q.count;
        auto offset = q.offset;
        h.query_s = q.query_s;
        h.query_e = q.query_e; // h.query_s + read_length/2;
        h.is_rc = q.is_rc;

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
                h.ref_s = r.position;
                h.ref_e = h.ref_s + r.strobe2_offset() + k;
//                h.count = count;
//                hits_per_ref[std::get<1>(r)].push_back(h);
                int diff = std::abs((h.query_e - h.query_s) - (h.ref_e - h.ref_s));
                if (diff <= min_diff ){
                    hits_per_ref[r.reference_index()].push_back(h);
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
    std::vector<Nam> open_nams;
    int nam_id_cnt = 0;
//    std::vector<nam> final_nams; // [ref_id] -> vector(struct nam)

    for (auto &it : hits_per_ref)
    {
        auto ref_id = it.first;
        std::vector<hit> hits = it.second;
        std::sort(hits.begin(), hits.end(), sort_hits);
        open_nams = std::vector<Nam> (); // Initialize vector
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
                Nam n;
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
}



static inline std::pair<float,int> find_nams(
    std::vector<Nam> &final_nams,
    robin_hood::unordered_map<unsigned int, std::vector<hit>> &hits_per_ref,
    const QueryRandstrobeVector &query_mers,
    const RefRandstrobeVector &ref_mers,
    RandstrobeMap &mers_index,
    int k,
    unsigned int filter_cutoff
) {
//    robin_hood::unordered_map< unsigned int, std::vector<hit>> hits_per_ref; // [ref_id] -> vector( struct hit)
//    std::vector<std::vector<hit>> hits_per_ref(10);
//    int read_length = read.length();
    std::pair<float,int> info (0.0f,0); // (nr_nonrepetitive_hits/total_hits, max_nam_n_hits)
    int nr_good_hits = 0, total_hits = 0;
    hit h;
    for (auto &q : query_mers)
//    for (size_t i = 0; i < query_mers.size(); ++i)
    {
        auto mer_hashv = q.hash;
        if (mers_index.find(mer_hashv) != mers_index.end()) { //  In  index
            total_hits ++;
            h.query_s = q.start;
            h.query_e = q.end; // h.query_s + read_length/2;
            h.is_rc = q.is_reverse;
            auto mer = mers_index[mer_hashv];
            auto offset = mer.offset();
            auto count = mer.count();
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
                    h.ref_s = r.position;
                    h.ref_e = h.ref_s + r.strobe2_offset() + k;
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
                        hits_per_ref[r.reference_index()].push_back(h);
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
    std::vector<Nam> open_nams;
//    std::vector<nam> final_nams; // [ref_id] -> vector(struct nam)

    for (auto &it : hits_per_ref)
    {
        auto ref_id = it.first;
        std::vector<hit> hits = it.second;

//    for(size_t i = 0; i < hits_per_ref.size(); ++i){
//        unsigned int ref_id = i;
//        auto hits = hits_per_ref[i];
        open_nams = std::vector<Nam> (); // Initialize vector
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
                Nam n;
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
}
