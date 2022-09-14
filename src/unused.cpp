#include <string>
#include <vector>
#include <deque>
#include <sstream>
#include <bitset>
#include "xxhash.h"
#include "index.hpp"

int main() {
}


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
