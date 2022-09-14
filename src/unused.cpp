#include <string>
#include <vector>
#include <deque>
#include <sstream>
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

