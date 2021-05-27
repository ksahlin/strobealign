//
//  index.cpp
//  cpp_course
//
//  Created by Kristoffer Sahlin on 4/21/21.
//

#include "index.hpp"
#include <iostream>
#include <math.h>       /* pow */
#include <bitset>



/**********************************************************
 *
 * hash kmer into uint64
 *
 * *******************************************************/
// copy from http://www.cse.yorku.ca/~oz/hash.html:

uint64_t hash(std::string kmer)
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
}//hash64


static unsigned char seq_nt4_table[256] = {
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
}; //seq_nt4_table

static inline uint64_t kmer_to_uint64(std::string &kmer, uint64_t kmask)
{
    uint64_t bkmer = 0;

    for (char i : kmer) {
        int c = seq_nt4_table[(uint8_t)i];
        bkmer = (bkmer << 2 | c) & kmask;

    }
    return bkmer;
}


mers_vector construct_flat_vector_three_pos(one_pos_index &tmp_index){
    mers_vector flat_vector;
    for (auto &it : tmp_index)  {
        for (auto &t : it.second) // it.second is the vector of k-mers, t is a tuple
        {
            flat_vector.push_back(t);
        }
    }
    //    flat_array sort
    std::sort(flat_vector.begin(), flat_vector.end());
    return flat_vector;
}


kmer_lookup index_vector_one_pos(mers_vector  &flat_vector){
    kmer_lookup mers_index;
    uint64_t offset = 0;
    uint64_t prev_offset = 0;
    unsigned int count = 0;

    uint64_t prev_k;
    std::tuple<uint64_t, unsigned int, unsigned int> t = flat_vector[0];
    prev_k = std::get<0>(t);
    uint64_t curr_k;

    for ( auto &t : flat_vector ) {
//        std::cout << t << std::endl;
        curr_k = std::get<0>(t);
        if (curr_k == prev_k){
            count ++;
        }
        else {
            std::tuple<unsigned int, unsigned int> s(prev_offset, count);
            mers_index[prev_k] = s;
            count = 1;
            prev_k = curr_k;
            prev_offset = offset;
        }
        offset ++;
    }

    // last k-mer
    std::tuple<unsigned int, unsigned int> s(prev_offset, count);
    mers_index[curr_k] = s;

    return mers_index;
}

mers_vector_reduced remove_kmer_hash_from_flat_vector(mers_vector &flat_vector){
    mers_vector_reduced flat_vector_reduced;
    for ( auto &t : flat_vector ) {
        std::tuple<unsigned int, unsigned int> s( std::get<1>(t), std::get<2>(t) );
        flat_vector_reduced.push_back(s);
    }
    return flat_vector_reduced;
}

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

// update queue and current minimum and position
static inline void update_window(std::deque <uint64_t> &q, uint64_t &q_min_val, int &q_min_pos, uint64_t new_strobe_hashval, int w, int i, bool &new_minimizer){
    uint64_t popped_val;
    popped_val = q.front();
    q.pop_front();
    q.push_back(new_strobe_hashval);
    if (popped_val == q_min_val){ // we popped the minimum value, find new brute force
        q_min_val = UINT64_MAX;
        q_min_pos = i;
        for (int j = 0; j <= q.size()-1; j++) {
            if (q[j] < q_min_val) {
                q_min_val = q[j];
                q_min_pos = i + j;
                new_minimizer = true;
            }
        }
    }
    else if ( new_strobe_hashval < q_min_val ) { // the new value added to queue is the new minimum
        q_min_val = new_strobe_hashval;
        q_min_pos = i + w;
        new_minimizer = true;
    }
}



static inline void make_string_to_hashvalues2(std::string &seq, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, int k, uint64_t kmask, int w) {
    // initialize the deque
    std::deque <uint64_t> q;
    int seq_length = seq.length();
    int q_size = 0;
    uint64_t q_min_val = UINT64_MAX;
    int q_min_pos = -1;


    robin_hood::hash<uint64_t> robin_hash;
//    std::vector<std::tuple<uint64_t, unsigned int, unsigned int> > kmers;
    unsigned int hash_count = 0;
    int l;
    int i;
    uint64_t x = 0;
    for (int i = l = 0; i < seq_length; i++) {
        int c = seq_nt4_table[(uint8_t) seq[i]];
        if (c < 4) { // not an "N" base
            x = (x << 2 | c) & kmask;                  // forward strand
            if (++l >= k) { // we find a k-mer
                uint64_t hash_k = robin_hash(x);

                // que not initialized yet
                if (q_size < w){
                    q.push_back(hash_k);
                    q_size ++;
                }
                // que just filled up
                else if (q_size == w - 1){
                    q.push_back(hash_k);
                    q_min_val = UINT64_MAX;
                    q_min_pos = -1;
                    for (int j = 0; j <= w-1; j++) {
                        if (q[j] < q_min_val) {
                            q_min_val = q[j];
                            q_min_pos = i - k + 1 + j;
                        }
                    }
                    string_hashes.push_back(q_min_val);
                    pos_to_seq_choord.push_back(q_min_pos);
                    hash_count ++;
                }
                // sliding the queue
                else{
                    bool new_minimizer = false;
                    update_window(q, q_min_val, q_min_pos, hash_k, w, i - k + 1, new_minimizer );
                    if (new_minimizer) {
                        string_hashes.push_back(q_min_val);
                        pos_to_seq_choord.push_back(q_min_pos);
                        hash_count++;
                    }
                }

            }
        } else {
            l = 0, x = 0; // if there is an "N", restart
        }
    }
//    std::cout << hash_count << " values produced from string of length " <<   seq_length << std::endl;
//    for(auto t: pos_to_seq_choord){
//        std::cout << t << " ";
//    }
//    std::cout << " " << std::endl;
}



static inline void get_next_strobe(std::vector<uint64_t> &string_hashes, uint64_t strobe_hashval, unsigned int &strobe_pos_next, uint64_t &strobe_hashval_next,  unsigned int w_start, unsigned int w_end, uint64_t q){
    uint64_t min_val = UINT64_MAX;
//    unsigned int min_pos;
//    min_pos = -1;
    for (auto i = w_start; i <= w_end; i++) {
        uint64_t res = (strobe_hashval + string_hashes[i]) & q ;
        if (res < min_val){
            min_val = res;
            strobe_pos_next = i;
            strobe_hashval_next = string_hashes[i];
        }
    }
}

mers_vector seq_to_kmers(int k, std::string &seq, unsigned int ref_index)
{
    mers_vector kmers;
    int l;
    int i;
    uint64_t mask=(1ULL<<2*k) - 1;
    uint64_t x = 0;
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    int cnt = 0;
    for (int i = l = 0; i <= seq.length()-1; i++) {
        int c = seq_nt4_table[(uint8_t)seq[i]];
        if (c < 4) { // not an "N" base
            x = (x << 2 | c) & mask;                  // forward strand
            if (++l >= k) { // we find a k-mer
                uint64_t hash_k = x;
                std::tuple<uint64_t, unsigned int, unsigned int> s (hash_k, ref_index, i-k+1);
                kmers.push_back(s);
                cnt ++;
                if ((cnt % 1000000) == 0 ){
                    std::cout << cnt << " kmers created." << std::endl;
                }
            }
        }
        else {
            l = 0, x = 0; // if there is an "N", restart
        }

    }
    return  kmers;
}

mers_vector seq_to_randstrobes2(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index, int w)
{
    mers_vector randstrobes2;

    if (seq.length() < w_max) {
        return randstrobes2;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    uint64_t kmask=(1ULL<<2*k) - 1;
    uint64_t q = pow (2, 16) - 1;
//    std::bitset<64> x(q);
//    std::cout << x << '\n';
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask, w);
    unsigned int seq_length = string_hashes.size();

//    std::cout << seq << std::endl;

    // create the randstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//            std::cout << i << " strobemers created." << std::endl;
//        }
        unsigned int strobe_pos_next;
        uint64_t strobe_hashval_next;

        if (i + w_max <= seq_length){
            unsigned int w_start = i+w_min;
            unsigned int w_end = i+w_max;
            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
        }
        else if ((i + w_min + 1 < seq_length) && (seq_length < i + w_max) ){
            unsigned int w_start = i+w_min;
            unsigned int w_end = seq_length -1;
            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
        }
        else{
//            std::cout << randstrobes2.size() << " randstrobes generated" << '\n';
            return randstrobes2;
        }

//        uint64_t hash_randstrobe2 = (string_hashes[i] << k) ^ strobe_hashval_next;
        uint64_t hash_randstrobe2 = (string_hashes[i]/2) + (strobe_hashval_next/3);

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
//        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe_pos_next];
        std::tuple<uint64_t, unsigned int, unsigned int> s (hash_randstrobe2, ref_index, seq_pos_strobe1);
        randstrobes2.push_back(s);


//        auto strobe1 = seq.substr(i, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next - (i+k), ' ') << std::string(k, 'X') << std::endl;

    }
//    std::cout << randstrobes2.size() << " randstrobes generated" << '\n';
    return randstrobes2;
}


mers_vector seq_to_randstrobes3(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index, int w)
{
    mers_vector randstrobes3;

    if (seq.length() < 2*w_max) {
        return randstrobes3;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    uint64_t kmask=(1ULL<<2*k) - 1;
    uint64_t q = pow (2, 16) - 1;
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    std::vector<uint64_t> string_hashes;
//    std::vector<uint64_t> pos_to_seq_choord;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask, w);
    unsigned int seq_length = string_hashes.size();

//    std::cout << seq << std::endl;

    // create the randstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//            std::cout << i << " randstrobes created." << std::endl;
//        }
        uint64_t strobe_hash;
        strobe_hash = string_hashes[i];

        unsigned int strobe_pos_next1;
        uint64_t strobe_hashval_next1;
        unsigned int strobe_pos_next2;
        uint64_t strobe_hashval_next2;

        if (i + 2*w_max <= seq_length){
            unsigned int w1_start = i+w_min;
            unsigned int w1_end = i+w_max;
            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next1, strobe_hashval_next1, w1_start, w1_end, q);

            unsigned int w2_start = i+w_max + w_min;
            unsigned int w2_end = i+2*w_max;
//            uint64_t conditional_next = strobe_hash ^ strobe_hashval_next1;
            get_next_strobe(string_hashes, strobe_hashval_next1, strobe_pos_next2, strobe_hashval_next2, w2_start, w2_end, q);
        }

        else if ((i + 2*w_min + 1 < seq_length) && (seq_length < i + 2*w_max) ){

            int overshot;
            overshot = i + 2*w_max - seq_length;
            unsigned int w1_start = i+w_min;
            unsigned int w1_end = i+w_max - overshot/2;
            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next1, strobe_hashval_next1, w1_start, w1_end, q);

            unsigned int w2_start = i+w_max - overshot/2 + w_min;
            unsigned int w2_end = i+2*w_max - overshot;
//            uint64_t conditional_next = strobe_hash ^ strobe_hashval_next1;
            get_next_strobe(string_hashes, strobe_hashval_next1, strobe_pos_next2, strobe_hashval_next2, w2_start, w2_end, q);
        }
        else{
            return randstrobes3;
        }

        uint64_t hash_randstrobe3 = (strobe_hash/3) + (strobe_hashval_next1/4) + (strobe_hashval_next2/5);

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
//        unsigned int seq_pos_strobe2 =  pos_to_seq_choord[strobe_pos_next1]; //seq_pos_strobe1 + (strobe_pos_next1 - i); //
//        unsigned int seq_pos_strobe3 =  pos_to_seq_choord[strobe_pos_next2]; //seq_pos_strobe1 + (strobe_pos_next2 - i); //
//        std::cout << i << " " << strobe_pos_next1 << " " << strobe_pos_next2 << " " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << seq_pos_strobe3 << " " << pos_to_seq_choord.size() << std::endl;

//         TODO: Take care of corner case (tmep if statement below. Some values in end of string produce a cororidnate of 0 for the last strobe. Probably an off-by-one error in the calculation of the strobe coord in the last strobe window
//        if (strobe_pos_next2 ==  seq_length){
////            std::cout << "OMGGGGGGG " << i << " " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << seq_pos_strobe3 << std::endl;
//            seq_pos_strobe3 = seq_length-1;
//        }
        std::tuple<uint64_t, unsigned int, unsigned int> s (hash_randstrobe3, ref_index, seq_pos_strobe1);
        randstrobes3.push_back(s);


//        auto strobe1 = seq.substr(i, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next1 - (i+k), ' ') << std::string(k, 'X') << std::string(strobe_pos_next2 - strobe_pos_next1 - k, ' ') << std::string(k, 'X') << std::endl;
//        std::cout << i << " " << strobe_pos_next1 << " " << strobe_pos_next2 << " " << seq_length << std::endl;
    }
    return randstrobes3;
}













