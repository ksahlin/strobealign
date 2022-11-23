#include "randstrobes.hpp"
#include <string>

#include <deque>
#include <bitset>
#include <algorithm>
#include <xxhash.h>

// a, A -> 0
// c, C -> 1
// g, G -> 2
// t, T, u, U -> 3
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
};

// update queue and current minimum and position
static inline void update_window(std::deque <uint64_t> &q, std::deque <unsigned int> &q_pos, uint64_t &q_min_val, int &q_min_pos, uint64_t new_strobe_hashval, int i, bool &new_minimizer){
//    uint64_t popped_val;
//    popped_val = q.front();
    q.pop_front();

    unsigned int popped_index;
    popped_index=q_pos.front();
    q_pos.pop_front();

    q.push_back(new_strobe_hashval);
    q_pos.push_back(i);
    if (q_min_pos == popped_index){ // we popped the previous minimizer, find new brute force
//    if (popped_val == q_min_val){ // we popped the minimum value, find new brute force
        q_min_val = UINT64_MAX;
        q_min_pos = i;
        for (int j = q.size() - 1; j >= 0; j--) { //Iterate in reverse to choose the rightmost minimizer in a window
//        for (int j = 0; j <= q.size()-1; j++) {
            if (q[j] < q_min_val) {
                q_min_val = q[j];
                q_min_pos = q_pos[j];
                new_minimizer = true;
            }
        }
    }
    else if ( new_strobe_hashval < q_min_val ) { // the new value added to queue is the new minimum
        q_min_val = new_strobe_hashval;
        q_min_pos = i;
        new_minimizer = true;
    }
}


static inline void make_string_to_hashvalues_open_syncmers_canonical(const std::string &seq, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, uint64_t kmask, int k, uint64_t smask, int s, int t)
{
    std::deque<uint64_t> qs;
    std::deque<unsigned int> qs_pos;
    int seq_length = seq.length();
    int qs_size = 0;
    uint64_t qs_min_val = UINT64_MAX;
    int qs_min_pos = -1;


//    robin_hood::hash<uint64_t> robin_hash;
    uint64_t mask = (1ULL<<2*k) - 1;
//    std::cerr << mask << std::endl;

//    std::vector<std::tuple<uint64_t, unsigned int, unsigned int> > kmers;
    int gap = 0;
    std::string subseq;
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
                // queue not initialized yet
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
                    if (qs_min_pos == qs_pos[t-1]) { // occurs at t:th position in k-mer
//                    if ( (qs_min_pos == qs_pos[t-1]) || ((gap > 10) && ((qs_min_pos == qs_pos[k - s]) || (qs_min_pos == qs_pos[0]))) ) { // occurs at first or last position in k-mer
                        uint64_t yk = std::min(xk[0], xk[1]);
//                        uint64_t hash_k = robin_hash(yk);
//                        uint64_t hash_k = yk;
//                        uint64_t hash_k =  hash64(yk, mask);
                        uint64_t hash_k = XXH64(&yk, 8,0);
//                        uint64_t hash_k =  sahlin_dna_hash(yk, mask);
                        string_hashes.push_back(hash_k);
                        pos_to_seq_choord.push_back(i - k + 1);
                        hash_count++;
//                        std::cerr << i - s + 1 << " " << i - k + 1 << " " << (xk[0] < xk[1]) << std::endl;
//                        std::cerr <<  "Sampled gap: " << gap (k-s+1) << std::endl;
                        gap = 0;
                    }
                }
                else{
                    bool new_minimizer = false;
                    update_window(qs, qs_pos, qs_min_val, qs_min_pos, hash_s, i - s + 1, new_minimizer );
                    if (qs_min_pos == qs_pos[t-1]) { // occurs at t:th position in k-mer
//                    if ( (qs_min_pos == qs_pos[t-1]) || ((gap > 10) && ((qs_min_pos == qs_pos[k - s]) || (qs_min_pos == qs_pos[0]))) ) { // occurs at first or last position in k-mer
//                        if ( (gap > k) && (gap < 200) ) { // open syncmers no window guarantee, fill in subsequence with closed syncmers
//                            subseq = seq.substr(i - k + 1 - gap + 1, gap +k);
//                            make_string_to_hashvalues_closed_syncmers_canonical(subseq, string_hashes, pos_to_seq_choord, kmask, k, smask, s, t, i - k + 1 - gap + 1);
//                        }

                        uint64_t yk = std::min(xk[0], xk[1]);
//                        uint64_t hash_k = robin_hash(yk);
//                        uint64_t hash_k = yk;
//                        uint64_t hash_k = hash64(yk, mask);
                        uint64_t hash_k = XXH64(&yk, 8, 0);
//                        uint64_t hash_k =  sahlin_dna_hash(yk, mask);
                        string_hashes.push_back(hash_k);
                        pos_to_seq_choord.push_back(i - k + 1);
//                        std::cerr << i - k + 1 << std::endl;
                        hash_count++;
//                        std::cerr << i - s + 1 << " " << i - k + 1 << " " << (xk[0] < xk[1]) << std::endl;
//                        std::cerr <<  "Gap: " << gap << " position:" << i - k + 1 << std::endl;
                        gap = 0;
                    }
                    gap ++;
                }
//                if (gap > 25){
//                    std::cerr <<  "Gap: " << gap << " position:" << i - k + 1 << std::endl;
//                    if (gap < 500 ) {
//                        std::cerr << seq.substr(i - k + 1 - gap + 1, gap +k) << std::endl;
//                    }
//                }
            }
        } else {
            // if there is an "N", restart
            qs_min_val = UINT64_MAX;
            qs_min_pos = -1;
            l = xs[0] = xs[1] = xk[0] = xk[1] = 0;
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

class RandstrobeIterator {
public:
    RandstrobeIterator(
        const std::vector<uint64_t> &string_hashes,
        const std::vector<unsigned int> &pos_to_seq_choord
    ) : string_hashes(string_hashes)
      , pos_to_seq_choord(pos_to_seq_choord)
    {
    }

    void get_next_strobe_dist_constraint(
        unsigned int &strobe_pos_next,  // Output
        uint64_t &strobe_hashval_next,  // Output
        unsigned int w_start,
        unsigned int w_end,
        uint64_t q,
        unsigned int seq_start,
        unsigned int seq_end_constraint,
        unsigned int strobe1_start
    ) {
        uint64_t strobe_hashval = string_hashes[strobe1_start];

        uint64_t min_val = UINT64_MAX;
        strobe_pos_next = strobe1_start; // Defaults if no nearby syncmer
        strobe_hashval_next = string_hashes[strobe1_start];
        std::bitset<64> b;

        for (auto i = w_start; i <= w_end; i++) {

            // Method 3' skew sample more for prob exact matching
            b = (strobe_hashval ^ string_hashes[i])  & q;
            uint64_t res = b.count();

            if (pos_to_seq_choord[i] > seq_end_constraint){
                return;
            }

            if (res < min_val){
                min_val = res;
                strobe_pos_next = i;
    //            std::cerr << strobe_pos_next << " " << min_val << std::endl;
                strobe_hashval_next = string_hashes[i];
            }
        }
    //    std::cerr << "Offset: " <<  strobe_pos_next - w_start << " val: " << min_val <<  ", P exact:" <<  1.0 - pow ( (float) (8-min_val)/9, strobe_pos_next - w_start) << std::endl;

    }

private:
    const std::vector<uint64_t> &string_hashes;
    const std::vector<unsigned int> &pos_to_seq_choord;
};


void seq_to_randstrobes2(
    ind_mers_vector& flat_vector,
    int k,
    int w_min,
    int w_max,
    const std::string &seq,
    int ref_index,
    int s,
    int t,
    uint64_t q,
    int max_dist
) {
    if (seq.length() < w_max) {
        return;
    }

    uint64_t kmask=(1ULL<<2*k) - 1;
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
//    make_string_to_hashvalues_random_minimizers(seq, string_hashes, pos_to_seq_choord, k, kmask, w);

//    int s = k-4;
//    int t = 3;
    uint64_t smask=(1ULL<<2*s) - 1;
    make_string_to_hashvalues_open_syncmers_canonical(seq, string_hashes, pos_to_seq_choord, kmask, k, smask, s, t);
//    make_string_to_hashvalues_open_syncmers(seq, string_hashes, pos_to_seq_choord, kmask, k, smask, s, t);

    unsigned int nr_hashes = string_hashes.size();
    if (nr_hashes == 0) {
        return;
    }

//        for (unsigned int i = 0; i < seq_length; i++) {
//        std::cerr << "REF POS INDEXED: " << pos_to_seq_choord[i] << " OTHER DIRECTION: " << seq.length() - pos_to_seq_choord[i] - k <<std::endl;
//    }
//    int tmp_cnt = 0;
//    for (auto a: pos_to_seq_choord){
//        std::cerr << " OK: " << tmp_cnt << " " << a << std::endl;
//        tmp_cnt ++;
//    }
//    std::cerr << seq << std::endl;

    RandstrobeIterator randstrobe_iter { string_hashes, pos_to_seq_choord };

    // create the randstrobes
    for (unsigned int i = 0; i < nr_hashes; i++) {

//        if ((i % 1000000) == 0 ){
//            std::cerr << i << " strobemers created." << std::endl;
//        }
        unsigned int strobe_pos_next;
        uint64_t strobe_hashval_next;
        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_end = seq_pos_strobe1 + max_dist;

        unsigned int w_start = i+w_min;
        if (i + w_max < nr_hashes){
            unsigned int w_end = i+w_max;
            randstrobe_iter.get_next_strobe_dist_constraint(strobe_pos_next, strobe_hashval_next, w_start, w_end, q, seq_pos_strobe1, seq_end,i);
        }
        else if (i + w_min + 1 < nr_hashes) {
            unsigned int w_end = nr_hashes - 1;
            randstrobe_iter.get_next_strobe_dist_constraint(strobe_pos_next, strobe_hashval_next, w_start, w_end, q, seq_pos_strobe1, seq_end, i);
        }
        else {
            return;
        }

//        uint64_t hash_randstrobe2 = (string_hashes[i] << k) ^ strobe_hashval_next;
        uint64_t hash_randstrobe2 = (string_hashes[i]) + (strobe_hashval_next);

        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe_pos_next];
//        std::cerr <<  "Seed length: " << seq_pos_strobe2 + k - seq_pos_strobe1 << std::endl;

        // UNTIL HERE roughly identical

        int packed = (ref_index << 8);
//        int offset_strobe =  seq_pos_strobe2 - seq_pos_strobe1;
        packed = packed + (seq_pos_strobe2 - seq_pos_strobe1);
        MersIndexEntry s {hash_randstrobe2, seq_pos_strobe1, packed};
        flat_vector.push_back(s);
//        std::cerr << seq_pos_strobe1 << " " << seq_pos_strobe2 << std::endl;
//        std::cerr << "FORWARD REF: " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << hash_randstrobe2 << std::endl;
//        std::cerr << "REFERENCE: " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << hash_randstrobe2 << std::endl;


//        auto strobe1 = seq.substr(seq_pos_strobe1, k);
//        auto strobe2 = seq.substr(seq_pos_strobe2, k);
//        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe_pos_next];
//        if (seq_pos_strobe2 > (seq_pos_strobe1+k)) {
////            std::cerr << seq_pos_strobe1 << " LOOOOL " << seq_pos_strobe2 << " " << seq_pos_strobe2 - (seq_pos_strobe1+k) << std::endl;
//            std::cerr << std::string(seq_pos_strobe1, ' ') << strobe1 << std::string(seq_pos_strobe2 - (seq_pos_strobe1+k), ' ') << strobe2 << std::endl;
//        }
//        std::cerr << i << " " << strobe_pos_next << " " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << seq_pos_strobe2 - (seq_pos_strobe1+k) << " " << hash_randstrobe2 << std::endl;
//        std::cerr << i << " " << strobe_pos_next << std::endl;

    }
//    std::cerr << randstrobes2.size() << " randstrobes generated" << '\n';
}


mers_vector_read seq_to_randstrobes2_read(
    int k,
    int w_min,
    int w_max,
    const std::string& seq,
    int s,
    int t,
    uint64_t q,
    int max_dist
) {
    // this function differs from  the function seq_to_randstrobes2 which creating randstrobes for the reference.
    // The seq_to_randstrobes2 stores randstobes only in one direction from canonical syncmers.
    // this function stores randstobes from both directions created from canonical syncmers.
    // Since creating canonical syncmers is the most time consuming step, we avoid perfomring it twice for the read and its RC here
    mers_vector_read randstrobes2;
    unsigned int read_length = seq.length();
    if (read_length < w_max) {
        return randstrobes2;
    }

    uint64_t kmask=(1ULL<<2*k) - 1;
//    std::bitset<64> x(q);
//    std::cerr << x << '\n';
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
//    make_string_to_hashvalues_random_minimizers(seq, string_hashes, pos_to_seq_choord, k, kmask, w);

//    int t = 3;
    uint64_t smask=(1ULL<<2*s) - 1;
    make_string_to_hashvalues_open_syncmers_canonical(seq, string_hashes, pos_to_seq_choord, kmask, k, smask, s, t);
//    make_string_to_hashvalues_open_syncmers(seq, string_hashes, pos_to_seq_choord, kmask, k, smask, s, t);

    unsigned int nr_hashes = string_hashes.size();
//    std::cerr << nr_hashes << " okay" << std::endl;
    if (nr_hashes == 0) {
        return randstrobes2;
    }

//    int tmp_cnt = 0;
//    for (auto a: pos_to_seq_choord){
//        std::cerr << " OK: " << tmp_cnt << " " << a << std::endl;
//        tmp_cnt ++;
//    }
//    std::cerr << seq << std::endl;

    // create the randstrobes FW direction!
    RandstrobeIterator randstrobe_fwd_iter { string_hashes, pos_to_seq_choord };
    for (unsigned int i = 0; i < nr_hashes; i++) {
        unsigned int strobe_pos_next;
        uint64_t strobe_hashval_next;
        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_end = seq_pos_strobe1 + max_dist;
        unsigned int w_start = i+w_min;
        if (i + w_max < nr_hashes){
            unsigned int w_end = i+w_max;
            // writes to strobe_pos_next, strobe_hashval_next
            randstrobe_fwd_iter.get_next_strobe_dist_constraint(strobe_pos_next, strobe_hashval_next, w_start, w_end, q, seq_pos_strobe1, seq_end, i);
        }
        else if (i + w_min + 1 < nr_hashes) {
            unsigned int w_end = nr_hashes -1;
            // writes to strobe_pos_next, strobe_hashval_next
            randstrobe_fwd_iter.get_next_strobe_dist_constraint(strobe_pos_next, strobe_hashval_next, w_start, w_end, q, seq_pos_strobe1, seq_end, i);
        }
        else {
            break;
        }

//        uint64_t hash_randstrobe2 = (string_hashes[i] << k) ^ strobe_hashval_next;
        uint64_t hash_randstrobe2 = (string_hashes[i]) + (strobe_hashval_next);

        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe_pos_next];

        // UNTIL HERE roughly identical

        // Output: from the above: seq_pos_strobe2, seq_pos_strobe1, hash_randstrobe2

        unsigned int offset_strobe =  seq_pos_strobe2 - seq_pos_strobe1;
        QueryMer s {hash_randstrobe2, seq_pos_strobe1, offset_strobe, false};
        randstrobes2.push_back(s);
//        std::cerr << "FORWARD: " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << hash_randstrobe2 << std::endl;

//        auto strobe1 = seq.substr(seq_pos_strobe1, k);
//        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe_pos_next];
//        if (seq_pos_strobe2 > (seq_pos_strobe1+k)) {
////            std::cerr << seq_pos_strobe1 << " LOOOOL " << seq_pos_strobe2 << " " << seq_pos_strobe2 - (seq_pos_strobe1+k) << std::endl;
//            std::cerr << std::string(seq_pos_strobe1, ' ') << strobe1 << std::string(seq_pos_strobe2 - (seq_pos_strobe1+k), ' ') << std::string(k, 'X') << std::endl;
//        }
//        std::cerr << i << " " << strobe_pos_next << " " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << seq_pos_strobe2 - (seq_pos_strobe1+k) << " " << hash_randstrobe2 << std::endl;
//        std::cerr << i << " " << strobe_pos_next << std::endl;

    }

    // create the randstrobes Reverse direction!
    std::reverse(string_hashes.begin(), string_hashes.end());
    std::reverse(pos_to_seq_choord.begin(), pos_to_seq_choord.end());
    for (unsigned int i = 0; i < nr_hashes; i++) {
        pos_to_seq_choord[i] = read_length - pos_to_seq_choord[i] - k;
    }

    RandstrobeIterator randstrobe_rc_iter { string_hashes, pos_to_seq_choord };
    for (unsigned int i = 0; i < nr_hashes; i++) {
        unsigned int strobe_pos_next;
        uint64_t strobe_hashval_next;
        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_end = seq_pos_strobe1 + max_dist;

        unsigned int w_start = i + w_min;
        if (i + w_max < nr_hashes){
            unsigned int w_end = i+w_max;
            randstrobe_rc_iter.get_next_strobe_dist_constraint(strobe_pos_next, strobe_hashval_next, w_start, w_end, q, seq_pos_strobe1, seq_end, i);

        }
        else if ((i + w_min + 1 < nr_hashes) && (nr_hashes <= i + w_max) ){
            unsigned int w_end = nr_hashes -1;
            randstrobe_rc_iter.get_next_strobe_dist_constraint(strobe_pos_next, strobe_hashval_next, w_start, w_end, q, seq_pos_strobe1, seq_end, i);
        }
        else{
            return randstrobes2;
        }

//        uint64_t hash_randstrobe2 = (string_hashes[i] << k) ^ strobe_hashval_next;
        uint64_t hash_randstrobe2 = (string_hashes[i]) + (strobe_hashval_next);

        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe_pos_next];
        unsigned int offset_strobe =  seq_pos_strobe2 - seq_pos_strobe1;
        QueryMer s {hash_randstrobe2, seq_pos_strobe1, offset_strobe, true};
        randstrobes2.push_back(s);

//        auto strobe1 = seq.substr(seq_pos_strobe1, k);
//        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe_pos_next];
//        if (seq_pos_strobe2 > (seq_pos_strobe1+k)) {
////            std::cerr << seq_pos_strobe1 << " LOOOOL " << seq_pos_strobe2 << " " << seq_pos_strobe2 - (seq_pos_strobe1+k) << std::endl;
//            std::cerr << std::string(seq_pos_strobe1, ' ') << strobe1 << std::string(seq_pos_strobe2 - (seq_pos_strobe1+k), ' ') << std::string(k, 'X') << std::endl;
//        }
//        std::cerr << i << " " << strobe_pos_next << " " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << seq_pos_strobe2 - (seq_pos_strobe1+k) << " " << hash_randstrobe2 << std::endl;
//        std::cerr << i << " " << strobe_pos_next << std::endl;

    }
    return randstrobes2;
}
