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

static inline uint64_t syncmer_kmer_hash(uint64_t packed) {
    // return robin_hash(yk);
    // return yk;
    // return hash64(yk, mask);
    // return sahlin_dna_hash(yk, mask);
    return XXH64(&packed, sizeof(uint64_t), 0);
}

static inline void make_string_to_hashvalues_open_syncmers_canonical(
    const std::string &seq,
    std::vector<uint64_t> &string_hashes,
    std::vector<unsigned int> &pos_to_seq_coordinate,
    const size_t k,
    const size_t s,
    const size_t t
) {
    const uint64_t kmask = (1ULL << 2*k) - 1;
    const uint64_t smask = (1ULL << 2*s) - 1;
    const uint64_t kshift = (k - 1) * 2;
    const uint64_t sshift = (s - 1) * 2;
    std::deque<uint64_t> qs;  // s-mer hashes
    uint64_t qs_min_val = UINT64_MAX;
    size_t qs_min_pos = -1;

    size_t l = 0;
    uint64_t xk[] = {0, 0};
    uint64_t xs[] = {0, 0};
    for (size_t i = 0; i < seq.length(); i++) {
        int c = seq_nt4_table[(uint8_t) seq[i]];
        if (c < 4) { // not an "N" base
            xk[0] = (xk[0] << 2 | c) & kmask;                  // forward strand
            xk[1] = xk[1] >> 2 | (uint64_t)(3 - c) << kshift;  // reverse strand
            xs[0] = (xs[0] << 2 | c) & smask;                  // forward strand
            xs[1] = xs[1] >> 2 | (uint64_t)(3 - c) << sshift;  // reverse strand
            if (++l < s) {
                continue;
            }
            // we find an s-mer
            uint64_t ys = std::min(xs[0], xs[1]);
//          uint64_t hash_s = robin_hash(ys);
            uint64_t hash_s = ys;
//          uint64_t hash_s = hash64(ys, mask);
//          uint64_t hash_s = XXH64(&ys, 8,0);
            qs.push_back(hash_s);
            // not enough hashes in the queue, yet
            if (qs.size() < k - s + 1) {
                continue;
            }
            if (qs.size() == k - s + 1) { // We are at the last s-mer within the first k-mer, need to decide if we add it
                for (size_t j = 0; j < qs.size(); j++) {
                    if (qs[j] < qs_min_val) {
                        qs_min_val = qs[j];
                        qs_min_pos = i - k + j + 1;
                    }
                }
            }
            else {
                // update queue and current minimum and position
                qs.pop_front();

                if (qs_min_pos == i - k) { // we popped the previous minimizer, find new brute force
                    qs_min_val = UINT64_MAX;
                    qs_min_pos = i - s + 1;
                    for (int j = qs.size() - 1; j >= 0; j--) { //Iterate in reverse to choose the rightmost minimizer in a window
                        if (qs[j] < qs_min_val) {
                            qs_min_val = qs[j];
                            qs_min_pos = i - k + j + 1;
                        }
                    }
                } else if (hash_s < qs_min_val) { // the new value added to queue is the new minimum
                    qs_min_val = hash_s;
                    qs_min_pos = i - s + 1;
                }
            }
            if (qs_min_pos == i - k + t) { // occurs at t:th position in k-mer
                uint64_t yk = std::min(xk[0], xk[1]);
                string_hashes.push_back(syncmer_kmer_hash(yk));
                pos_to_seq_coordinate.push_back(i - k + 1);
            }
        } else {
            // if there is an "N", restart
            qs_min_val = UINT64_MAX;
            qs_min_pos = -1;
            l = xs[0] = xs[1] = xk[0] = xk[1] = 0;
            qs.clear();
        }
    }
}

struct Randstrobe {
    uint64_t hash;
    unsigned int strobe1_pos;
    unsigned int strobe2_pos;
};

class RandstrobeIterator {
public:
    RandstrobeIterator(
        const std::vector<uint64_t> &string_hashes,
        const std::vector<unsigned int> &pos_to_seq_choord,
        int w_min,
        int w_max,
        uint64_t q,
        int max_dist
    ) : string_hashes(string_hashes)
      , pos_to_seq_choord(pos_to_seq_choord)
      , w_min(w_min)
      , w_max(w_max)
      , q(q)
      , max_dist(max_dist)
    {
    }

    Randstrobe next() {
        return get(strobe1_start++);
    }

    bool has_next() {
        return (strobe1_start + w_max < string_hashes.size())
            || (strobe1_start + w_min + 1 < string_hashes.size());
    }

private:
    Randstrobe get(unsigned int strobe1_start) const {
        unsigned int strobe_pos_next;
        uint64_t strobe_hashval_next;
        unsigned int w_end;
        if (strobe1_start + w_max < string_hashes.size()) {
            w_end = strobe1_start + w_max;
        } else if (strobe1_start + w_min + 1 < string_hashes.size()) {
            w_end = string_hashes.size() - 1;
        }

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[strobe1_start];
        unsigned int seq_end_constraint = seq_pos_strobe1 + max_dist;

        unsigned int w_start = strobe1_start + w_min;
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
                break;
            }

            if (res < min_val){
                min_val = res;
                strobe_pos_next = i;
    //            std::cerr << strobe_pos_next << " " << min_val << std::endl;
                strobe_hashval_next = string_hashes[i];
            }
        }
    //    std::cerr << "Offset: " <<  strobe_pos_next - w_start << " val: " << min_val <<  ", P exact:" <<  1.0 - pow ( (float) (8-min_val)/9, strobe_pos_next - w_start) << std::endl;

        uint64_t hash_randstrobe2 = string_hashes[strobe1_start] + strobe_hashval_next;

        return Randstrobe { hash_randstrobe2, seq_pos_strobe1, pos_to_seq_choord[strobe_pos_next] };
    }

    const std::vector<uint64_t> &string_hashes;
    const std::vector<unsigned int> &pos_to_seq_choord;
    const int w_min;
    const int w_max;
    const uint64_t q;
    const unsigned int max_dist;
    unsigned int strobe1_start = 0;
};


/* Generate randstrobes for a reference sequence. The randstrobes are appended
 * to the given flat_vector, which allows to call this function repeatedly for
 * multiple reference sequences (use a different ref_index each time).
 */
void randstrobes_reference(
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

    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_coordinate;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
//    make_string_to_hashvalues_random_minimizers(seq, string_hashes, pos_to_seq_choord, k, kmask, w);

    make_string_to_hashvalues_open_syncmers_canonical(seq, string_hashes, pos_to_seq_coordinate, k, s, t);

    unsigned int nr_hashes = string_hashes.size();
    if (nr_hashes == 0) {
        return;
    }

    RandstrobeIterator randstrobe_iter { string_hashes, pos_to_seq_coordinate, w_min, w_max, q, max_dist };
    while (randstrobe_iter.has_next()) {
        auto randstrobe = randstrobe_iter.next();
        int packed = (ref_index << 8);
        packed = packed + (randstrobe.strobe2_pos - randstrobe.strobe1_pos);
        MersIndexEntry s {randstrobe.hash, randstrobe.strobe1_pos, packed};
        flat_vector.push_back(s);
    }
}

/*
 * Generate randstrobes for a query sequence (read).
 *
 * This function stores randstrobes for both directions created from canonical
 * syncmers. Since creating canonical syncmers is the most time consuming step,
 * we avoid performing it twice for the read and its reverse complement here.
 */
mers_vector_read randstrobes_query(
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
    auto read_length = seq.length();
    if (read_length < w_max) {
        return randstrobes2;
    }

    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_coordinate;

    make_string_to_hashvalues_open_syncmers_canonical(seq, string_hashes, pos_to_seq_coordinate, k, s, t);

    unsigned int nr_hashes = string_hashes.size();
    if (nr_hashes == 0) {
        return randstrobes2;
    }

    RandstrobeIterator randstrobe_fwd_iter { string_hashes, pos_to_seq_coordinate, w_min, w_max, q, max_dist };
    while (randstrobe_fwd_iter.has_next()) {
        auto randstrobe = randstrobe_fwd_iter.next();
        unsigned int offset_strobe = randstrobe.strobe2_pos - randstrobe.strobe1_pos;
        QueryMer s {randstrobe.hash, randstrobe.strobe1_pos, offset_strobe, false};
        randstrobes2.push_back(s);
    }

    std::reverse(string_hashes.begin(), string_hashes.end());
    std::reverse(pos_to_seq_coordinate.begin(), pos_to_seq_coordinate.end());
    for (unsigned int i = 0; i < nr_hashes; i++) {
        pos_to_seq_coordinate[i] = read_length - pos_to_seq_coordinate[i] - k;
    }

    RandstrobeIterator randstrobe_rc_iter { string_hashes, pos_to_seq_coordinate, w_min, w_max, q, max_dist };
    while (randstrobe_rc_iter.has_next()) {
        auto randstrobe = randstrobe_rc_iter.next();
        unsigned int offset_strobe = randstrobe.strobe2_pos - randstrobe.strobe1_pos;
        QueryMer s {randstrobe.hash, randstrobe.strobe1_pos, offset_strobe, true};
        randstrobes2.push_back(s);
    }
    return randstrobes2;
}
