#include <string>

#include <deque>
#include <bitset>
#include <algorithm>
#include <cassert>
#include <xxhash.h>

#include "randstrobes.hpp"

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

static inline syncmer_hash_t syncmer_kmer_hash(uint64_t packed) {
    // return robin_hash(yk);
    // return yk;
    // return hash64(yk, mask);
    // return sahlin_dna_hash(yk, mask);
    return XXH64(&packed, sizeof(uint64_t), 0);
}

std::ostream& operator<<(std::ostream& os, const Syncmer& syncmer) {
    os << "Syncmer(hash=" << syncmer.hash << ", position=" << syncmer.position << ")";
    return os;
}

Syncmer SyncmerIterator::next() {
    for ( ; i < seq.length(); ++i) {
//    for (size_t i = 0; i < seq.length(); i++) {
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
                auto syncmer = Syncmer{syncmer_kmer_hash(yk), i - k + 1};
                i++;
                return syncmer;
            }
        } else {
            // if there is an "N", restart
            qs_min_val = UINT64_MAX;
            qs_min_pos = -1;
            l = xs[0] = xs[1] = xk[0] = xk[1] = 0;
            qs.clear();
        }
    }
    return Syncmer{0, 0}; // end marker
}

std::pair<std::vector<syncmer_hash_t>, std::vector<unsigned int>> make_string_to_hashvalues_open_syncmers_canonical(
    const std::string &seq,
    const size_t k,
    const size_t s,
    const size_t t
) {
    std::vector<syncmer_hash_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_coordinate;
    SyncmerIterator syncmer_iterator{seq, k, s, t};
    Syncmer syncmer;
    while (!(syncmer = syncmer_iterator.next()).is_end()) {
        string_hashes.push_back(syncmer.hash);
        pos_to_seq_coordinate.push_back(syncmer.position);
    }
    return make_pair(string_hashes, pos_to_seq_coordinate);
}

std::ostream& operator<<(std::ostream& os, const Randstrobe& randstrobe) {
    os << "Randstrobe(hash=" << randstrobe.hash << ", strobe1_pos=" << randstrobe.strobe1_pos << ", strobe2_pos=" << randstrobe.strobe2_pos << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const QueryRandstrobe& randstrobe) {
    os << "QueryRandstrobe(hash=" << randstrobe.hash
        << ", start=" << randstrobe.start
        << ", end=" << randstrobe.end
        << ", is_reverse=" << randstrobe.is_reverse
        << ")";

    return os;
}

Randstrobe RandstrobeIterator::get(unsigned int strobe1_start) const {
    unsigned int w_end = std::min(static_cast<size_t>(strobe1_start + w_max), string_hashes.size() - 1);

    unsigned int seq_pos_strobe1 = pos_to_seq_coordinate[strobe1_start];
    unsigned int seq_end_constraint = seq_pos_strobe1 + max_dist;

    unsigned int w_start = strobe1_start + w_min;
    uint64_t strobe_hashval = string_hashes[strobe1_start];
    uint64_t min_val = UINT64_MAX;
    unsigned int strobe_pos_next = strobe1_start; // Defaults if no nearby syncmer
    uint64_t strobe_hashval_next = string_hashes[strobe1_start];
    std::bitset<64> b;

    for (auto i = w_start; i <= w_end; i++) {
        assert(i < string_hashes.size());
        // Method 3' skew sample more for prob exact matching
        b = (strobe_hashval ^ string_hashes[i])  & q;
        uint64_t res = b.count();

        if (pos_to_seq_coordinate[i] > seq_end_constraint) {
            break;
        }

        if (res < min_val){
            min_val = res;
            strobe_pos_next = i;
            strobe_hashval_next = string_hashes[i];
        }
    }
    uint64_t hash_randstrobe2 = string_hashes[strobe1_start] + strobe_hashval_next;

    return Randstrobe { hash_randstrobe2, seq_pos_strobe1, pos_to_seq_coordinate[strobe_pos_next] };
}

Randstrobe RandstrobeIterator2::next() {
    while (syncmers.size() <= w_max) {
        Syncmer syncmer = syncmer_iterator.next();
        if (syncmer.is_end()) {
            break;
        }
        syncmers.push_back(syncmer);
    }
    if (syncmers.size() <= w_min) {
        return RandstrobeIterator2::end();
    }
    auto strobe1 = syncmers[0];
    auto max_position = strobe1.position + max_dist;
    uint64_t min_val = UINT64_MAX;
    Syncmer strobe2 = syncmers[0]; // Defaults if no nearby syncmer

    for (auto i = w_min; i < syncmers.size() && syncmers[i].position <= max_position; i++) {
        assert(i <= w_max);
        // Method 3' skew sample more for prob exact matching
        std::bitset<64> b;
        b = (strobe1.hash ^ syncmers[i].hash) & q;
        uint64_t res = b.count();
        if (res < min_val) {
            min_val = res;
            strobe2 = syncmers[i];
        }
    }
    syncmers.pop_front();
    return Randstrobe{strobe1.hash + strobe2.hash, static_cast<unsigned int>(strobe1.position), static_cast<unsigned int>(strobe2.position)};
}

/*
 * Generate randstrobes for a query sequence and its reverse complement.
 */
QueryRandstrobeVector randstrobes_query(const std::string& seq, const IndexParameters& parameters) {
    QueryRandstrobeVector randstrobes;
    if (seq.length() < parameters.w_max) {
        return randstrobes;
    }

    // Generate syncmers for the forward sequence
    auto [string_hashes, pos_to_seq_coordinate] = make_string_to_hashvalues_open_syncmers_canonical(
        seq, parameters.k, parameters.s, parameters.t_syncmer
    );
    if (string_hashes.empty()) {
        return randstrobes;
    }

    // Generate randstrobes for the forward sequence
    RandstrobeIterator randstrobe_fwd_iter{
        string_hashes, pos_to_seq_coordinate, parameters.w_min, parameters.w_max, parameters.q, parameters.max_dist
    };
    while (randstrobe_fwd_iter.has_next()) {
        auto randstrobe = randstrobe_fwd_iter.next();
        randstrobes.push_back(
            QueryRandstrobe{
                randstrobe.hash, randstrobe.strobe1_pos, randstrobe.strobe2_pos + parameters.k, false
            }
        );
    }

    // For the reverse complement, we can re-use the syncmers of the forward
    // sequence because canonical syncmers are invariant under reverse
    // complementing. Only the coordinates need to be adjusted.
    std::reverse(string_hashes.begin(), string_hashes.end());
    std::reverse(pos_to_seq_coordinate.begin(), pos_to_seq_coordinate.end());
    for (size_t i = 0; i < string_hashes.size(); i++) {
        pos_to_seq_coordinate[i] = seq.length() - pos_to_seq_coordinate[i] - parameters.k;
    }

    // Randstrobes cannot be re-used for the reverse complement:
    // If in the forward direction, syncmer[i] and syncmer[j] were paired up, it
    // is not necessarily the case that syncmer[j] is going to be paired with
    // syncmer[i] in the reverse direction because i is fixed in the forward
    // direction and j is fixed in the reverse direction.
    RandstrobeIterator randstrobe_rc_iter{
        string_hashes, pos_to_seq_coordinate, parameters.w_min, parameters.w_max, parameters.q, parameters.max_dist
    };
    while (randstrobe_rc_iter.has_next()) {
        auto randstrobe = randstrobe_rc_iter.next();
        randstrobes.push_back(
            QueryRandstrobe{
                randstrobe.hash, randstrobe.strobe1_pos, randstrobe.strobe2_pos + parameters.k, true
            }
        );
    }
    return randstrobes;
}
