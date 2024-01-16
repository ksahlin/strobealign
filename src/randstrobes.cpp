#include <string>
#include <deque>
#include <bitset>
#include <algorithm>
#include <cassert>

#include "hash.hpp"
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
    return xxh64(packed);
}

static inline syncmer_hash_t syncmer_smer_hash(uint64_t packed) {
    // return ys;
    // return robin_hash(ys);
    // return hash64(ys, mask);
    return xxh64(packed);
}

static inline randstrobe_hash_t randstrobe_hash(syncmer_hash_t hash1, syncmer_hash_t hash2, size_t aux_len) {
    // Make the function symmetric
    if (hash1 > hash2) {
        std::swap(hash1, hash2);
    }
    return ((hash1 >> aux_len) << aux_len) ^ (hash2 >> (64 - aux_len));
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
            uint64_t hash_s = syncmer_smer_hash(ys);
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

std::vector<Syncmer> canonical_syncmers(
    const std::string_view seq,
    SyncmerParameters parameters
) {
    std::vector<Syncmer> syncmers;
    SyncmerIterator syncmer_iterator{seq, parameters};
    Syncmer syncmer;
    while (!(syncmer = syncmer_iterator.next()).is_end()) {
        syncmers.push_back(syncmer);
    }
    return syncmers;
}

std::ostream& operator<<(std::ostream& os, const Randstrobe& randstrobe) {
    os << "Randstrobe(hash=" << randstrobe.hash << ", strobe1_pos=" << randstrobe.strobe1_pos << ", strobe2_pos="
       << randstrobe.strobe2_pos << ", main_is_first=" << randstrobe.main_is_first << ")";
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

Randstrobe make_randstrobe(Syncmer strobe1, Syncmer strobe2, int aux_len) {
    bool main_is_first = strobe1.hash < strobe2.hash;
    return Randstrobe{
        randstrobe_hash(strobe1.hash, strobe2.hash, aux_len),
        static_cast<uint32_t>(strobe1.position),
        static_cast<uint32_t>(strobe2.position),
        main_is_first
    };
}

Randstrobe RandstrobeIterator::get(unsigned int strobe1_index) const {
    unsigned int w_end = std::min(static_cast<size_t>(strobe1_index + w_max), syncmers.size() - 1);

    auto strobe1 = syncmers[strobe1_index];
    auto max_position = strobe1.position + max_dist;
    unsigned int w_start = strobe1_index + w_min;
    uint64_t min_val = std::numeric_limits<uint64_t>::max();
    Syncmer strobe2 = strobe1;

    for (auto i = w_start; i <= w_end && syncmers[i].position <= max_position; i++) {
        assert(i < syncmers.size());

        // Method 3' skew sample more for prob exact matching
        std::bitset<64> b = (strobe1.hash ^ syncmers[i].hash) & q;
        uint64_t res = b.count();

        if (res < min_val) {
            min_val = res;
            strobe2 = syncmers[i];
        }
    }

    return make_randstrobe(strobe1, strobe2, aux_len);
}

Randstrobe RandstrobeGenerator::next() {
    while (syncmers.size() <= w_max) {
        Syncmer syncmer = syncmer_iterator.next();
        if (syncmer.is_end()) {
            break;
        }
        syncmers.push_back(syncmer);
    }
    if (syncmers.size() <= w_min) {
        return RandstrobeGenerator::end();
    }
    auto strobe1 = syncmers[0];
    auto max_position = strobe1.position + max_dist;
    uint64_t min_val = std::numeric_limits<uint64_t>::max();
    Syncmer strobe2 = strobe1; // Default if no nearby syncmer

    for (auto i = w_min; i < syncmers.size() && syncmers[i].position <= max_position; i++) {
        assert(i <= w_max);
        // Method 3' skew sample more for prob exact matching
        std::bitset<64> b = (strobe1.hash ^ syncmers[i].hash) & q;
        uint64_t res = b.count();

        if (res < min_val) {
            min_val = res;
            strobe2 = syncmers[i];
        }
    }
    syncmers.pop_front();

    return make_randstrobe(strobe1, strobe2, aux_len);
}

/*
 * Generate randstrobes for a query sequence and its reverse complement.
 */
QueryRandstrobeVector randstrobes_query(const std::string_view seq, const IndexParameters& parameters) {
    QueryRandstrobeVector randstrobes;
    if (seq.length() < parameters.randstrobe.w_max) {
        return randstrobes;
    }

    // Generate syncmers for the forward sequence
    auto syncmers = canonical_syncmers(seq, parameters.syncmer);
    if (syncmers.empty()) {
        return randstrobes;
    }

    // Generate randstrobes for the forward sequence
    RandstrobeIterator randstrobe_fwd_iter{syncmers, parameters.randstrobe};
    while (randstrobe_fwd_iter.has_next()) {
        auto randstrobe = randstrobe_fwd_iter.next();
        const unsigned int partial_start = randstrobe.main_is_first ? randstrobe.strobe1_pos : randstrobe.strobe2_pos;
        randstrobes.push_back(
            QueryRandstrobe {
                randstrobe.hash, randstrobe.strobe1_pos, randstrobe.strobe2_pos + parameters.syncmer.k,
                partial_start, partial_start + parameters.syncmer.k,  false
            }
        );
    }

    // For the reverse complement, we can re-use the syncmers of the forward
    // sequence because canonical syncmers are invariant under reverse
    // complementing. Only the coordinates need to be adjusted.
    std::reverse(syncmers.begin(), syncmers.end());
    for (size_t i = 0; i < syncmers.size(); i++) {
        syncmers[i].position = seq.length() - syncmers[i].position - parameters.syncmer.k;
    }

    // Randstrobes cannot be re-used for the reverse complement:
    // If in the forward direction, syncmer[i] and syncmer[j] were paired up, it
    // is not necessarily the case that syncmer[j] is going to be paired with
    // syncmer[i] in the reverse direction because i is fixed in the forward
    // direction and j is fixed in the reverse direction.
    RandstrobeIterator randstrobe_rc_iter{syncmers, parameters.randstrobe};
    while (randstrobe_rc_iter.has_next()) {
        auto randstrobe = randstrobe_rc_iter.next();
        bool main_is_first = randstrobe.main_is_first;
        const unsigned int partial_start = main_is_first ? randstrobe.strobe1_pos : randstrobe.strobe2_pos;
        randstrobes.push_back(
            QueryRandstrobe {
                randstrobe.hash, randstrobe.strobe1_pos, randstrobe.strobe2_pos + parameters.syncmer.k,
                partial_start, partial_start + parameters.syncmer.k, true
            }
        );
    }
    return randstrobes;
}
