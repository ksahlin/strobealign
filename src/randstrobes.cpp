#include <string>
#include <deque>
#include <bitset>
#include <algorithm>
#include <cassert>
#include <array>

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

/*
 * This function combines two individual syncmer hashes into a single hash
 * for the randstrobe.
 *
 * One syncmer is designated as the "main", the other is the "auxiliary".
 * The combined hash is obtained by setting the top bits to the bits of
 * the main hash and the bottom bits to the bits of the auxiliary
 * hash. Since entries in the index are sorted by randstrobe hash, this allows
 * us to search for the main syncmer only by masking out the lower bits.
 */
static inline randstrobe_hash_t randstrobe_hash(
    syncmer_hash_t hash1,
    syncmer_hash_t hash2,
    syncmer_hash_t hash3,
    randstrobe_hash_t main_hash_mask,
    randstrobe_hash_t strobe2_mask
) {
    randstrobe_hash_t strobe3_mask = (hash3 & ~(main_hash_mask ^ strobe2_mask));
    return ((hash1 & main_hash_mask) | (hash2 & strobe2_mask) | (hash3 & strobe3_mask)) & RANDSTROBE_HASH_MASK;
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
            if (++l < static_cast<size_t>(parameters.s)) {
                continue;
            }
            // we find an s-mer
            uint64_t ys = std::min(xs[0], xs[1]);
            uint64_t hash_s = syncmer_smer_hash(ys);
            qs.push_back(hash_s);
            // not enough hashes in the queue, yet
            if (qs.size() < static_cast<size_t>(parameters.k - parameters.s + 1)) {
                continue;
            }
            if (qs.size() == static_cast<size_t>(parameters.k - parameters.s + 1)) { // We are at the last s-mer within the first k-mer, need to decide if we add it
                for (size_t j = 0; j < qs.size(); j++) {
                    if (qs[j] <= qs_min_val) {
                        qs_min_val = qs[j];
                    }
                }
            }
            else {
                // update queue and current minimum and position
                uint64_t front = qs.front();
                qs.pop_front();

                if (front == qs_min_val) {
                    // we popped a minimum, find new brute force
                    qs_min_val = UINT64_MAX;
                    for (size_t j = 0; j < qs.size(); j++) {
                        if (qs[j] <= qs_min_val) {
                            qs_min_val = qs[j];
                        }
                    }
                } else if (hash_s < qs_min_val) { // the new value added to queue is the new minimum
                    qs_min_val = hash_s;
                }
            }
            if (qs[parameters.t_syncmer - 1] == qs_min_val) { // occurs at t:th position in k-mer
                uint64_t yk = std::min(xk[0], xk[1]);
                auto syncmer = Syncmer{syncmer_kmer_hash(yk), i - parameters.k + 1};
                i++;
                return syncmer;
            }
        } else {
            // if there is an "N", restart
            qs_min_val = UINT64_MAX;
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
       << randstrobe.strobe2_pos << ", strobe3_pos=" << randstrobe.strobe3_pos << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const QueryRandstrobe& randstrobe) {
    os << "QueryRandstrobe(hash=" << randstrobe.hash
        << ", start=" << randstrobe.start
        << ", end=" << randstrobe.end
        << ")";

    return os;
}

Randstrobe make_randstrobe(
    Syncmer strobe1,
    Syncmer strobe2,
    Syncmer strobe3,
    randstrobe_hash_t main_hash_mask,
    randstrobe_hash_t strobe2_mask
) {
    return Randstrobe{
        randstrobe_hash(strobe1.hash, strobe2.hash, strobe3.hash, main_hash_mask, strobe2_mask),
        randstrobe_hash(strobe3.hash, strobe2.hash, strobe1.hash, main_hash_mask, strobe2_mask),
        static_cast<uint32_t>(strobe1.position),
        static_cast<uint32_t>(strobe2.position),
        static_cast<uint32_t>(strobe3.position)
    };
}

Randstrobe RandstrobeIterator::get(unsigned int strobe1_index) const {
    auto strobe1 = syncmers[strobe1_index];
    auto curr_hash = strobe1.hash;
    uint strobe2_index = get_next_strobe_index(strobe1_index, strobe1_index, curr_hash, 2);
    auto strobe2 = syncmers[strobe2_index];
    curr_hash ^= strobe2.hash;
    uint strobe3_index = get_next_strobe_index(strobe1_index, strobe2_index, curr_hash, 3);
    auto strobe3 = syncmers[strobe3_index];

    return make_randstrobe(strobe1, strobe2, strobe3, parameters.main_hash_mask, parameters.strobe2_mask);
}

uint RandstrobeIterator::get_next_strobe_index(unsigned int first_strobe_index,
                                               unsigned int curr_strobe_index,
                                               uint64_t curr_hash,
                                               uint strobe_count) const {
    auto first_strobe = syncmers[first_strobe_index];
    auto max_position = first_strobe.position + parameters.max_dist;
    unsigned int w_start = first_strobe_index + parameters.w_min + (strobe_count - 2) * parameters.w_max;
    unsigned int w_end = std::min(static_cast<size_t>(first_strobe_index + (strobe_count - 1) * parameters.w_max), syncmers.size() - 1);
    uint64_t min_val = std::numeric_limits<uint64_t>::max();
    uint next_strobe_index = curr_strobe_index;

    for (auto i = w_start; i <= w_end && syncmers[i].position <= max_position; i++) {
        assert(i < syncmers.size());

        // Method 3' skew sample more for prob exact matching
        std::bitset<64> b = (curr_hash ^ syncmers[i].hash) & parameters.q;
        uint64_t res = b.count();

        if (res < min_val) {
            min_val = res;
            next_strobe_index = i;
        }
    }
    return next_strobe_index;

//    return make_randstrobe(strobe1, strobe2, parameters.main_hash_mask);
}

Randstrobe RandstrobeGenerator::next() {
    while (syncmers.size() <= 2 * parameters.w_max) {
        Syncmer syncmer = syncmer_iterator.next();
        if (syncmer.is_end()) {
            break;
        }
        syncmers.push_back(syncmer);
    }
    if (syncmers.size() <= parameters.w_min + parameters.w_max) {
        return RandstrobeGenerator::end();
    }
    auto strobe1 = syncmers[0];
    auto curr_hash = strobe1.hash;
    uint strobe2_index = get_next_strobe_index(0, 0, curr_hash, 2);
    auto strobe2 = syncmers[strobe2_index];
    curr_hash ^= strobe2.hash;
    uint strobe3_index = get_next_strobe_index(0, strobe2_index, curr_hash, 3);
    auto strobe3 = syncmers[strobe3_index];
    syncmers.pop_front();
    return make_randstrobe(strobe1, strobe2, strobe3, parameters.main_hash_mask, parameters.strobe2_mask);
}

uint RandstrobeGenerator::get_next_strobe_index(unsigned int first_strobe_index,
                                                unsigned int curr_strobe_index,
                                                uint64_t curr_hash,
                                                uint strobe_count) const {
    auto first_strobe = syncmers[first_strobe_index];
    auto max_position = first_strobe.position + parameters.max_dist;
    unsigned int w_start = first_strobe_index + parameters.w_min + (strobe_count - 2) * parameters.w_max;
    unsigned int w_end = std::min(static_cast<size_t>(first_strobe_index + (strobe_count - 1) * parameters.w_max), syncmers.size() - 1);
    uint64_t min_val = std::numeric_limits<uint64_t>::max();
    uint next_strobe_index = curr_strobe_index;

    for (auto i = w_start; i <= w_end && syncmers[i].position <= max_position; i++) {
        assert(i < syncmers.size());

        // Method 3' skew sample more for prob exact matching
        std::bitset<64> b = (curr_hash ^ syncmers[i].hash) & parameters.q;
        uint64_t res = b.count();

        if (res < min_val) {
            min_val = res;
            next_strobe_index = i;
        }
    }
    return next_strobe_index;
}

/*
 * Generate randstrobes for a query sequence and its reverse complement.
 */
std::array<std::vector<QueryRandstrobe>, 2> randstrobes_query(const std::string_view seq, const IndexParameters& parameters) {
    std::array<std::vector<QueryRandstrobe>, 2> randstrobes;
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
        randstrobes[0].push_back(
            QueryRandstrobe {
                randstrobe.hash, randstrobe.hash_revcomp, randstrobe.strobe1_pos, randstrobe.strobe3_pos + parameters.syncmer.k
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
        randstrobes[1].push_back(
            QueryRandstrobe {
                randstrobe.hash, randstrobe.hash_revcomp, randstrobe.strobe1_pos, randstrobe.strobe3_pos + parameters.syncmer.k
            }
        );
    }
    return randstrobes;
}
