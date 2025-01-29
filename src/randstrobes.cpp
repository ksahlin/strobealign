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
    syncmer_hash_t hash1, syncmer_hash_t hash2, randstrobe_hash_t main_hash_mask
) {
    return ((hash1 & main_hash_mask) | (hash2 & ~main_hash_mask)) & RANDSTROBE_HASH_MASK;
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
       << randstrobe.strobe2_pos << ", first_strobe_is_main=" << randstrobe.first_strobe_is_main << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const QueryRandstrobe& randstrobe) {
    os << "QueryRandstrobe(hash=" << randstrobe.hash
        << ", start=" << randstrobe.start
        << ", end=" << randstrobe.end
        << ")";

    return os;
}

Randstrobe make_randstrobe(Syncmer strobe1, Syncmer strobe2, randstrobe_hash_t main_hash_mask) {
    return Randstrobe{
        randstrobe_hash(strobe1.hash, strobe2.hash, main_hash_mask),
        static_cast<uint32_t>(strobe1.position),
        static_cast<uint32_t>(strobe2.position),
        true
    };
}

Randstrobe RandstrobeIterator::get(unsigned int strobe1_index) const {
    unsigned int w_end = std::min(static_cast<size_t>(strobe1_index + parameters.w_max), syncmers.size() - 1);

    auto strobe1 = syncmers[strobe1_index];
    auto max_position = strobe1.position + parameters.max_dist;
    unsigned int w_start = strobe1_index + parameters.w_min;
    uint64_t min_val = std::numeric_limits<uint64_t>::max();
    Syncmer strobe2 = strobe1;

    for (auto i = w_start; i <= w_end && syncmers[i].position <= max_position; i++) {
        assert(i < syncmers.size());

        // Method 3' skew sample more for prob exact matching
        std::bitset<64> b = (strobe1.hash ^ syncmers[i].hash) & parameters.q;
        uint64_t res = b.count();

        if (res < min_val) {
            min_val = res;
            strobe2 = syncmers[i];
        }
    }

    return make_randstrobe(strobe1, strobe2, parameters.main_hash_mask);
}

Randstrobe RandstrobeGenerator::next() {
    while (syncmers.size() <= parameters.w_max) {
        Syncmer syncmer = syncmer_iterator.next();
        if (syncmer.is_end()) {
            break;
        }
        syncmers.push_back(syncmer);
    }
    if (syncmers.size() <= parameters.w_min) {
        return RandstrobeGenerator::end();
    }
    auto strobe1 = syncmers[0];
    auto max_position = strobe1.position + parameters.max_dist;
    uint64_t min_val = std::numeric_limits<uint64_t>::max();
    Syncmer strobe2 = strobe1; // Default if no nearby syncmer

    for (auto i = parameters.w_min; i < syncmers.size() && syncmers[i].position <= max_position; i++) {
        assert(i <= parameters.w_max);
        // Method 3' skew sample more for prob exact matching
        std::bitset<64> b = (strobe1.hash ^ syncmers[i].hash) & parameters.q;
        uint64_t res = b.count();

        if (res < min_val) {
            min_val = res;
            strobe2 = syncmers[i];
        }
    }
    syncmers.pop_front();

    return make_randstrobe(strobe1, strobe2, parameters.main_hash_mask);
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
        const unsigned int partial_start = randstrobe.first_strobe_is_main ? randstrobe.strobe1_pos : randstrobe.strobe2_pos;
        randstrobes[0].push_back(
            QueryRandstrobe {
                randstrobe.hash, randstrobe.strobe1_pos, randstrobe.strobe2_pos + parameters.syncmer.k,
                partial_start, partial_start + parameters.syncmer.k
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
        bool first_strobe_is_main = randstrobe.first_strobe_is_main;
        const unsigned int partial_start = first_strobe_is_main ? randstrobe.strobe1_pos : randstrobe.strobe2_pos;
        randstrobes[1].push_back(
            QueryRandstrobe {
                randstrobe.hash, randstrobe.strobe1_pos, randstrobe.strobe2_pos + parameters.syncmer.k,
                partial_start, partial_start + parameters.syncmer.k
            }
        );
    }
    return randstrobes;
}
