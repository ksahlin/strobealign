#include "indexparameters.hpp"
#include <vector>
#include <cmath>
#include <sstream>
#include <iostream>
#include "io.hpp"


bool SyncmerParameters::operator==(const SyncmerParameters& other) const {
    return this->s == other.s
        && this->k == other.k
        && this->t_syncmer == other.t_syncmer;
}

bool RandstrobeParameters::operator==(const RandstrobeParameters& other) const {
    return this->q == other.q
        && this->max_dist == other.max_dist
        && this->w_min == other.w_min
        && this->w_max == other.w_max
        && this->main_hash_mask == other.main_hash_mask
        && this->strobe2_mask == other.strobe2_mask;
}

/* Pre-defined index parameters that work well for a certain
 * "canonical" read length (and similar read lengths)  */
struct Profile {
    int canonical_read_length;
    int r_threshold;
    int k;
    int s_offset;
    int w_min;
    int w_max;
};

static auto max{std::numeric_limits<int>::max()};

static std::vector<Profile> profiles = {
        Profile{ 50,  70, 18, -4,  1,  4},
        Profile{ 75,  90, 20, -4,  1,  6},
        Profile{100, 110, 20, -4,  2,  6},
        Profile{125, 135, 20, -4,  3,  8},
        Profile{150, 175, 20, -4,  5, 11},
        Profile{250, 375, 22, -4,  6, 16},
        Profile{400, max, 23, -6,  5, 15},
    };

/* Create an IndexParameters instance based on a given read length.
 * k, s, l, u, c and max_seed_len can be used to override determined parameters
 * by setting them to a value other than IndexParameters::DEFAULT.
 */
IndexParameters IndexParameters::from_read_length(int read_length, int k, int s, int w_min, int w_max, int c, int max_seed_len, int aux_len) {
    const int default_c = 8;
    size_t canonical_read_length = 50;
    for (const auto& p : profiles) {
        if (read_length <= p.r_threshold) {
            if (k == DEFAULT) {
                k = p.k;
            }
            if (s == DEFAULT) {
                s = k + p.s_offset;
            }
            if (w_min == DEFAULT) {
                w_min = p.w_min;
            }
            if (w_max == DEFAULT) {
                w_max = p.w_max;
            }
            canonical_read_length = p.canonical_read_length;
            break;
        }
    }

    int max_dist;
    if (max_seed_len == DEFAULT) {
        max_dist = std::clamp(static_cast<int>(canonical_read_length) - 70, k, 255);
    } else {
        max_dist = max_seed_len - k; // convert to distance in start positions
    }
    int q = std::pow(2, c == DEFAULT ? default_c : c) - 1;
    if (aux_len == DEFAULT) {
        aux_len = 8;
    }

    return IndexParameters(canonical_read_length, k, s, w_min, w_max, q, max_dist, aux_len);
}

void IndexParameters::write(std::ostream& os) const {
    write_int_to_ostream(os, canonical_read_length);
    write_int_to_ostream(os, syncmer.k);
    write_int_to_ostream(os, syncmer.s);
    write_int_to_ostream(os, randstrobe.w_min);
    write_int_to_ostream(os, randstrobe.w_max);
    write_int_to_ostream(os, randstrobe.q);
    write_int_to_ostream(os, randstrobe.max_dist);
    write_int_to_ostream(os, randstrobe.main_hash_mask);
    write_int_to_ostream(os, randstrobe.strobe2_mask);
}

IndexParameters IndexParameters::read(std::istream& is) {
    size_t canonical_read_length = read_int_from_istream(is);
    int k = read_int_from_istream(is);
    int s = read_int_from_istream(is);
    const SyncmerParameters syncmer_parameters{k, s};

    uint32_t w_min = read_int_from_istream(is);
    uint32_t w_max = read_int_from_istream(is);
    uint64_t q = read_int_from_istream(is);
    int max_dist = read_int_from_istream(is);
    uint64_t main_hash_mask = read_uint64_from_istream(is);
    uint64_t strobe2_mask = read_uint64_from_istream(is);
    const RandstrobeParameters randstrobe_parameters{q, max_dist, w_min, w_max, main_hash_mask, strobe2_mask};

    return IndexParameters(canonical_read_length, syncmer_parameters, randstrobe_parameters);
}

bool IndexParameters::operator==(const IndexParameters& other) const {
    return this->canonical_read_length == other.canonical_read_length
        && this->syncmer == other.syncmer
        && this->randstrobe == other.randstrobe;
}

/*
 * Return a parameter-specific filename extension such as ".r100.sti"
 * If any of the parameters deviate from the defaults for the current
 * canonical read length, the returned extension is just ".sti".
 */
std::string IndexParameters::filename_extension() const {
    std::stringstream sstream;
    if (*this == from_read_length(canonical_read_length)) {
        // nothing was overridden
        sstream << ".r" << canonical_read_length;
    }
    sstream << ".sti";
    return sstream.str();
}

std::ostream& operator<<(std::ostream& os, const SyncmerParameters& parameters) {
    os << "SyncmerParameters("
        << "k=" << parameters.k
        << ", s=" << parameters.s
        << ", t_syncmer=" << parameters.t_syncmer
        << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const RandstrobeParameters& parameters) {
    os << "RandstrobeParameters("
        << "q=" << parameters.q
        << ", max_dist=" << parameters.max_dist
        << ", w_min=" << parameters.w_min
        << ", w_max=" << parameters.w_max
        << ", main_hash_mask=0x" << std::hex << parameters.main_hash_mask << std::dec
        << ", strobe2_mask=0x" << std::hex << parameters.strobe2_mask << std::dec
        << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const IndexParameters& parameters) {
    os << "IndexParameters("
        << "r=" << parameters.canonical_read_length
        << ", k=" << parameters.syncmer.k
        << ", s=" << parameters.syncmer.s
        << ", t_syncmer=" << parameters.syncmer.t_syncmer
        << ", q=" << parameters.randstrobe.q
        << ", max_dist=" << parameters.randstrobe.max_dist
        << ", w_min=" << parameters.randstrobe.w_min
        << ", w_max=" << parameters.randstrobe.w_max
        << ", main_hash_mask=0x" << std::hex << parameters.randstrobe.main_hash_mask << std::dec
        << ", strobe2_mask=0x" << std::hex << parameters.randstrobe.strobe2_mask << std::dec
        << ")";
    return os;
}
