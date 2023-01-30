#include "indexparameters.hpp"
#include <vector>
#include <cmath>
#include <sstream>
#include <iostream>
#include "io.hpp"


/* Pre-defined index parameters that work well for a certain
 * "canonical" read length (and similar read lengths)  */
struct Profile {
    int canonical_read_length;
    int r_threshold;
    int k;
    int s_offset;
    int l;
    int u;
};

static auto max{std::numeric_limits<int>::max()};

static std::vector<Profile> profiles = {
        Profile{ 50,  90, 20, -4, -4,  2},
        Profile{100, 110, 20, -4, -2,  2},
        Profile{125, 135, 20, -4, -1,  4},
        Profile{150, 175, 20, -4,  1,  7},
        Profile{250, 275, 20, -4,  4, 13},
        Profile{300, 375, 22, -4,  2, 12},
        Profile{400, max, 23, -6,  2, 12},
    };

/* Create an IndexParameters instance based on a given read length.
 * c, k, s, l, u and max_seed_len can be used to override determined parameters
 * by setting them to a value other than IndexParameters::DEFAULT.
 */
IndexParameters IndexParameters::from_read_length(int read_length, int c, int k, int s, int l, int u, int max_seed_len) {
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
            if (l == DEFAULT) {
                l = p.l;
            }
            if (u == DEFAULT) {
                u = p.u;
            }
            canonical_read_length = p.canonical_read_length;
            break;
        }
    }

    int max_dist;
    if (max_seed_len == DEFAULT) {
        max_dist = std::max(static_cast<int>(canonical_read_length) - 70, k);
        max_dist = std::min(255, max_dist);
    } else {
        max_dist = max_seed_len - k; // convert to distance in start positions
    }
    int q = std::pow(2, c == DEFAULT ? default_c : c) - 1;
    return IndexParameters(canonical_read_length, k, s, l, u, q, max_dist);
}

void IndexParameters::write(std::ostream& os) const {
    write_int_to_ostream(os, canonical_read_length);
    write_int_to_ostream(os, k);
    write_int_to_ostream(os, s);
    write_int_to_ostream(os, l);
    write_int_to_ostream(os, u);
    write_int_to_ostream(os, q);
    write_int_to_ostream(os, max_dist);
}

IndexParameters IndexParameters::read(std::istream& is) {
    size_t canonical_read_length = read_int_from_istream(is);
    int k = read_int_from_istream(is);
    int s = read_int_from_istream(is);
    int l = read_int_from_istream(is);
    int u = read_int_from_istream(is);
    int q = read_int_from_istream(is);
    int max_dist = read_int_from_istream(is);
    return IndexParameters(canonical_read_length, k, s, l, u, q, max_dist);
}

bool IndexParameters::operator==(const IndexParameters& other) const {
    return
        this->canonical_read_length == other.canonical_read_length
        && this->k == other.k
        && this->s == other.s
        && this->l == other.l
        && this->u == other.u
        && this->q == other.q
        && this->max_dist == other.max_dist
        && this->t_syncmer == other.t_syncmer
        && this->w_min == other.w_min
        && this->w_max == other.w_max;
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

std::ostream& operator<<(std::ostream& os, const IndexParameters& parameters) {
    os << "IndexParameters("
        << "r=" << parameters.canonical_read_length
        << ", k=" << parameters.k
        << ", s=" << parameters.s
        << ", l=" << parameters.l
        << ", u=" << parameters.u
        << ", q=" << parameters.q
        << ", max_dist=" << parameters.max_dist
        << ", t_syncmer=" << parameters.t_syncmer
        << ", w_min=" << parameters.w_min
        << ", w_max=" << parameters.w_max
        << ")";
    return os;
}
