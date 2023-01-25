#include "indexparameters.hpp"
#include <vector>
#include <limits>
#include <cmath>
#include <iostream>
#include "io.hpp"

/* Create an IndexParameters instance based on a given read length.
 * c, k, s and max_seed_len can be used to override determined parameters
 * by setting them to a value other than -1.
 */
IndexParameters IndexParameters::from_read_length(int read_length, int c, int k, int s, int max_seed_len) {
    int l, u;
    const int default_c = 8;
    struct settings {
        int r_threshold;
        int k;
        int s_offset;
        int l;
        int u;
    };
    std::vector<settings> d = {
        settings {75, 20, -4, -4, 2},
        settings {125, 20, -4, -2, 2},
        settings {175, 20, -4, 1, 7},
        settings {275, 20, -4, 4, 13},
        settings {375, 22, -4, 2, 12},
        settings {std::numeric_limits<int>::max(), 23, -6, 2, 12},
    };
    for (const auto& v : d) {
        if (read_length <= v.r_threshold) {
            if (k == -1) {
                k = v.k;
            }
            if (s == -1) {
                s = k + v.s_offset;
            }
            l = v.l;
            u = v.u;
            break;
        }
    }

    int max_dist;
    if (max_seed_len == -1) {
        max_dist = std::max(read_length - 70, k);
        max_dist = std::min(255, max_dist);
    } else {
        max_dist = max_seed_len - k; // convert to distance in start positions
    }
    int q = std::pow(2, c == -1 ? default_c : c) - 1;
    return IndexParameters(k, s, l, u, q, max_dist);
}

void IndexParameters::write(std::ostream& os) const {
    write_int_to_ostream(os, k);
    write_int_to_ostream(os, s);
    write_int_to_ostream(os, l);
    write_int_to_ostream(os, u);
    write_int_to_ostream(os, q);
    write_int_to_ostream(os, max_dist);
}

IndexParameters IndexParameters::read(std::istream& is) {
    int k = read_int_from_istream(is);
    int s = read_int_from_istream(is);
    int l = read_int_from_istream(is);
    int u = read_int_from_istream(is);
    int q = read_int_from_istream(is);
    int max_dist = read_int_from_istream(is);
    return IndexParameters(k, s, l, u, q, max_dist);
}

bool IndexParameters::operator==(const IndexParameters& other) const {
    return
        this->k == other.k
        && this->s == other.s
        && this->l == other.l
        && this->u == other.u
        && this->q == other.q
        && this->max_dist == other.max_dist
        && this->t_syncmer == other.t_syncmer
        && this->w_min == other.w_min
        && this->w_max == other.w_max;
}

// Return a parameter-specific filename extension. Example: ".r100.sti"
std::string IndexParameters::filename_extension() const {
    return ".sti";
}

std::ostream& operator<<(std::ostream& os, const IndexParameters& parameters) {
    os << "IndexParameters("
        << "k=" << parameters.k
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
