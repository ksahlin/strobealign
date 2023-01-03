#ifndef REVCOMP_HPP
#define REVCOMP_HPP

#include <string>
#include <algorithm>

static inline std::string reverse_complement(const std::string &read) {
    auto read_rev = read;
    std::reverse(read_rev.begin(), read_rev.end()); // reverse
    for (size_t j = 0; j < read_rev.length(); ++j) { // complement
        if (read_rev[j] == 'A') read_rev[j] = 'T';
        else if (read_rev[j] == 'T') read_rev[j] = 'A';
        else if (read_rev[j] == 'C') read_rev[j] = 'G';
        else if (read_rev[j] == 'G') read_rev[j] = 'C';
    }
    return read_rev;
}


/* A (nucleotide) sequence and its reverse complement.
 * The reverse complement is computed on first access only
 * (and cached).
 */
class Read {
public:
    const std::string& seq;

    Read(const std::string& s) : seq(s) {
    }

    /* Return reverse complemented sequence */
    std::string rc() const {
        if (!has_reverse_complement) {
            rc_sequence = reverse_complement(seq);
            has_reverse_complement = true;
        }
        return rc_sequence;
    }

    std::string::size_type size() const {
        return seq.size();
    }

private:
    mutable std::string rc_sequence;
    mutable bool has_reverse_complement = false;
};

#endif
