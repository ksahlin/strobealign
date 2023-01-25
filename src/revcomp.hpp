#ifndef REVCOMP_HPP
#define REVCOMP_HPP

#include <string>
#include <algorithm>

// a, A -> T
// c, C -> G
// g, G -> C
// t, T, u, U -> A
static unsigned char revcomp_table[256] = {
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G',  'N', 'N', 'N', 'C',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'A', 'A', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G',  'N', 'N', 'N', 'C',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'A', 'A', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
};

static inline std::string reverse_complement(const std::string &sequence) {
    std::string result;
    result.reserve(sequence.size());
    for (size_t i = 0; i < sequence.length(); ++i) {
        auto ch = sequence[sequence.length() - i - 1];
        result.push_back(revcomp_table[static_cast<int>(ch)]);
    }
    return result;
}

/* A (nucleotide) sequence and its reverse complement. */
class Read {
public:
    const std::string& seq;
    const std::string rc;

    Read(const std::string& s)
      : seq(s)
      , rc(reverse_complement(s))
    {
    }

    std::string::size_type size() const {
        return seq.size();
    }
};

#endif
