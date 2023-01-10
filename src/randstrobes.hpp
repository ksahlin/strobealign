#ifndef RANDSTROBES_H
#define RANDSTROBES_H

#include <vector>
#include <string>
#include <tuple>
#include <inttypes.h>

// only used during index generation
struct MersIndexEntry {
    uint64_t hash;
    uint32_t position;
    int32_t packed; // packed representation of ref_index and strobe offset

    bool operator< (const MersIndexEntry& other) const {
        return hash < other.hash;
    }
};

typedef std::vector<MersIndexEntry> ind_mers_vector;


struct QueryMer {
    uint64_t hash;
    unsigned int start;
    unsigned int end;
    bool is_reverse;
};


typedef std::vector<QueryMer> mers_vector_read;

void randstrobes_reference(ind_mers_vector& flat_vector, int k, int w_min, int w_max, const std::string &seq, int ref_index, int s, int t, uint64_t q, int max_dist);
mers_vector_read randstrobes_query(             int k, int w_min, int w_max, const std::string &seq,                int s, int t, uint64_t q, int max_dist);

#endif
