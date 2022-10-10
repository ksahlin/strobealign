#ifndef RANDSTROBES_H
#define RANDSTROBES_H

#include <vector>
#include <string>
#include <tuple>
#include <inttypes.h>

typedef std::vector<std::tuple<uint64_t, uint32_t, int32_t >> ind_mers_vector; //only used during index generation
typedef std::vector<std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, bool>> mers_vector_read;

void seq_to_randstrobes2(ind_mers_vector& flat_vector, int n, int k, int w_min, int w_max, const std::string &seq, int ref_index,          int s, int t, uint64_t q, int max_dist);
mers_vector_read seq_to_randstrobes2_read(             int n, int k, int w_min, int w_max, const std::string &seq, unsigned int ref_index, int s, int t, uint64_t q, int max_dist);

#endif
