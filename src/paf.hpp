#ifndef STROBEALIGN_PAF_HPP
#define STROBEALIGN_PAF_HPP

#include <string>
#include "refs.hpp"
#include "chain.hpp"

void output_hits_paf_PE(
    std::string &paf_output, const Chain &n, const std::string &query_name, const References& references, int read_len, uint8_t mapq
);

void output_hits_paf(
    std::string &paf_output, const std::vector<Chain> &all_nams, const std::string& query_name, const References& references, int read_len, uint8_t mapq
);

#endif
