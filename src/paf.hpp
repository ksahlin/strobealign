#ifndef STROBEALIGN_PAF_HPP
#define STROBEALIGN_PAF_HPP

#include <string>
#include "refs.hpp"
#include "nam.hpp"

void output_hits_paf_PE(
    std::string &paf_output, const Nam &n, const std::string &query_name, const References& references, int k, int read_len
);

void output_hits_paf(
    std::string &paf_output, const std::vector<Nam> &all_nams, const std::string& query_name, const References& references, int k, int read_len
);

#endif
