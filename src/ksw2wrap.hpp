#ifndef STROBEALIGN_KSW2WRAP_HPP
#define STROBEALIGN_KSW2WRAP_HPP

#include <string>
#include <cstdint>
#include "aligner.hpp"

aln_info ksw_extend(const std::string& query, const std::string& ref, int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend, bool reverse_cigar = false);

#endif
