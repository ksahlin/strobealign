#ifndef READLEN_HPP
#define READLEN_HPP

#include "kseq++.hpp"
#include <zlib.h>

int estimate_read_length(const std::string& filename1, const std::string& filename2);

#endif
