#ifndef STROBEALIGN_READLEN_HPP
#define STROBEALIGN_READLEN_HPP

#include "kseq++/kseq++.hpp"
#include <zlib.h>
#include "pc.hpp"

uint64_t estimate_read_length(InputBuffer& input_buffer);

#endif
