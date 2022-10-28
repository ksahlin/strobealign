#ifndef FASTQ_HPP
#define FASTQ_HPP

#include <zlib.h>
#include <string>

#include "exceptions.hpp"
#include "kseq++.hpp"

typedef klibpp::KStream<gzFile_s*, int (*)(gzFile_s*, void*, unsigned int), klibpp::mode::In_> input_stream_t;

input_stream_t open_fastq(std::string& filename);

#endif
