#include "fastq.hpp"

input_stream_t open_fastq(std::string& filename) {
    gzFile fp = gzopen(filename.c_str(), "r");
    if (fp == nullptr) {
        throw InvalidFile("Could not open FASTQ file");
    }
    // 16384 is the default buffer size. We need to pass it to be able to pass gzclose.
    return klibpp::make_ikstream(fp, gzread, 16384, gzclose);
}
