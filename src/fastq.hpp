#ifndef FASTQ_HPP
#define FASTQ_HPP

#include <zlib.h>
#include <string>

#include "exceptions.hpp"
#include "kseq++.hpp"

// File that can be rewound (once only!)
class RewindableFile {

public:
    typedef klibpp::KStream<RewindableFile*, int (*)(RewindableFile*, void*, unsigned int), klibpp::mode::In_> stream_type;
    // if filename == "", then the result is a null file (i.e., every read fails)
    explicit RewindableFile(const std::string& filename);
    ~RewindableFile();

    stream_type& stream() { return stream_; }
    int read(void* buffer, int length);

    // Reset to the beginning of the file. Can only be called once!
    void rewind();

protected:
    gzFile file;
    std::vector<std::vector<unsigned char>> saved_buffer;
    // if rewindable is false, the file cannot be rewound anymore and is consuming from saved_buffer (if it is not empty)
    bool rewindable;
    stream_type stream_;
};

int rewind_read(RewindableFile* file, void* buffer, unsigned int length);

typedef std::unique_ptr<RewindableFile> input_stream_t;

input_stream_t open_fastq(std::string& filename);

#endif
