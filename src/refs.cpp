#include "refs.hpp"
#include <vector>
#include <fstream>
#include <cstdio>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include <cstring>
#include <zlib.h>

/* Convert string to uppercase in-place */
void to_uppercase(std::string& s) {
        std::transform(s.begin(), s.end(), s.begin(),
        [](unsigned char c) {
            return c & ~32;
        });
    }

namespace {

class GZFile {
public:
    GZFile(const std::string& filename, const char* mode) {
        file_ = gzopen(filename.c_str(), mode);
        if (!file_) {
            throw InvalidFasta("Could not open file " + filename);
        }
    }
    ~GZFile() {
        gzclose(file_);
    }
    int getc() {
        return gzgetc(file_);
    }
    int read(void* buf, unsigned int len) {
        return gzread(file_, buf, len);
    }

private:
    gzFile file_;
};

void fill_buffer(GZFile& file, char* buffer, size_t& buffer_size, size_t& buffer_pos, size_t max_buffer_size) {
    buffer_pos = 0;
    buffer_size = file.read(buffer, max_buffer_size);
}

std::string read_seq_name(GZFile& gf, char* buffer, size_t& buffer_size, const size_t max_buffer_size, size_t& buffer_pos) {
    char* eol = static_cast<char*>(std::memchr(buffer + buffer_pos, '\n', buffer_size - buffer_pos));
    if (eol) {
        // Normal (simple) case. Given the large size of the buffer, we expect this to be the case most of the time.
        char* space = static_cast<char*>(std::memchr(buffer + buffer_pos, ' ', eol - buffer + buffer_pos));
        char* name_end = space ? space : eol;
        std::string name{buffer + buffer_pos, name_end};
        buffer_pos = eol - buffer + 1;
        return name;
    }

    // We need to read more data, but first we save the data we already have
    std::string name;
    char* space = static_cast<char*>(std::memchr(buffer + buffer_pos, ' ', buffer_size - buffer_pos));
    bool found_space = false;
    if (space) {
        name.assign(buffer + buffer_pos, space);
        found_space = true;
    } else {
        name.assign(buffer + buffer_pos, buffer + buffer_size);
    }
    fill_buffer(gf, buffer, buffer_size, buffer_pos, max_buffer_size);
    if (buffer_size == 0) return name;
    eol = static_cast<char*>(std::memchr(buffer, '\n', buffer_size));
    if (!eol) { // No EOL found. This file is malformed
        throw InvalidFasta("Malformed FASTA file");
    }
    buffer_pos = eol - buffer + 1;
    if (!found_space) {
        space = static_cast<char*>(std::memchr(buffer, ' ', eol - buffer));
        if (!space) space = eol;
        name.append(buffer, space);
    }
    return name;
}

std::string read_sequence(GZFile& gf, char* buffer, size_t& buffer_size, const size_t max_buffer_size, size_t& buffer_pos) {
    if (buffer_pos == buffer_size) {
        fill_buffer(gf, buffer, buffer_size, buffer_pos, max_buffer_size);
        if (buffer_size == 0) return "";
    }
    std::string seq;
    while (buffer_size > 0 && buffer[buffer_pos] != '>') {
        char* eol = static_cast<char*>(std::memchr(buffer + buffer_pos, '\n', buffer_size - buffer_pos));
        if (eol) {
            seq.append(buffer + buffer_pos, eol);
            buffer_pos = eol - buffer + 1;
        } else {
            seq.append(buffer + buffer_pos, buffer + buffer_size);
            fill_buffer(gf, buffer, buffer_size, buffer_pos, max_buffer_size);
        }
    }
    to_uppercase(seq);
    return seq;
}
}


References References::from_fasta(const std::string& filename) {
    std::vector<std::string> sequences;
    ref_names names;

    GZFile gf(filename, "r");

    std::string name, seq;
    char buffer[128 * 1024];
    size_t buffer_size = 0;
    size_t buffer_pos = 0;
    buffer_size = gf.read(buffer, sizeof(buffer));
    if (buffer_size <= 0) {
        throw InvalidFasta("Could not read from file " + filename);
    }
    if (buffer[0] != '>') {
        throw InvalidFasta("File " + filename + " does not start with >");
    }

    do {
        // The read_sequence function should ensure that the buffer is not empty and at a '>' character
        assert(buffer_pos < buffer_size && buffer[buffer_pos] == '>');
        ++buffer_pos; // Skip the >
        name = read_seq_name(gf, buffer, buffer_size, sizeof(buffer), buffer_pos);
        if (name.empty()) {
            break;
        }
        seq = read_sequence(gf, buffer, buffer_size, sizeof(buffer), buffer_pos);
        if (!seq.empty()) {
            sequences.push_back(std::move(seq));
            names.push_back(std::move(name));
        }
    } while (buffer_size > 0);

    return References(std::move(sequences), std::move(names));
}

void References::add(std::string&& name, std::string&& sequence) {
    names.push_back(name);
    sequences.push_back(sequence);
    lengths.push_back(sequence.size());
}
