#include "fastq.hpp"

RewindableFile::RewindableFile(const std::string& filename)
    : file(nullptr),
    rewindable(true),
    stream_(klibpp::make_ikstream(this, rewind_read, 16384)) {
    if (filename != "") {
        file = gzopen(filename.c_str(), "r");
        if (file == nullptr) {
            throw InvalidFile("Could not open FASTQ file: " + filename);
        }
    }
    stream_ = klibpp::make_ikstream(this, rewind_read, 16384);
}

RewindableFile::~RewindableFile() {
    if (file) gzclose(file);
}

void RewindableFile::rewind() {
    if (!rewindable) {
        throw std::runtime_error("Cannot rewind non-rewindable file");
    }
    rewindable = false;
    stream_ = klibpp::make_ikstream(this, rewind_read, 16384);
}

int RewindableFile::read(void* buffer, const int length) {
    if (file == nullptr) {
        return -1;
    }
    if (length == 0) return 0;
    if (!rewindable && !saved_buffer.empty()) {
        const int can_read = std::min(length, static_cast<int>(saved_buffer[0].size()));
        std::memcpy(buffer, saved_buffer[0].data(), can_read);
        if (unsigned(can_read) == saved_buffer[0].size()) {
            saved_buffer.erase(saved_buffer.begin());
        } else {
            saved_buffer[0].erase(saved_buffer[0].begin(), saved_buffer[0].begin() + can_read);
        }
        return can_read +
            this->read(static_cast<char*>(buffer) + can_read, length - can_read);
    }
    const int bytes_read = gzread(file, buffer, length);
    if (bytes_read < 0) {
        throw std::runtime_error("Error reading FASTQ file");
    }
    if (rewindable) {
        saved_buffer.push_back(std::vector<unsigned char>(
                    static_cast<unsigned char*>(buffer),
                    static_cast<unsigned char*>(buffer) + bytes_read));
    }
    return bytes_read;
}

int rewind_read(RewindableFile* file, void* buffer, unsigned int length) {
    return file->read(buffer, length);
}

input_stream_t open_fastq(std::string& filename) {
    if (filename == "-") {
        filename = "/dev/stdin";
    }
    return std::unique_ptr<RewindableFile>(new RewindableFile(filename));
}
