#include "fastq.hpp"

namespace {
    bool check_ext(const std::string& filename, const std::string& target_ext)
    {
        auto ext_pos = filename.find_last_of(".");
        if(ext_pos == std::string::npos) {
            return false;
        }
        auto ext = filename.substr(ext_pos, filename.size() - ext_pos);
        if(ext == target_ext) {
            return true;
        }

        return false;
    }

    inline bool is_gzip(const std::string& filename)
    {
        return check_ext(filename, ".gz");
    }

    inline bool is_raw(const std::string& filename)
    {
        return check_ext(filename, ".fq") || check_ext(filename, ".fastq");
    }

    std::unique_ptr<Reader> make_reader(const std::string& filename)
    {
        std::unique_ptr<Reader> io;
        if(is_gzip(filename)) {
            io = std::make_unique<IsalGzipReader>(filename);
        } else if(is_raw(filename)) {
            io = std::make_unique<UncompressedReader>(filename);
        } else {
            io = std::make_unique<GzipReader>(filename);
        }
        return io;
    }
};

RewindableFile::RewindableFile(const std::string& filename)
    : reader(nullptr),
    rewindable(true),
    stream_(klibpp::make_ikstream(this, rewind_read, 16384)) {

    reader = make_reader(filename);

    stream_ = klibpp::make_ikstream(this, rewind_read, 16384);
}

RewindableFile::~RewindableFile() {
}

void RewindableFile::rewind() {
    if (!rewindable) {
        throw std::runtime_error("Cannot rewind non-rewindable file");
    }
    rewindable = false;
    stream_ = klibpp::make_ikstream(this, rewind_read, 16384);
}

int RewindableFile::read(void* buffer, const int length) {
    if (reader == nullptr) {
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
    const auto bytes_read = reader->read(buffer, length);
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
