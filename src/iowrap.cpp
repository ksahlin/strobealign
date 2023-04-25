#include "iowrap.hpp"
#include "exceptions.hpp"
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <fcntl.h>
#include <system_error>
#include <cstring>

void GeneralIO::open(const std::string& filename)
{
    if(filename != "") {
        file = gzopen(filename.c_str(), "r");
        if(file == nullptr) {
            throw InvalidFile("Could not open FASTQ file: " + filename);
        }
    }
}

int64_t GeneralIO::read(void* buffer, size_t length)
{
    return gzread(file, buffer, length);
}

void RawIO::open(const std::string& filename)
{
    fd = ::open(filename.c_str(), 0);
    if(fd < 0) {
        throw InvalidFile("Could not open FASTQ file: " + filename);
    }
}

void RawIO::close()
{
    if(fd != -1) {
        ::close(fd);
    }
    fd = -1;
}

int64_t RawIO::read(void* buffer, size_t length)
{
    ssize_t ret = ::read(fd, buffer, length);
    return ret;
}

void IsalIO::initialize()
{
    isal_inflate_init(&state);
    state.crc_flag = ISAL_GZIP_NO_HDR_VER;

    isal_gzip_header_init(&gz_hdr);
}

void IsalIO::open(const std::string& filename)
{
    fd = ::open(filename.c_str(), 0);
    if(fd < 0) {
        throw InvalidFile("Could not open FASTQ file: " + filename);
    }

    struct stat _stat;
    if(fstat(fd, &_stat) < 0) {
        throw std::system_error(errno, std::generic_category(), filename);
    }
    filesize = _stat.st_size;

    mmap_mem = mmap(NULL, filesize, PROT_READ, MAP_SHARED, fd, 0);
    if(mmap_mem == MAP_FAILED) {
        mmap_mem = NULL;
        if(errno == ENODEV) {
            throw std::system_error(errno, std::generic_category(), "mmap is not supported on this file: " + filename);
        }
        if(errno == ENOMEM) {
            throw std::system_error(errno, std::generic_category(), "There not enough memory to open file: " + filename);
        } else {
            throw std::system_error(errno, std::generic_category(), "mmap failed to open file: " + filename);
        }
    }
    mmap_size       = filesize;
    compressed_data = reinterpret_cast<uint8_t*>(mmap_mem);
    compressed_size = mmap_size;

    // decompress gz header
    state.next_in  = compressed_data;
    state.avail_in = std::min(decompress_chunk_size, compressed_size);
    auto pre       = state.avail_in;

    int ret = isal_read_gzip_header(&state, &gz_hdr);
    if(ret != ISAL_DECOMP_OK) {
        throw std::runtime_error("Invalid gzip header found");
    }
    size_t processed_size = pre - state.avail_in;
    compressed_data += processed_size;
    compressed_size -= processed_size;

    thread_reader = std::thread(&IsalIO::decompress, this, std::min(decompress_chunk_size, compressed_size));
}

void IsalIO::close()
{
    if(thread_reader.joinable()) thread_reader.join();

    if(mmap_mem != nullptr) {
        munmap(mmap_mem, mmap_size);
    }
    mmap_mem  = nullptr;
    mmap_size = 0;

    if(fd != -1)
        ::close(fd);
    fd = -1;
}

int64_t IsalIO::read(void* buffer, size_t length)
{
    size_t actual_count = 0;
    while(length > 0) {
        size_t size_from_data = std::min(length, uncompressed_data.size() - uncompressed_data_copied);
        memcpy(buffer, uncompressed_data.data() + uncompressed_data_copied, size_from_data);
        buffer = (uint8_t*)buffer + size_from_data;
        length -= size_from_data;
        uncompressed_data_copied += size_from_data;
        actual_count += size_from_data;

        if(uncompressed_data_copied == uncompressed_data.size()) {
            if(thread_reader.joinable()) {
                thread_reader.join();
                std::swap(uncompressed_data, uncompressed_data_work);
                uncompressed_data_copied = 0;
                if(uncompressed_data.size() == 0) {
                    break;
                }
            }
            thread_reader = std::thread(&IsalIO::decompress, this, std::min(decompress_chunk_size, compressed_size));
        }
    }
    return actual_count;
}

void IsalIO::decompress(size_t count)
{
    uncompressed_data_work.resize(std::max(count, previous_member_size));
    uint8_t* ptr = uncompressed_data_work.data();

    while(ptr - uncompressed_data_work.data() < count) {
        if(compressed_size == 0) break;

        size_t actual_input_size     = std::min(decompress_chunk_size, compressed_size);
        size_t output_size_available = uncompressed_data_work.data() + uncompressed_data_work.size() - ptr;

        state.next_in  = compressed_data;
        state.avail_in = actual_input_size;

        state.next_out  = ptr;
        state.avail_out = output_size_available;

        auto ret = isal_inflate(&state);
        if(ret != ISAL_DECOMP_OK) {
            throw std::runtime_error("Error encountered whikle decompressing");
            exit(-1);
        }

        size_t processed = actual_input_size - state.avail_in;
        size_t decomp    = (size_t)(uintptr_t)(state.next_out - ptr);

        compressed_data += processed;
        compressed_size -= processed;
        ptr += decomp;
    }
    uncompressed_data_work.resize(ptr - uncompressed_data_work.data());
    previous_member_size = uncompressed_data_work.size();
}