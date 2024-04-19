#ifndef STROBEALIGN_IOWRAP_HPP
#define STROBEALIGN_IOWRAP_HPP

#include <cstdint>
#include <string>

#include <isa-l/igzip_lib.h>
#include <zlib.h>
#include <vector>

#include <thread>

class Reader {
   public:
    Reader(const std::string&) { }

    virtual ~Reader() { }

    virtual int64_t read(void* buffer, size_t length) = 0;
    [[maybe_unused]] virtual std::string name() const = 0;

   protected:
    virtual void open(const std::string& filename) = 0;
};

class GzipReader : public Reader {
   public:
    GzipReader(const std::string& filename) : Reader(filename), file() { open(filename); }

    virtual ~GzipReader() {
        if (file) {
            gzclose(file);
        }
    }

    int64_t read(void* buffer, size_t length) override;
    [[maybe_unused]] std::string name() const override { return "GzipReader"; }

   private:
    gzFile file;
    void open(const std::string& filename) override;
};

class UncompressedReader : public Reader {
   public:
    UncompressedReader(const std::string& filename)
        : Reader(filename)
        , fd(-1)
        , preload_size(64ull * 1024 * 1024)
        , read_buffer()
        , read_buffer_work()
        , read_buffer_copied(0)
        , thread_reader() {
        open(filename);
    }

    virtual ~UncompressedReader() {
        if (fd != -1) {
            close();
        }
    }

    int64_t read(void* buffer, size_t length) override;
    [[maybe_unused]] std::string name() const override { return "UncompressReader"; }

   private:
    int fd;

    void preload(size_t size);

    size_t preload_size;

    std::vector<uint8_t> read_buffer;
    std::vector<uint8_t> read_buffer_work;
    size_t read_buffer_copied;

    std::thread thread_reader;

    void open(const std::string& filename) override;
    void close();
};

class IsalGzipReader : public Reader {
   public:
    IsalGzipReader(const std::string& filename)
        : Reader(filename)
        , fd(-1)
        , mmap_mem(nullptr)
        , filesize(-1)
        , mmap_size(-1)
        , uncompressed_data()
        , uncompressed_data_work()
        , uncompressed_data_copied(0)
        , compressed_data(nullptr)
        , compressed_size(0)
        , decompress_chunk_size(2ull * 1024 * 1024)
        , previous_member_size(8ull * 1024 * 1024)
        , thread_reader() {
        initialize();
        open(filename);
    }

    virtual ~IsalGzipReader() {
        if (fd != -1) {
            close();
        }
    }

    int64_t read(void* buffer, size_t length) override;
    [[maybe_unused]] std::string name() const override { return "IsalGzipReader"; }

   private:
    int fd;
    void* mmap_mem;
    ssize_t filesize;
    ssize_t mmap_size;
    std::vector<uint8_t> uncompressed_data;
    std::vector<uint8_t> uncompressed_data_work;
    size_t uncompressed_data_copied;

    uint8_t* compressed_data;
    size_t compressed_size;

    size_t decompress_chunk_size;
    size_t previous_member_size;

    std::thread thread_reader;

    inflate_state state;
    isal_gzip_header gz_hdr;

    void initialize();
    void open(const std::string& filename) override;
    void close();

    void decompress(size_t count);
};

#endif
