#include "iowrap.hpp"
#include "exceptions.hpp"

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
