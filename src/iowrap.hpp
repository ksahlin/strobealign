#ifndef IO_WRAP_HPP
#define IO_WRAP_HPP

#include <string>
#include <cstdint>

#include <zlib.h>

class AbstructIO
{
public:
    AbstructIO(const std::string&)
    {}

    virtual ~AbstructIO() {}

    virtual int64_t read(void* buffer, size_t length) = 0;

protected:
    virtual void open(const std::string& filename) = 0;
};


class GeneralIO : public AbstructIO
{
public:
    GeneralIO(const std::string& filename)
        : AbstructIO(filename)
    {
        open(filename);
    }

    virtual ~GeneralIO()
    {
        if(file) {
            gzclose(file);
        }
    }

    int64_t read(void* buffer, size_t length) override;

private:
    gzFile file;
    void open(const std::string& filename) override;
};

#endif
