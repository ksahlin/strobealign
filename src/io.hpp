#ifndef IO_HPP
#define IO_HPP

#include <iostream>
#include <vector>

void write_int_to_ostream(std::ostream& os, int value);
int read_int_from_istream(std::istream& is);

/* Write a vector to an output stream, preceded by its length */
template <typename T>
void write_vector(std::ostream& os, const std::vector<T>& v) {
    auto size = uint64_t(v.size());
    os.write(reinterpret_cast<char*>(&size), sizeof(size));
    os.write(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(T));
}

/* Read a vector written by write_vector */
template <typename T>
void read_vector(std::istream& is, std::vector<T>& v) {
    uint64_t size;
    v.clear();
    is.read(reinterpret_cast<char*>(&size), sizeof(size));
    v.resize(size);
    is.read(reinterpret_cast<char*>(v.data()), size * sizeof(T));
}

#endif
