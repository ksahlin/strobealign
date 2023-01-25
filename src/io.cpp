#include "io.hpp"

void write_int_to_ostream(std::ostream& os, int value) {
    int val;
    val = value;
    os.write(reinterpret_cast<const char*>(&val), sizeof(val));
}

int read_int_from_istream(std::istream& is) {
    int val;
    is.read(reinterpret_cast<char*>(&val), sizeof(val));
    return val;
}
