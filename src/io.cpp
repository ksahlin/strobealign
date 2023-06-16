#include "io.hpp"

void write_int_to_ostream(std::ostream& os, int32_t value) {
    int32_t val;
    val = value;
    os.write(reinterpret_cast<const char*>(&val), sizeof(val));
}

int32_t read_int_from_istream(std::istream& is) {
    int32_t val;
    is.read(reinterpret_cast<char*>(&val), sizeof(val));
    return val;
}
