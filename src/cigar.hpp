#ifndef STROBEALIGN_CIGAR_HPP
#define STROBEALIGN_CIGAR_HPP

#include <cstdint>
#include <string>
#include <vector>

enum CIGAR {
    CIGAR_MATCH = 0,
    CIGAR_INS = 1,
    CIGAR_DEL = 2,
    CIGAR_N_SKIP = 3,
    CIGAR_SOFTCLIP = 4,
    CIGAR_HARDCLIP = 5,
    CIGAR_PAD = 6,
    CIGAR_EQ = 7,
    CIGAR_X = 8
};

class Cigar {
public:
    explicit Cigar() { }
    explicit Cigar(uint32_t* ops, size_t n) {
        m_ops.assign(ops, ops + n);
    }

    explicit Cigar(const std::string& cig);
    void push(uint8_t op, int len);

    void operator+=(const Cigar& other) {
        for (auto op_len : other.m_ops) {
            push(op_len & 0xf, op_len >> 4);
        }
    }

    /* Return a new Cigar that uses =/X instead of M */
    Cigar to_eqx(const std::string& query, const std::string& ref) const;

    std::string to_string() const;

    std::vector<uint32_t> m_ops;
};

std::string compress_cigar(const std::string& ops);

#endif
