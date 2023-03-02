#ifndef STROBEALIGN_CIGAR_HPP
#define STROBEALIGN_CIGAR_HPP

#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include <cassert>


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

    explicit Cigar(std::vector<uint32_t> ops) : m_ops(std::move(ops)) { }

    Cigar(Cigar& other) : m_ops(other.m_ops) { }

    explicit Cigar(uint32_t* ops, size_t n) {
        m_ops.assign(ops, ops + n);
    }

    explicit Cigar(const std::string& cig);

    explicit Cigar(Cigar&& other) noexcept {
        *this = std::move(other);
    }

    Cigar& operator=(Cigar&& other) {
        if (this != &other) {
            m_ops = std::move(other.m_ops);
        }
        return *this;
    }

    bool empty() const { return m_ops.empty(); }

    void push(uint8_t op, int len) {
        assert(op < 16);
        if (m_ops.empty() || (m_ops.back() & 0xf) != op) {
            m_ops.push_back(len << 4 | op);
        } else {
            m_ops.back() += len << 4;
        }
    }

    void operator+=(const Cigar& other) {
        for (auto op_len : other.m_ops) {
            push(op_len & 0xf, op_len >> 4);
        }
    }

    /* This works only if I, D, X, = are the only operations used */
    int edit_distance() const {
        auto dist = 0;
        for (auto op_len : m_ops) {
            auto op = op_len & 0xf;
            auto len = op_len >> 4;
            if (op == CIGAR_INS || op == CIGAR_DEL || op == CIGAR_X) {
                dist += len;
            }
        }
        return dist;
    }

    void reverse() {
        std::reverse(m_ops.begin(), m_ops.end());
    }

    /* Return a new Cigar that uses =/X instead of M */
    Cigar to_eqx(const std::string& query, const std::string& ref) const;

    std::string to_string() const;

    std::vector<uint32_t> m_ops;
};

std::string compress_cigar(const std::string& ops);

#endif
