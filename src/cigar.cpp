#include <cassert>
#include <sstream>
#include "cigar.hpp"


Cigar Cigar::to_eqx(const std::string& query, const std::string& ref) const {
    size_t i = 0, j = 0;
    Cigar cigar;
    for (auto op_len : m_ops) {
        auto op = op_len & 0xf;
        auto len = op_len >> 4;
        if (op == CIGAR_MATCH) {
            for (size_t u = 0; u < len; ++u) {
                if (query[i] == ref[j]) {
                    cigar.push(CIGAR_EQ, 1);
                } else {
                    cigar.push(CIGAR_X, 1);
                }
                i++;
                j++;
            }
        } else if (op == CIGAR_INS) {
            cigar.push(op, len);
            i += len;
        } else if (op == CIGAR_DEL) {
            cigar.push(op, len);
            j += len;
        }
    }
    return cigar;
}

std::string Cigar::to_string() const {
    std::stringstream s;
    for (auto op_len : m_ops) {
        s << (op_len >> 4) << "MIDNSHP=X"[op_len & 0xf];
    }
    return s.str();
}

Cigar::Cigar(const std::string& cig) {
    int number = -1;  // -1 means "not given"
    for (auto c : cig) {
        if (isdigit(c)) {
            if (number == -1) {
                number = c - '0';
            } else {
                number = number * 10 + (c - '0');
            }
        } else {
            int op;
            switch(c) {
                case ' ': continue;
                case 'M': op = CIGAR_MATCH; break;
                case 'I': op = CIGAR_INS; break;
                case 'D': op = CIGAR_DEL; break;
                case 'N': op = CIGAR_N_SKIP; break;
                case 'S': op = CIGAR_SOFTCLIP; break;
                case 'H': op = CIGAR_HARDCLIP; break;
                case 'P': op = CIGAR_PAD; break;
                case '=': op = CIGAR_EQ; break;
                case 'X': op = CIGAR_X; break;
                default: throw std::invalid_argument("Invalid CIGAR operator"); break;
            }
            if (number == -1) {
                push(op, 1);
            } else if (number > 0) {
                push(op, number);
                number = -1;
            }
        }
    }
    if (number != -1) {
        throw std::invalid_argument("CIGAR must not end with a number");
    }
}

std::string compress_cigar(const std::string& ops) {
    char prev = 0;
    int count = 0;
    std::stringstream cigar;
    bool first = true;
    for (auto op : ops) {
        if (!first && op != prev) {
            cigar << count << prev;
            count = 0;
        }
        count++;
        prev = op;
        first = false;
    }
    if (!first) {
        cigar << count << prev;
    }
    return cigar.str();
}
