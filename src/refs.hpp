#ifndef FASTA_HPP
#define FASTA_HPP

#include <cstdint>
#include <string>
#include <stdexcept>
#include <numeric>
#include <vector>

typedef std::vector<unsigned int> ref_lengths;
typedef std::vector<std::string> ref_names;

class References {
public:
    References() { }
    References(
        std::vector<std::string>&& sequences,
        ref_names&& names,
        ref_lengths&& lengths
    ) : sequences(std::move(sequences)), names(std::move(names)), lengths(std::move(lengths)) {

        _total_length = std::accumulate(this->lengths.begin(), this->lengths.end(), (size_t)0);
        if (sequences.size() != names.size() || names.size() != lengths.size()) {
            throw std::invalid_argument("lengths do not match");
        }
    }

    static References from_fasta(const std::string& filename);

    size_t size() const {
        return sequences.size();
    }

    size_t total_length() const {
        return _total_length;
    }

    std::vector<std::string> sequences;
    ref_names names;
    ref_lengths lengths;
private:
    size_t _total_length;
};

#endif
