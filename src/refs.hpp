#ifndef FASTA_HPP
#define FASTA_HPP

#include <cstdint>
#include <string>
#include <stdexcept>
#include <numeric>
#include <vector>


class References {
public:
    References() { }
    References(
        std::vector<std::string>&& sequences,
        std::vector<std::string>&& names,
        std::vector<unsigned int>&& lengths
    ) : sequences(std::move(sequences)), names(std::move(names)), lengths(std::move(lengths)) {

        _total_length = std::accumulate(this->lengths.begin(), this->lengths.end(), (size_t)0);
        if (sequences.size() != names.size() || names.size() != lengths.size()) {
            throw std::invalid_argument("lengths do not match");
        }
    }

    static References from_fasta(std::string path);

    size_t size() const {
        return sequences.size();
    }

    size_t total_length() const {
        return _total_length;
    }

    std::vector<std::string> sequences;
    std::vector<std::string> names;
    std::vector<unsigned int> lengths;
private:
    size_t _total_length;
};

#endif
