#ifndef REFS_HPP
#define REFS_HPP

#include <cstdint>
#include <string>
#include <stdexcept>
#include <numeric>
#include <vector>
#include "exceptions.hpp"

class References {
    typedef std::vector<unsigned int> ref_lengths;
    typedef std::vector<std::string> ref_names;

public:
    References() { }
    References(
        std::vector<std::string> sequences_,
        ref_names names_
    ) : sequences(std::move(sequences_)), names(std::move(names_)) {

        if (sequences.size() != names.size()) {
            throw std::invalid_argument("lengths do not match");
        }
        lengths.reserve(sequences.size());
        for (auto& seq : sequences) {
            lengths.push_back(seq.size());
        }
        _total_length = std::accumulate(this->lengths.begin(), this->lengths.end(), (size_t)0);
    }

    void add(std::string&& name, std::string&& sequence);

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

void to_uppercase(std::string& s);

#endif
