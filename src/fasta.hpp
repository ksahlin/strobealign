#ifndef FASTA_HPP
#define FASTA_HPP

#include <cstdint>
#include <string>
#include "aln.hpp"  // idx_to_acc

uint64_t read_references(std::vector<std::string> &seqs, std::vector<unsigned int> &lengths, idx_to_acc &acc_map, std::string fn);


class References {
public:
    References(
        std::vector<std::string>&& sequences,
        std::vector<std::string>&& names,
        std::vector<unsigned int>&& lengths
    ) : sequences(std::move(sequences)), names(std::move(names)), lengths(std::move(lengths)) { }

    static References from_fasta(std::string path);

    const std::vector<std::string> sequences;
    const std::vector<std::string> names;
    const std::vector<unsigned int> lengths;
};

#endif
