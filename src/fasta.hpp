#ifndef FASTA_HPP
#define FASTA_HPP

#include <cstdint>
#include <string>
#include "aln.hpp"  // idx_to_acc

uint64_t read_references(std::vector<std::string> &seqs, std::vector<unsigned int> &lengths, idx_to_acc &acc_map, std::string fn);

#endif
