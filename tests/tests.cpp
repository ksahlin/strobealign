#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "fasta.hpp"

TEST_CASE("read_references") {
    std::vector<std::string> seqs;
    std::vector<unsigned int> lengths;
    idx_to_acc reference_names;
    uint64_t n = read_references(seqs, lengths, reference_names, "tests/phix.fasta");

    CHECK(n == 5386);
    CHECK(reference_names.size() == 1);
    CHECK(lengths.size() == 1);
    CHECK(seqs.size() == 1);
    CHECK(reference_names[0] == "NC_001422.1");
    CHECK(lengths[0] == 5386);
}
