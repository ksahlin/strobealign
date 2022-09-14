#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "refs.hpp"

TEST_CASE("References::from_fasta") {
    auto references = References::from_fasta("tests/phix.fasta");
    CHECK(references.names.size() == 1);
    CHECK(references.lengths.size() == 1);
    CHECK(references.sequences.size() == 1);
    CHECK(references.names[0] == "NC_001422.1");
    CHECK(references.lengths[0] == 5386);
    CHECK(references.total_length() == 5386);
}
