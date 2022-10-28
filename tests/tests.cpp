#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "refs.hpp"
#include "exceptions.hpp"
#include "readlen.hpp"
#include "index.hpp"


TEST_CASE("References::from_fasta") {
    auto references = References::from_fasta("tests/phix.fasta");
    CHECK(references.names.size() == 1);
    CHECK(references.lengths.size() == 1);
    CHECK(references.sequences.size() == 1);
    CHECK(references.names[0] == "NC_001422.1");
    CHECK(references.lengths[0] == 5386);
    CHECK(references.total_length() == 5386);
}

TEST_CASE("References::from_fasta parse error") {
    REQUIRE_THROWS_AS(References::from_fasta("tests/phix.1.fastq"), InvalidFasta);
}

TEST_CASE("Reference FASTA not found") {
    REQUIRE_THROWS_AS(References::from_fasta("does-not-exist.fasta"), InvalidFasta);
}

TEST_CASE("estimate_read_length") {
    CHECK(estimate_read_length("tests/phix.1.fastq", "") == 296);
    CHECK(estimate_read_length("tests/phix.1.fastq", "tests/phix.2.fastq") == 296);
}

TEST_CASE("IndexParameters==") {
    IndexParameters a = IndexParameters::from_read_length(150, 8);
    IndexParameters b = IndexParameters::from_read_length(150, 8);
    CHECK(a == b);
}

TEST_CASE("sti file same parameters") {
    auto references = References::from_fasta("tests/phix.fasta");
    auto parameters = IndexParameters::from_read_length(300, 8);
    StrobemerIndex index(references, parameters);
    index.populate(0.0002);
    index.write("tmpindex.sti");

    auto other_parameters = IndexParameters::from_read_length(30, 8);
    StrobemerIndex other_index(references, other_parameters);

    REQUIRE_THROWS_AS(other_index.read("tmpindex.sti"), InvalidIndexFile);
    std::remove("tmpindex.sti");
}
