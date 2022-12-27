#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "refs.hpp"
#include "exceptions.hpp"
#include "readlen.hpp"
#include "index.hpp"
#include "fastq.hpp"
#include "aln.hpp"


TEST_CASE("estimate_read_length") {
    CHECK(estimate_read_length("tests/phix.1.fastq", "") == 289);
    CHECK(estimate_read_length("tests/phix.1.fastq", "tests/phix.2.fastq") == 289);
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

TEST_CASE("Reads file missing") {
    std::string filename("does-not-exist.fastq");
    REQUIRE_THROWS_AS(open_fastq(filename), InvalidFile);
}

TEST_CASE("has_shared_substring") {
    std::string ref{"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"};
    std::string read{"GGGGGGGGGGGGGGGGG"};
    CHECK(!has_shared_substring(read, ref, 20));
}

TEST_CASE("read_/write_vector") {
    std::string filename{"tmp-vector"};
    const std::vector<int> expected{2, 3, 5, 7, 11};
    {
        std::ofstream ofs{filename, std::ios::binary};
        write_vector(ofs, expected);
    }
    std::vector<int> y;
    {
        std::ifstream ifs{filename, std::ios::binary};
        read_vector(ifs, y);
    }
    CHECK(y == expected);
}
