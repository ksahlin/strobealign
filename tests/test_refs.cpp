#include <fstream>
#include <zlib.h>

#include "doctest.h"
#include "refs.hpp"

TEST_CASE("References::add") {
    References references;
    references.add(std::string("thename"), std::string("ACGT"));

    CHECK(references.names.size() == 1);
    CHECK(references.names[0] == "thename");
    CHECK(references.sequences[0] == "ACGT");
    CHECK(references.lengths[0] == 4);
}

TEST_CASE("References::from_fasta") {
    auto references = References::from_fasta("tests/phix.fasta");
    CHECK(references.names.size() == 1);
    CHECK(references.names[0] == "NC_001422.1");
    CHECK(references.sequences.size() == 1);
    CHECK(references.lengths.size() == 1);
    CHECK(references.lengths[0] == 5386);
    CHECK(references.total_length() == 5386);
}

TEST_CASE("References::from_fasta parse error") {
    REQUIRE_THROWS_AS(References::from_fasta("tests/phix.1.fastq"), InvalidFasta);
}

TEST_CASE("Reference FASTA not found") {
    REQUIRE_THROWS_AS(References::from_fasta("does-not-exist.fasta"), InvalidFasta);
}

TEST_CASE("Reference uppercase") {
    {
        std::ofstream ofs("tmpref.fasta");
        ofs
            << ">ref1\n"
            << "acgt\n\n"
            << ">ref2\n"
            << "aacc\ngg\n\ntt\n"
            << ">empty\n"
            << ">empty_at_end_of_file";
    }
    auto refs = References::from_fasta("tmpref.fasta");
    std::remove("tmpref.fasta");
    CHECK(refs.sequences.size() == 2);
    CHECK(refs.sequences[0].size() == 4);
    CHECK(refs.sequences[0] == "ACGT");
    CHECK(refs.sequences[1].size() == 8);
    CHECK(refs.sequences[1] == "AACCGGTT");
    CHECK(refs.names.size() == 2);
    CHECK(refs.lengths.size() == 2);
}

// Copy & pasted from the above, but with a gzipped file
TEST_CASE("Reference gzipped") {
    {
        gzFile writer = gzopen("tmpref.fasta.gz", "w");
        gzputs(writer, ">ref1\n");
        gzputs(writer, "acgt\n\n");
        gzputs(writer, ">ref2\n");
        gzputs(writer, "aacc\ngg\n\ntt\n");
        gzputs(writer, ">empty\n");
        gzputs(writer, ">empty_at_end_of_file");
        gzclose(writer);
    }
    auto refs = References::from_fasta("tmpref.fasta.gz");
    std::remove("tmpref.fasta");
    CHECK(refs.sequences.size() == 2);
    CHECK(refs.sequences[0].size() == 4);
    CHECK(refs.sequences[0] == "ACGT");
    CHECK(refs.sequences[1].size() == 8);
    CHECK(refs.sequences[1] == "AACCGGTT");
    CHECK(refs.names.size() == 2);
    CHECK(refs.lengths.size() == 2);
}
