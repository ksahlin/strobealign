#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "refs.hpp"
#include "exceptions.hpp"
#include "readlen.hpp"
#include "index.hpp"
#include "fastq.hpp"
#include "sam.hpp"


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

TEST_CASE("Reads file missing") {
    std::string filename("does-not-exist.fastq");
    REQUIRE_THROWS_AS(open_fastq(filename), InvalidFile);
}


TEST_CASE("unmapped SAM record") {
    klibpp::KSeq kseq;
    kseq.name = "read1";
    kseq.seq = "ACGT";
    kseq.qual = ">#BB";
    std::string sam_string;
    Sam sam(sam_string, References());

    sam.add_unmapped(kseq);

    CHECK(sam_string == "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t>#BB\n");
}
