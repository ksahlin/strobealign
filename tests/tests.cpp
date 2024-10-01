#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "refs.hpp"
#include "exceptions.hpp"
#include "readlen.hpp"
#include "index.hpp"
#include "fastq.hpp"
#include "tmpdir.hpp"
#include "io.hpp"
#include "revcomp.hpp"


TEST_CASE("estimate_read_length") {
    InputBuffer input_buffer1("tests/phix.1.fastq", "", 5, false);
    CHECK(estimate_read_length(input_buffer1) == 289);

    InputBuffer input_buffer2("tests/phix.1.fastq", "tests/phix.2.fastq", 5, false);
    CHECK(estimate_read_length(input_buffer2) == 287);
}

TEST_CASE("IndexParameters==") {
    IndexParameters a = IndexParameters::from_read_length(150);
    IndexParameters b = IndexParameters::from_read_length(150);
    CHECK(a == b);
}

TEST_CASE("IndexParameters identical for slight read-length differences") {
    IndexParameters a = IndexParameters::from_read_length(150);
    IndexParameters b = IndexParameters::from_read_length(151);
    IndexParameters c = IndexParameters::from_read_length(149);
    CHECK(a == b);
    CHECK(a == c);
}

TEST_CASE("Parameters in sti file do not match") {
    TemporaryDirectory tmp_dir;
    std::string ref_path = (tmp_dir.path() / "ref.fasta").string();
    std::filesystem::copy("tests/phix.fasta", ref_path);
    auto references = References::from_fasta(ref_path);
    auto parameters = IndexParameters::from_read_length(300);
    StrobemerIndex index(references, parameters);
    index.populate(0.0002, 1);
    std::string sti_path = (tmp_dir.path() / "index.sti").string();
    index.write(sti_path);

    auto other_parameters = IndexParameters::from_read_length(50);
    StrobemerIndex other_index(references, other_parameters);

    REQUIRE_THROWS_AS(other_index.read(sti_path), InvalidIndexFile);
}

TEST_CASE("Missing sti file") {
    TemporaryDirectory tmp_dir;
    auto references = References::from_fasta("tests/phix.fasta");
    auto parameters = IndexParameters::from_read_length(300);
    StrobemerIndex index(references, parameters);
    REQUIRE_THROWS_AS(index.read((tmp_dir.path() / "index.sti").string()), InvalidIndexFile);
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
    TemporaryDirectory tmp_dir;
    std::string filename = (tmp_dir.path() / "vector").string();
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

TEST_CASE("RandstrobeGenerator is sane") {
    auto references = References::from_fasta("tests/phix.fasta");
    auto& seq = references.sequences[0];
    auto parameters = IndexParameters::from_read_length(300);

    RandstrobeGenerator iter{seq, parameters.syncmer, parameters.randstrobe};

    Randstrobe randstrobe;
    while ((randstrobe = iter.next()) != iter.end()) {
        CHECK(randstrobe.hash > 0);
        CHECK(randstrobe.strobe1_pos >= 0);
        CHECK(randstrobe.strobe2_pos >= randstrobe.strobe2_pos);
        CHECK(randstrobe.strobe1_pos < seq.length());
        CHECK(randstrobe.strobe2_pos < seq.length());
    }
}

TEST_CASE("RandstrobeIterator is sane") {
    auto references = References::from_fasta("tests/phix.fasta");
    auto& seq = references.sequences[0];
    auto parameters = IndexParameters::from_read_length(300);

    auto syncmers = canonical_syncmers(seq, parameters.syncmer);
    RandstrobeIterator iter{syncmers, parameters.randstrobe};

    while (iter.has_next()) {
        auto randstrobe = iter.next();
        CHECK(randstrobe.hash > 0);
        CHECK(randstrobe.strobe1_pos >= 0);
        CHECK(randstrobe.strobe2_pos >= randstrobe.strobe2_pos);
        CHECK(randstrobe.strobe1_pos < seq.length());
        CHECK(randstrobe.strobe2_pos < seq.length());
    }
}

TEST_CASE("both randstrobes iterator implementations give same results") {
    auto references = References::from_fasta("tests/phix.fasta");
    auto& seq = references.sequences[0];
    auto parameters = IndexParameters::from_read_length(300);

    auto syncmers = canonical_syncmers(seq, parameters.syncmer);
    RandstrobeIterator iter1{syncmers, parameters.randstrobe};
    RandstrobeGenerator iter2{seq, parameters.syncmer, parameters.randstrobe};

    while (iter1.has_next()) {
        auto randstrobe1 = iter1.next();
        auto randstrobe2 = iter2.next();
        CHECK(randstrobe2 != iter2.end());
        CHECK(randstrobe1 == randstrobe2);
    }
    CHECK(iter2.next() == iter2.end());
}

// This tests an assumption that we need to hold for count_randstrobes()
TEST_CASE("syncmer and randstrobe iterators return (nearly) same no. of items") {
    auto references = References::from_fasta("tests/phix.fasta");
    auto& seq = references.sequences[0];
    auto parameters = IndexParameters::from_read_length(100);

    uint64_t randstrobe_count = 0;
    auto randstrobe_iter = RandstrobeGenerator(seq, parameters.syncmer, parameters.randstrobe);
    Randstrobe randstrobe;
    while ((randstrobe = randstrobe_iter.next()) != randstrobe_iter.end()) {
        randstrobe_count++;
    }

    uint64_t syncmer_count = 0;
    auto syncmer_iterator = SyncmerIterator(seq, parameters.syncmer);
    Syncmer syncmer;
    while (!(syncmer = syncmer_iterator.next()).is_end()) {
        syncmer_count++;
    }

    CHECK(randstrobe_count + parameters.randstrobe.w_min == syncmer_count);
}

TEST_CASE("reverse complement") {
    CHECK(reverse_complement("") == "");
    CHECK(reverse_complement("A") == "T");
    CHECK(reverse_complement("C") == "G");
    CHECK(reverse_complement("G") == "C");
    CHECK(reverse_complement("T") == "A");
    CHECK(reverse_complement("TG") == "CA");
    CHECK(reverse_complement("AC") == "GT");
    CHECK(reverse_complement("ACG") == "CGT");
    CHECK(reverse_complement("AACGT") == "ACGTT");
}
