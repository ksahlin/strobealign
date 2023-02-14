#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "refs.hpp"
#include "exceptions.hpp"
#include "readlen.hpp"
#include "index.hpp"
#include "fastq.hpp"
#include "aln.hpp"
#include "tmpdir.hpp"
#include "io.hpp"
#include "revcomp.hpp"
#include "aln.hpp"


TEST_CASE("estimate_read_length") {
    InputBuffer input_buffer1("tests/phix.1.fastq", "", 5, false);
    CHECK(estimate_read_length(input_buffer1) == 289);

    InputBuffer input_buffer2("tests/phix.1.fastq", "tests/phix.2.fastq", 5, false);
    CHECK(estimate_read_length(input_buffer2) == 289);
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
    std::string ref_path = tmp_dir.path() / "ref.fasta";
    std::filesystem::copy("tests/phix.fasta", ref_path);
    auto references = References::from_fasta(ref_path);
    auto parameters = IndexParameters::from_read_length(300);
    StrobemerIndex index(references, parameters);
    index.populate(0.0002, 1);
    std::string sti_path = tmp_dir.path() / "index.sti";
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
    REQUIRE_THROWS_AS(index.read(tmp_dir.path() / "index.sti"), InvalidIndexFile);
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
    std::string filename = tmp_dir.path() / "vector";
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

TEST_CASE("both randstrobes iterator implementations give same results") {
    auto references = References::from_fasta("tests/phix.fasta");
    auto& seq = references.sequences[0];
    auto parameters = IndexParameters::from_read_length(300, 8);

    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_coordinate;
    std::tie(string_hashes, pos_to_seq_coordinate) = make_string_to_hashvalues_open_syncmers_canonical(seq, parameters.k, parameters.s, parameters.t_syncmer);
    RandstrobeIterator iter1{string_hashes, pos_to_seq_coordinate, parameters.w_min, parameters.w_max, parameters.q, parameters.max_dist };
    RandstrobeIterator2 iter2(seq, parameters.k, parameters.s, parameters.t_syncmer, parameters.w_min, parameters.w_max, parameters.q, parameters.max_dist);

    while (iter1.has_next()) {
        auto randstrobe1 = iter1.next();
        auto randstrobe2 = iter2.next();
        CHECK(randstrobe2 != iter2.end());
        CHECK(randstrobe1 == randstrobe2);
    }
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

TEST_CASE("hamming_align") {
    int soft_left, soft_right;

    // empty sequences
    auto info = hamming_align(
        "", "",
        7, 5,
        soft_left, soft_right
    );
    CHECK(info.cigar == "");
    CHECK(info.ed == 0);
    CHECK(info.global_ed == 0);
    CHECK(info.length == 0);
    CHECK(info.sw_score == 0);
    CHECK(soft_left == 0);
    CHECK(soft_right == 0);
    CHECK(info.ref_offset == soft_left);

    info = hamming_align(
        "AAXGGG",
        "AAYGGG",
        1, 1,
        soft_left, soft_right
    );
    CHECK(info.cigar == "2=1X3=");
    CHECK(info.ed == 1);
    CHECK(info.global_ed == 1);
    CHECK(info.length == 6);
    CHECK(info.sw_score == 4);
    CHECK(soft_left == 0);
    CHECK(soft_right == 0);
    CHECK(info.ref_offset == soft_left);

    info = hamming_align(
        "AXGGG",
        "AYGGG",
        1, 4,
        soft_left, soft_right
    );
    CHECK(info.cigar == "2S3=");
    CHECK(info.ed == 0);
    CHECK(info.global_ed == 2);
    CHECK(info.length == 3);
    CHECK(info.sw_score == 3);
    CHECK(soft_left == 2);
    CHECK(soft_right == 0);
    CHECK(info.ref_offset == soft_left);

    info = hamming_align(
        "NAACCG",
        "TAACCG",
        3, 7,
        soft_left, soft_right
    );
    CHECK(info.cigar == "1S5=");
    CHECK(info.ed == 0);
    CHECK(info.global_ed == 1);
    CHECK(info.length == 5);
    CHECK(info.sw_score == 5 * 3);
    CHECK(soft_left == 1);
    CHECK(soft_right == 0);
    CHECK(info.ref_offset == soft_left);

    info = hamming_align(
        "AACCGN",
        "AACCGT",
        3, 7,
        soft_left, soft_right
    );
    CHECK(info.cigar == "5=1S");
    CHECK(info.ed == 0);
    CHECK(info.global_ed == 1);
    CHECK(info.length == 5);
    CHECK(info.sw_score == 5 * 3);
    CHECK(soft_left == 0);
    CHECK(soft_right == 1);
    CHECK(info.ref_offset == soft_left);

    // negative total score, soft clipping on both ends
    info = hamming_align(
        "NAAAAAAAAAAAAAA",
        "TAAAATTTTTTTTTT",
        3, 7,
        soft_left, soft_right
    );
    CHECK(info.cigar == "1S4=10S");
    CHECK(info.ed == 0);
    CHECK(info.global_ed == 11);
    CHECK(info.length == 4);
    CHECK(info.sw_score == 4 * 3);
    CHECK(soft_left == 1);
    CHECK(soft_right == 10);
    CHECK(info.ref_offset == soft_left);

    info = hamming_align(
        "NAAAAAAAAAAAAAA",
        "TAAAATTTAAAAAAT",
        3, 7,
        soft_left, soft_right
    );
    CHECK(info.cigar == "8S6=1S");
    CHECK(info.ed == 0);
    CHECK(info.global_ed == 9);
    CHECK(info.length == 6);
    CHECK(info.sw_score == 6 * 3);
    CHECK(soft_left == 8);
    CHECK(soft_right == 1);
    CHECK(info.ref_offset == soft_left);

    info = hamming_align(
        "AAAAAAAAAAAAAAA",
        "TAAAAAATTTAAAAT",
        3, 7,
        soft_left, soft_right
    );
    CHECK(info.cigar == "1S6=8S");
    CHECK(info.ed == 0);
    CHECK(info.global_ed == 9);
    CHECK(info.length == 6);
    CHECK(info.sw_score == 6 * 3);
    CHECK(soft_left == 1);
    CHECK(soft_right == 8);
    CHECK(info.ref_offset == soft_left);
}

TEST_CASE("highest_scoring_segment") {
    auto x = highest_scoring_segment("", "", 5, 7);
    CHECK(x.first == 0);
    CHECK(x.second == 0);
    x = highest_scoring_segment("AAAAAAAAAA", "AAAAAATTTT", 5, 7);
    CHECK(x.first == 0);
    CHECK(x.second == 6);
    x = highest_scoring_segment("AAAAAAAAAA", "TTTTAAAAAA", 5, 7);
    CHECK(x.first == 4);
    CHECK(x.second == 10);
    CHECK(highest_scoring_segment("AAAAAAAAAA", "AAAAAATTTT", 5, 7) == std::make_pair(0ul, 6ul));
    CHECK(highest_scoring_segment("AAAAAAAAAA", "TTTTAAAAAA", 5, 7) == std::make_pair(4ul, 10ul));
    CHECK(highest_scoring_segment("AAAAAAAAAA", "TTAAAAAATT", 5, 7) == std::make_pair(2ul, 8ul));
    CHECK(highest_scoring_segment("AAAAAAAAAAAAAAA", "TAAAAAATTTAAAAT", 5, 7) == std::make_pair(1ul, 7ul));
}
