#include "doctest.h"
#include "aligner.hpp"

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
    CHECK(info.ref_start == soft_left);

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
    CHECK(info.ref_start == soft_left);

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
    CHECK(info.ref_start == soft_left);

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
    CHECK(info.ref_start == soft_left);

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
    CHECK(info.ref_start == soft_left);

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
    CHECK(info.ref_start == soft_left);

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
    CHECK(info.ref_start == soft_left);

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
    CHECK(info.ref_start == soft_left);
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
