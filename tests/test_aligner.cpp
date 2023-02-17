#include "doctest.h"
#include "aligner.hpp"

TEST_CASE("hamming_align") {
    // empty sequences
    auto info = hamming_align(
        "", "",
        7, 5
    );
    CHECK(info.cigar == "");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 0);
    CHECK(info.ref_start == 0);
    CHECK(info.ref_span() == 0);
    CHECK(info.query_start == 0);
    CHECK(info.query_end == 0);

    info = hamming_align(
        "AAXGGG",
        "AAYGGG",
        1, 1
    );
    CHECK(info.cigar == "2=1X3=");
    CHECK(info.edit_distance == 1);
    CHECK(info.sw_score == 4);
    CHECK(info.ref_start == 0);
    CHECK(info.ref_span() == 6);
    CHECK(info.query_start == 0);
    CHECK(info.query_end == 6);

    info = hamming_align(
        "AXGGG",
        "AYGGG",
        1, 4
    );
    CHECK(info.cigar == "2S3=");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 3);
    CHECK(info.ref_start == 2);
    CHECK(info.ref_span() == 3);
    CHECK(info.query_start == 2);
    CHECK(info.query_end == 5);

    info = hamming_align(
        "NAACCG",
        "TAACCG",
        3, 7
    );
    CHECK(info.cigar == "1S5=");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 5 * 3);
    CHECK(info.ref_start == 1);
    CHECK(info.ref_span() == 5);
    CHECK(info.query_start == 1);
    CHECK(info.query_end == 6);

    info = hamming_align(
        "AACCGN",
        "AACCGT",
        3, 7
    );
    CHECK(info.cigar == "5=1S");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 5 * 3);
    CHECK(info.ref_start == 0);
    CHECK(info.ref_span() == 5);
    CHECK(info.query_start == 0);
    CHECK(info.query_end == 5);

    // negative total score, soft clipping on both ends
    info = hamming_align(
        "NAAAAAAAAAAAAAA",
        "TAAAATTTTTTTTTT",
        3, 7
    );
    CHECK(info.cigar == "1S4=10S");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 4 * 3);
    CHECK(info.ref_start == 1);
    CHECK(info.ref_span() == 4);
    CHECK(info.query_start == 1);
    CHECK(info.query_end == 5);

    info = hamming_align(
        "NAAAAAAAAAAAAAA",
        "TAAAATTTAAAAAAT",
        3, 7
    );
    CHECK(info.cigar == "8S6=1S");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 6 * 3);
    CHECK(info.ref_start == 8);
    CHECK(info.ref_span() == 6);
    CHECK(info.query_start == 8);
    CHECK(info.query_end == 14);

    info = hamming_align(
        "AAAAAAAAAAAAAAA",
        "TAAAAAATTTAAAAT",
        3, 7
    );
    CHECK(info.cigar == "1S6=8S");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 6 * 3);
    CHECK(info.ref_start == 1);
    CHECK(info.ref_span() == 6);
    CHECK(info.query_start == 1);
    CHECK(info.query_end == 7);
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
