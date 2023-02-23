#include "doctest.h"
#include "aligner.hpp"

TEST_CASE("hamming_align") {
    // empty sequences
    auto info = hamming_align(
        "", "",
        7, 5, 0
    );
    CHECK(info.cigar.to_string() == "");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 0);
    CHECK(info.ref_start == 0);
    CHECK(info.ref_span() == 0);
    CHECK(info.query_start == 0);
    CHECK(info.query_end == 0);

    info = hamming_align(
        "AAXGGG",
        "AAYGGG",
        1, 1, 0
    );
    CHECK(info.cigar.to_string() == "2=1X3=");
    CHECK(info.edit_distance == 1);
    CHECK(info.sw_score == 4);
    CHECK(info.ref_start == 0);
    CHECK(info.ref_span() == 6);
    CHECK(info.query_start == 0);
    CHECK(info.query_end == 6);

    info = hamming_align(
        "AXGGG",
        "AYGGG",
        1, 4, 0
    );
    CHECK(info.cigar.to_string() == "2S3=");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 3);
    CHECK(info.ref_start == 2);
    CHECK(info.ref_span() == 3);
    CHECK(info.query_start == 2);
    CHECK(info.query_end == 5);

    info = hamming_align(
        "NAACCG",
        "TAACCG",
        3, 7, 0
    );
    CHECK(info.cigar.to_string() == "1S5=");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 5 * 3);
    CHECK(info.ref_start == 1);
    CHECK(info.ref_span() == 5);
    CHECK(info.query_start == 1);
    CHECK(info.query_end == 6);

    info = hamming_align(
        "AACCGN",
        "AACCGT",
        3, 7, 0
    );
    CHECK(info.cigar.to_string() == "5=1S");
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
        3, 7, 0
    );
    CHECK(info.cigar.to_string() == "1S4=10S");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 4 * 3);
    CHECK(info.ref_start == 1);
    CHECK(info.ref_span() == 4);
    CHECK(info.query_start == 1);
    CHECK(info.query_end == 5);

    info = hamming_align(
        "NAAAAAAAAAAAAAA",
        "TAAAATTTAAAAAAT",
        3, 7, 0
    );
    CHECK(info.cigar.to_string() == "8S6=1S");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 6 * 3);
    CHECK(info.ref_start == 8);
    CHECK(info.ref_span() == 6);
    CHECK(info.query_start == 8);
    CHECK(info.query_end == 14);

    info = hamming_align(
        "AAAAAAAAAAAAAAA",
        "TAAAAAATTTAAAAT",
        3, 7, 0
    );
    CHECK(info.cigar.to_string() == "1S6=8S");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 6 * 3);
    CHECK(info.ref_start == 1);
    CHECK(info.ref_span() == 6);
    CHECK(info.query_start == 1);
    CHECK(info.query_end == 7);
}

TEST_CASE("highest_scoring_segment") {
    auto x = highest_scoring_segment("", "", 5, 7, 0);
    CHECK(std::get<0>(x) == 0);
    CHECK(std::get<1>(x) == 0);
    x = highest_scoring_segment("AAAAAAAAAA", "AAAAAATTTT", 5, 7, 0);
    CHECK(std::get<0>(x) == 0);
    CHECK(std::get<1>(x) == 6);
    x = highest_scoring_segment("AAAAAAAAAA", "TTTTAAAAAA", 5, 7, 0);
    CHECK(std::get<0>(x) == 4);
    CHECK(std::get<1>(x) == 10);
    CHECK(highest_scoring_segment("AAAAAAAAAA", "AAAAAATTTT", 5, 7, 0) == std::make_tuple(0ul, 6ul, 30));
    CHECK(highest_scoring_segment("AAAAAAAAAA", "TTTTAAAAAA", 5, 7, 0) == std::make_tuple(4ul, 10ul, 30));
    CHECK(highest_scoring_segment("AAAAAAAAAA", "TTAAAAAATT", 5, 7, 0) == std::make_tuple(2ul, 8ul, 30));
    CHECK(highest_scoring_segment("AAAAAAAAAAAAAAA", "TAAAAAATTTAAAAT", 5, 7, 0) == std::make_tuple(1ul, 7ul, 30));
}

TEST_CASE("highest_scoring_segment with soft clipping") {
    auto x = highest_scoring_segment("", "", 2, 4, 5);
    CHECK(std::get<0>(x) == 0);
    CHECK(std::get<1>(x) == 0);
    CHECK(std::get<2>(x) == 10);

    x = highest_scoring_segment("TAAT", "TAAA", 2, 4, 5);
    CHECK(std::get<0>(x) == 0);
    CHECK(std::get<1>(x) == 4);
    CHECK(std::get<2>(x) == 3 * 2 - 4 + 10);

       x = highest_scoring_segment("AAA", "AAA", 2, 4, 5);
    CHECK(std::get<0>(x) == 0);
    CHECK(std::get<1>(x) == 3);
    CHECK(std::get<2>(x) == 3 * 2 + 10);

    x = highest_scoring_segment("TAAT", "AAAA", 2, 4, 5);
    CHECK(std::get<0>(x) == 0);
    CHECK(std::get<1>(x) == 4);
    CHECK(std::get<2>(x) == 2 * 2 - 2 * 4 + 10);

    x = highest_scoring_segment("ATAATA", "AAAAAA", 2, 4, 5);
    CHECK(std::get<0>(x) == 0);
    CHECK(std::get<1>(x) == 6);
    CHECK(std::get<2>(x) == 4 * 2 - 2 * 4 + 10);

    x = highest_scoring_segment("TTAATA", "AAAAAA", 2, 4, 5);
    CHECK(std::get<0>(x) == 2);
    CHECK(std::get<1>(x) == 6);
    CHECK(std::get<2>(x) == 3 * 2 - 1 * 4 + 5);

}
