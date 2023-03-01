#include "doctest.h"
#include "cigar.hpp"

TEST_CASE("compress_cigar") {
    CHECK(compress_cigar("") == "");
    CHECK(compress_cigar("M") == "1M");
    CHECK(compress_cigar("MM") == "2M");
    CHECK(compress_cigar("MMI") == "2M1I");
    CHECK(compress_cigar("MMII") == "2M2I");
    CHECK(compress_cigar("MI") == "1M1I");
    CHECK(compress_cigar("MII") == "1M2I");
}

TEST_CASE("Parse CIGAR") {
    std::vector<std::string> cigars = {"", "1M", "10M2I1D99=1X4P5S"};
    for (auto& s : cigars) {
        CHECK(Cigar(s).to_string() == s);
    }
    // Not standard, only for convenience
    CHECK(Cigar("M").to_string() == "1M");
    CHECK(Cigar("M M").to_string() == "2M");
    CHECK(Cigar("MMII").to_string() == "2M2I");
    CHECK(Cigar("M 2M X").to_string() == "3M1X");
    CHECK(Cigar("1M 2D 1M").to_string() == "1M2D1M");
}

TEST_CASE("Cigar construction and push") {
    Cigar c1;
    CHECK(c1.to_string() == "");

    c1.push(CIGAR_MATCH, 1);
    CHECK(c1.to_string() == "1M");

    c1.push(CIGAR_MATCH, 1);
    CHECK(c1.to_string() == "2M");

    c1.push(CIGAR_INS, 3);
    CHECK(c1.to_string() == "2M3I");

    uint32_t ops[4] = {
        3 << 4 | CIGAR_MATCH,
        5 << 4 | CIGAR_X,
        7 << 4 | CIGAR_INS,
        13 << 4 | CIGAR_DEL
    };
    Cigar c2{ops, 4};
    CHECK(c2.m_ops.size() == 4);
    CHECK(c2.to_string() == "3M5X7I13D");
}

TEST_CASE("Cigar =/X conversion") {
    Cigar c1{"1M"};
    CHECK(c1.to_eqx("A", "A").to_string() == "1=");
    CHECK(c1.to_eqx("A", "G").to_string() == "1X");

    Cigar c2{"3M"};
    CHECK(c2.to_eqx("AAA", "AAA").to_string() == "3=");
    CHECK(c2.to_eqx("AAA", "ATA").to_string() == "1=1X1=");
    CHECK(c2.to_eqx("AAA", "ATT").to_string() == "1=2X");

    Cigar c{"2M 1D 4M 1I 3M"};
    CHECK(c.to_eqx("ACTTTGCATT", "ACGTATGAAA").to_string() == "2=1D1=1X2=1I1=2X");
}

TEST_CASE("concatenate Cigar") {
    Cigar c{"3M"};
    c += Cigar{"2M1X"};
    CHECK(c.to_string() == "5M1X");
}

TEST_CASE("edit distance") {
    CHECK(Cigar("3=1X4D5I7=").edit_distance() == 10);
}

TEST_CASE("reverse") {
    Cigar c{"3=1X4D5I7="};
    c.reverse();
    CHECK(c.to_string() == "7=5I4D1X3=");
}
