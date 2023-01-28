#include <vector>
#include "doctest.h"
#include "pc.hpp"
#include "kseq++.hpp"

TEST_CASE("InputBuffer interleaved") {
    InputBuffer ibuf("tests/interleaved.fq", "", 1, true);
    std::vector<klibpp::KSeq> records1;
    std::vector<klibpp::KSeq> records2;
    std::vector<klibpp::KSeq> records3;
    AlignmentStatistics stats;
    int total_pe = 0;
    int total_se = 0;
    while (true) {
        ibuf.read_records(records1,
                            records2,
                            records3,
                            stats);
        CHECK(records1.size() == records2.size());
        total_pe += records1.size();
        total_se += records3.size();
        if (records1.empty() && records3.empty()) {
            break;
        }
    }
    CHECK(total_pe == 2);
    CHECK(total_se == 4);
}

TEST_CASE("InputBuffer paired") {
    InputBuffer ibuf("tests/phix.1.fastq", "tests/phix.2.fastq", 3, false);
    std::vector<klibpp::KSeq> records1;
    std::vector<klibpp::KSeq> records2;
    std::vector<klibpp::KSeq> records3;
    AlignmentStatistics stats;
    int total_pe = 0;
    int total_se = 0;
    while (true) {
        ibuf.read_records(records1,
                            records2,
                            records3,
                            stats);
        CHECK(records1.size() == records2.size());
        total_pe += records1.size();
        total_se += records3.size();
        if (records1.empty() && records3.empty()) {
            break;
        }
    }
    CHECK(total_pe == 44);
    CHECK(total_se == 0);
}

TEST_CASE("InputBuffer single-end (with rewind)") {
    InputBuffer ibuf("tests/phix.1.fastq", "", 3, false);
    std::vector<klibpp::KSeq> records1;
    std::vector<klibpp::KSeq> records2;
    std::vector<klibpp::KSeq> records3;
    AlignmentStatistics stats;
    int total_pe = 0;
    int total_se = 0;
    while (true) {
        ibuf.read_records(records1,
                            records2,
                            records3,
                            stats);
        CHECK(records1.size() == records2.size());
        total_pe += records1.size();
        total_se += records3.size();
        if (records1.empty() && records3.empty()) {
            break;
        }
    }
    CHECK(total_pe == 0);
    CHECK(total_se == 44);

    ibuf.rewind_reset();
    total_pe = 0;
    total_se = 0;
    while (true) {
        ibuf.read_records(records1,
                            records2,
                            records3,
                            stats);
        CHECK(records1.size() == records2.size());
        total_pe += records1.size();
        total_se += records3.size();
        if (records1.empty() && records3.empty()) {
            break;
        }
    }
    CHECK(total_pe == 0);
    CHECK(total_se == 44);
}

TEST_CASE("RewindableFile") {
    RewindableFile rf("tests/phix.1.fastq");
    char buf1[1024];
    char buf2[1024];
    rf.read(buf1, 512);
    rf.read(buf1 + 512, 256);
    rf.read(buf1 + 512 + 256, 256);
    rf.rewind();
    rf.read(buf2, 1024);
    CHECK(memcmp(buf1, buf2, 1024) == 0);
}
