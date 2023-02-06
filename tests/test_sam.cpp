#include "doctest.h"
#include "sam.hpp"
#include "revcomp.hpp"

TEST_CASE("Formatting unmapped SAM record") {
    klibpp::KSeq kseq;
    kseq.name = "read1";
    kseq.seq = "ACGT";
    kseq.qual = ">#BB";
    std::string sam_string;

    SUBCASE("without RG") {
        Sam sam(sam_string, References());
        sam.add_unmapped(kseq);

        CHECK(sam_string == "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t>#BB\n");
    }

    SUBCASE("with RG") {
        Sam sam(sam_string, References(), "rg1");
        sam.add_unmapped(kseq);

        CHECK(sam_string == "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t>#BB\tRG:Z:rg1\n");
    }
}

TEST_CASE("Pair with one unmapped SAM record") {
    References references;
    references.add("contig1", "ACGT");
    std::string sam_string;
    Sam sam(sam_string, references);

    alignment aln1;
    aln1.ref_id = 0;
    aln1.is_unaligned = false;
    aln1.ref_start = 2;
    aln1.is_rc = true;
    aln1.ed = 17;
    aln1.aln_score = 9;
    aln1.cigar = "2M";

    alignment aln2;
    aln2.is_unaligned = true;

    klibpp::KSeq record1;
    klibpp::KSeq record2;
    record1.name = "readname";
    record1.seq = "AACC";
    record1.qual = "#!B<";
    record2.name = "readname";
    record2.seq = "GGTT";
    record2.qual = "IHB#";
    std::string read1_rc = reverse_complement(record1.seq);
    std::string read2_rc = reverse_complement(record2.seq);

    int mapq1 = 55;
    int mapq2 = 57;
    bool is_proper = false;
    bool is_primary = true;

    sam.add_pair(
        aln1,
        aln2,
        record1,
        record2,
        read1_rc,
        read2_rc,
        mapq1,
        mapq2,
        is_proper,
        is_primary
    );
    // 89: PAIRED,MUNMAP,REVERSE,READ1
    // 165: PAIRED,UNMAP,MREVERSE,READ2
    CHECK(sam_string ==
      "readname\t89\tcontig1\t3\t55\t2M\t*\t0\t0\tGGTT\t<B!#\tNM:i:17\tAS:i:9\n"
      "readname\t165\t*\t0\t0\t*\tcontig1\t3\t0\tGGTT\tIHB#\n"
    );
}

TEST_CASE("TLEN zero when reads map to different contigs") {
    References references;
    references.add("contig1", "ACGT");
    references.add("contig2", "GGAA");
    std::string sam_string;

    alignment aln1;
    aln1.ref_id = 0;
    aln1.is_unaligned = false;
    aln1.ref_start = 2;
    aln1.is_rc = false;
    aln1.ed = 17;
    aln1.aln_score = 9;
    aln1.cigar = "2M";

    alignment aln2;
    aln2.is_unaligned = false;
    aln2.ref_id = 1;
    aln2.ref_start = 3;
    aln2.is_rc = false;
    aln2.ed = 2;
    aln2.aln_score = 4;
    aln2.cigar = "3M";

    klibpp::KSeq record1;
    klibpp::KSeq record2;
    record1.name = "readname";
    record1.seq = "AACC";
    record1.qual = "#!B<";
    record2.name = "readname";
    record2.seq = "GGTT";
    record2.qual = "IHB#";
    std::string read1_rc = reverse_complement(record1.seq);
    std::string read2_rc = reverse_complement(record2.seq);

    int mapq1 = 55;
    int mapq2 = 57;
    bool is_proper = false;
    bool is_primary = true;

    Sam sam(sam_string, references);

    sam.add_pair(
        aln1,
        aln2,
        record1,
        record2,
        read1_rc,
        read2_rc,
        mapq1,
        mapq2,
        is_proper,
        is_primary
    );
    // 65: PAIRED,READ1
    // 129: PAIRED,READ2
    CHECK(sam_string ==
    "readname\t65\tcontig1\t3\t55\t2M\tcontig2\t4\t0\tAACC\t#!B<\tNM:i:17\tAS:i:9\n"
    "readname\t129\tcontig2\t4\t57\t3M\tcontig1\t3\t0\tGGTT\tIHB#\tNM:i:2\tAS:i:4\n"
    );
}
