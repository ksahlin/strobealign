#include "doctest.h"
#include "sam.hpp"
#include "revcomp.hpp"

TEST_CASE("Formatting unmapped SAM record") {
    klibpp::KSeq kseq;
    kseq.name = "read1";
    kseq.seq = "ACGT";
    kseq.qual = ">#BB";
    std::string sam_string;
    References references;

    SUBCASE("without RG") {
        Sam sam(sam_string, references);
        sam.add_unmapped(kseq);

        CHECK(sam_string == "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t>#BB\n");
    }

    SUBCASE("with RG") {
        Sam sam(sam_string, references, CigarOps::EQX, "rg1");
        sam.add_unmapped(kseq);

        CHECK(sam_string == "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t>#BB\tRG:Z:rg1\n");
    }

    SUBCASE("without qualities") {
        Sam sam(sam_string, references);
        klibpp::KSeq kseq_noqual;
        kseq_noqual.name = "read1";
        kseq_noqual.seq = "ACGT";
        kseq_noqual.qual = "";
        sam.add_unmapped(kseq_noqual);

        CHECK(sam_string == "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n");
    }
}

TEST_CASE("Sam::add") {
    References references;
    references.add("contig1", "AACCGGTT");

    klibpp::KSeq record;
    record.name = "readname";
    record.seq = "AACGT";
    record.qual = ">#BB";

    Alignment aln;
    aln.ref_id = 0;
    aln.is_unaligned = false;
    aln.is_rc = true;
    aln.ref_start = 2;
    aln.ed = 3;
    aln.sw_score = 9;
    aln.mapq = 55;
    aln.cigar = Cigar("2S2=1X3=3S");

    std::string read_rc = reverse_complement(record.seq);
    bool is_primary = true;
    Details details;
    SUBCASE("Cigar =/X") {
        std::string sam_string;
        Sam sam(sam_string, references);
        sam.add(aln, record, read_rc, is_primary, details);
        CHECK(sam_string ==
            "readname\t16\tcontig1\t3\t55\t2S2=1X3=3S\t*\t0\t0\tACGTT\tBB#>\tNM:i:3\tAS:i:9\n"
        );
    }
    SUBCASE("Cigar M") {
        std::string sam_string;
        Sam sam(sam_string, references, CigarOps::M);
        sam.add(aln, record, read_rc, is_primary, details);
        CHECK(sam_string ==
            "readname\t16\tcontig1\t3\t55\t2S6M3S\t*\t0\t0\tACGTT\tBB#>\tNM:i:3\tAS:i:9\n"
        );
    }
}

TEST_CASE("Pair with one unmapped SAM record") {
    References references;
    references.add("contig1", "ACGT");
    std::string sam_string;
    Sam sam(sam_string, references);

    Alignment aln1;
    aln1.ref_id = 0;
    aln1.is_unaligned = false;
    aln1.ref_start = 2;
    aln1.is_rc = true;
    aln1.ed = 17;
    aln1.sw_score = 9;
    aln1.cigar = Cigar("2M");

    Alignment aln2;
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
    std::array<Details, 2> details;

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
        is_primary,
        details
    );
    // 89: PAIRED,MUNMAP,REVERSE,READ1
    // 165: PAIRED,UNMAP,MREVERSE,READ2
    CHECK(sam_string ==
      "readname\t89\tcontig1\t3\t55\t2M\t=\t3\t0\tGGTT\t<B!#\tNM:i:17\tAS:i:9\n"
      "readname\t165\tcontig1\t3\t0\t*\t=\t3\t0\tGGTT\tIHB#\n"
    );
}

TEST_CASE("TLEN zero when reads map to different contigs") {
    References references;
    references.add("contig1", "ACGT");
    references.add("contig2", "GGAA");
    std::string sam_string;

    Alignment aln1;
    aln1.ref_id = 0;
    aln1.is_unaligned = false;
    aln1.ref_start = 2;
    aln1.is_rc = false;
    aln1.ed = 17;
    aln1.sw_score = 9;
    aln1.cigar = Cigar("2M");

    Alignment aln2;
    aln2.is_unaligned = false;
    aln2.ref_id = 1;
    aln2.ref_start = 3;
    aln2.is_rc = false;
    aln2.ed = 2;
    aln2.sw_score = 4;
    aln2.cigar = Cigar("3M");

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
    std::array<Details, 2> details;

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
        is_primary,
        details
    );
    // 65: PAIRED,READ1
    // 129: PAIRED,READ2
    CHECK(sam_string ==
    "readname\t65\tcontig1\t3\t55\t2M\tcontig2\t4\t0\tAACC\t#!B<\tNM:i:17\tAS:i:9\n"
    "readname\t129\tcontig2\t4\t57\t3M\tcontig1\t3\t0\tGGTT\tIHB#\tNM:i:2\tAS:i:4\n"
    );
}
