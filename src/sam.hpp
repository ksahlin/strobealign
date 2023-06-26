#ifndef STROBEALIGN_SAM_HPP
#define STROBEALIGN_SAM_HPP

#include <string>
#include <array>
#include "kseq++/kseq++.hpp"
#include "refs.hpp"
#include "cigar.hpp"


struct Alignment {
    Cigar cigar;
    int ref_start;
    int ed;
    int global_ed;
    int sw_score;
    int aln_score;
    int ref_id;
    int mapq;
    int aln_length;
    bool is_rc;
    bool is_unaligned = false;
};

enum SamFlags {
    PAIRED = 1,
    PROPER_PAIR = 2,
    UNMAP = 4,
    MUNMAP = 8,
    REVERSE = 0x10,
    MREVERSE = 0x20,
    READ1 = 0x40,
    READ2 = 0x80,
    SECONDARY = 0x100,
    QCFAIL = 0x200,
    DUP = 0x400,
    SUPPLEMENTARY = 0x800,
};

enum struct CigarOps {
    EQX = 0,  // use = and X CIGAR operations
    M = 1,    // use M CIGAR operations
};

/* Details about aligning a single or paired-end read */
struct Details {
    bool nam_rescue{false}; // find_nams_rescue() was needed
    uint64_t nams{0};  // No. of NAMs found
    uint64_t mate_rescue{false}; // No. of times rescue by local alignment was attempted
    uint64_t tried_alignment{0}; // No. of computed alignments (get_alignment or rescue_mate)
};

class Sam {

public:
    Sam(
        std::string& sam_string,
        const References& references,
        CigarOps cigar_ops = CigarOps::EQX,
        const std::string& read_group_id = "",
        bool output_unmapped = true,
        bool show_details = false
    )
        : sam_string(sam_string)
        , references(references)
        , cigar_ops(cigar_ops)
        , output_unmapped(output_unmapped)
        , show_details(show_details)
    {
            if (read_group_id.empty()) {
                tail = "\n";
            } else {
                tail = "\tRG:Z:" + read_group_id + "\n";
            }
        }

    /* Add an alignment */
    void add(const Alignment& sam_aln, const klibpp::KSeq& record, const std::string& sequence_rc, bool is_primary, const Details& details);
    void add_pair(const Alignment& sam_aln1, const Alignment& sam_aln2, const klibpp::KSeq& record1, const klibpp::KSeq& record2, const std::string& read1_rc, const std::string& read2_rc, int mapq1, int mapq2, bool is_proper, bool is_primary, const std::array<Details, 2>& details);
    void add_unmapped(const klibpp::KSeq& record, int flags = UNMAP);
    void add_unmapped_pair(const klibpp::KSeq& r1, const klibpp::KSeq& r2);
    void add_unmapped_mate(const klibpp::KSeq& record, int flags, const std::string& mate_reference_name, int mate_pos);

private:
    void add_record(const std::string& query_name, int flags, const std::string& reference_name, int pos, int mapq, const Cigar& cigar, const std::string& mate_reference_name, int mate_pos, int template_len, const std::string& query_sequence, const std::string& query_sequence_rc, const std::string& qual, int ed, int aln_score, const Details& details);
    void append_qual(const std::string& qual) {
        sam_string.append("\t");
        sam_string.append(qual.empty() ? "*" : qual);
    }
    void append_details(const Details& details);
    void append_paired_details(const Details& details);
    void append_rg_and_newline();

    std::string cigar_string(const Cigar& cigar) const;
    std::string& sam_string;
    const References& references;
    const CigarOps cigar_ops;
    std::string tail;
    bool output_unmapped;
    bool show_details;
};

bool is_proper_pair(const Alignment& sam_aln1, const Alignment& sam_aln2, float mu, float sigma);

#endif
