#ifndef STROBEALIGN_SAM_HPP
#define STROBEALIGN_SAM_HPP

#include <string>
#include <array>
#include "kseq++/kseq++.hpp"
#include "refs.hpp"
#include "cigar.hpp"


struct Alignment {
    int ref_id;
    int ref_start;
    Cigar cigar;
    int edit_distance;
    int global_ed;
    int score;
    int length;
    bool is_rc;
    bool is_unaligned{false};
    // Whether a gapped alignment function was used to obtain this alignment
    // (even if true, the alignment can still be without gaps)
    bool gapped{false};
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
    uint64_t nam_inconsistent{0};
    uint64_t mate_rescue{0}; // No. of times rescue by local alignment was attempted
    uint64_t tried_alignment{0}; // No. of computed alignments (get_alignment or rescue_mate)
    uint64_t gapped{0};  // No. of gapped alignments computed (in get_alignment)
    uint64_t best_alignments{0}; // No. of best alignments with same score
    std::chrono::duration<double> align_duration{0};  // Time spent on extension alignment

    Details& operator+=(const Details& other) {
        nam_rescue = nam_rescue || other.nam_rescue;
        nams += other.nams;
        nam_inconsistent += other.nam_inconsistent;
        mate_rescue += other.mate_rescue;
        tried_alignment += other.tried_alignment;
        gapped += other.gapped;
        best_alignments += other.best_alignments;
        align_duration += other.align_duration;
        return *this;
    }
};

class Sam {

public:
    Sam(
        std::string& sam_string,
        const References& references,
        CigarOps cigar_ops = CigarOps::EQX,
        const std::string& read_group_id = "",
        bool output_unmapped = true,
        bool show_details = false,
        bool fastq_comments = false
    )
        : sam_string(sam_string)
        , references(references)
        , cigar_ops(cigar_ops)
        , output_unmapped(output_unmapped)
        , show_details(show_details)
        , fastq_comments(fastq_comments)
    {
            if (read_group_id.empty()) {
                tail = "";
            } else {
                tail = "\tRG:Z:" + read_group_id;
            }
        }

    /* Add an alignment */
    void add(const Alignment& alignment, const klibpp::KSeq& record, const std::string& sequence_rc, uint8_t mapq, bool is_primary, const Details& details);
    void add_pair(const Alignment& alignment1, const Alignment& alignment2, const klibpp::KSeq& record1, const klibpp::KSeq& record2, const std::string& read1_rc, const std::string& read2_rc, uint8_t mapq1, uint8_t mapq2, bool is_proper, bool is_primary, const std::array<Details, 2>& details);
    void add_unmapped(const klibpp::KSeq& record, uint16_t flags = UNMAP);
    void add_unmapped_pair(const klibpp::KSeq& r1, const klibpp::KSeq& r2);
    void add_unmapped_mate(const klibpp::KSeq& record, uint16_t flags, const std::string& mate_reference_name, uint32_t mate_pos);

private:
    void add_record(const std::string& query_name, const std::string& comment, uint16_t flags, const std::string& reference_name, uint32_t pos, uint8_t mapq, const Cigar& cigar, const std::string& mate_reference_name, uint32_t mate_pos, int32_t template_len, const std::string& query_sequence, const std::string& query_sequence_rc, const std::string& qual, int ed, int aln_score, const Details& details);

    void append_seq(const std::string& seq) {
        sam_string.append(seq.empty() ? "*" : seq);
    }

    void append_qual(const std::string& qual) {
        sam_string.append("\t");
        sam_string.append(qual.empty() ? "*" : qual);
    }
    void append_details(const Details& details);
    void append_paired_details(const Details& details);
    void append_rg();

    std::string cigar_string(const Cigar& cigar) const;
    std::string& sam_string;
    const References& references;
    const CigarOps cigar_ops;
    std::string tail;
    bool output_unmapped;
    bool show_details;
    bool fastq_comments;
};

bool is_proper_pair(const Alignment& alignment1, const Alignment& alignment2, float mu, float sigma);

#endif
