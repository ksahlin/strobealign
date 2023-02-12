#include "sam.hpp"
#include <algorithm>
#include <ostream>

#define SAM_UNMAPPED_MAPQ 0
#define SAM_UNMAPPED_MAPQ_STRING "0"

using namespace klibpp;

/*
 * SAM columns
 *
 * 1 QNAME query template name
 * 2 FLAG  bitwise flag
 * 3 RNAME reference sequence name
 * 4 POS   1-based leftmost mapping position (set to 0 for unmapped read w/o pos)
 * 5 MAPQ  mapping quality (255 if unavailable)
 * 6 CIGAR
 * 7 RNEXT reference name of mate/next read
 * 8 PNEXT position of mate/next read
 * 9 TLEN template length
 * 10 SEQ
 * 11 QUAL
 * 12 optional fields
*/

/* Strip the /1 or /2 suffix from a read name */
std::string strip_suffix(const std::string& name) {
    auto len = name.length();
    if (
        len >= 2
        && name[len - 2] == '/'
        && (name[len - 1] == '1' || name[len - 1] == '2')
    ) {
        return std::string{name, 0, len - 2};
    } else {
        return name;
    }
}

void Sam::append_tail() {
    sam_string.append(tail);
}

void Sam::add_unmapped(const KSeq& record, int flags) {
    if (!output_unmapped) {
        return;
    }
    assert((flags & ~(UNMAP|PAIRED|MUNMAP|READ1|READ2)) == 0);
    assert(flags & UNMAP);
    sam_string.append(strip_suffix(record.name));
    sam_string.append("\t");
    sam_string.append(std::to_string(flags));
    sam_string.append("\t*\t0\t" SAM_UNMAPPED_MAPQ_STRING "\t*\t*\t0\t0\t");
    sam_string.append(record.seq);
    sam_string.append("\t");
    sam_string.append(record.qual);
    append_tail();
}

void Sam::add_unmapped_mate(const KSeq& record, int flags, const std::string& mate_rname, int mate_pos) {
    assert(flags & (UNMAP|PAIRED));
    sam_string.append(strip_suffix(record.name));
    sam_string.append("\t");
    sam_string.append(std::to_string(flags));
    sam_string.append("\t*\t0\t" SAM_UNMAPPED_MAPQ_STRING "\t*\t");
    sam_string.append(mate_rname);
    sam_string.append("\t");
    sam_string.append(std::to_string(mate_pos + 1));
    sam_string.append("\t0\t");
    sam_string.append(record.seq);
    sam_string.append("\t");
    sam_string.append(record.qual);
    append_tail();
}

void Sam::add_unmapped_pair(const KSeq& r1, const KSeq& r2) {
    add_unmapped(r1, PAIRED | UNMAP | MUNMAP | READ1);
    add_unmapped(r2, PAIRED | UNMAP | MUNMAP | READ2);
}

// Add single-end alignment
void Sam::add(
    const alignment& sam_aln,
    const KSeq& record,
    const std::string& sequence_rc,
    bool is_secondary
) {
    assert(!sam_aln.is_unaligned);
    int flags = 0;
    if (!sam_aln.is_unaligned && sam_aln.is_rc) {
        flags |= REVERSE;
    }
    if (is_secondary) {
        flags |= SECONDARY;
    }
    add_record(record.name, flags, references.names[sam_aln.ref_id], sam_aln.ref_start, sam_aln.mapq, sam_aln.cigar, "*", -1, 0, record.seq, sequence_rc, record.qual, sam_aln.ed, sam_aln.aln_score);
}

// Add one individual record
void Sam::add_record(
    const std::string& query_name,
    int flags,
    const std::string& reference_name,
    int pos,
    int mapq,
    const std::string& cigar,
    const std::string& mate_name,
    int mate_ref_start,
    int template_len,
    const std::string& query_sequence,
    const std::string& query_sequence_rc,
    const std::string& qual,
    int ed,
    int aln_score
) {
    sam_string.append(strip_suffix(query_name));
    sam_string.append("\t");
    sam_string.append(std::to_string(flags));
    sam_string.append("\t");
    sam_string.append(reference_name);
    sam_string.append("\t");
    sam_string.append(std::to_string(pos + 1));  // convert to 1-based coordinate
    sam_string.append("\t");
    sam_string.append(std::to_string(mapq));
    sam_string.append("\t");
    sam_string.append(cigar);
    sam_string.append("\t");

    sam_string.append(mate_name);
    sam_string.append("\t");
    sam_string.append(std::to_string(mate_ref_start + 1));
    sam_string.append("\t");
    sam_string.append(std::to_string(template_len));
    sam_string.append("\t");

    if (flags & REVERSE) {
        sam_string.append(query_sequence_rc);
    } else {
        sam_string.append(query_sequence);
    }
    sam_string.append("\t");

    if (!(flags & UNMAP)) {
        if (flags & REVERSE) {
            auto qual_rev = qual;
            std::reverse(qual_rev.begin(), qual_rev.end());
            sam_string.append(qual_rev);
        } else {
            sam_string.append(qual);
        }
        sam_string.append("\t");
        sam_string.append("NM:i:");
        sam_string.append(std::to_string(ed));
        sam_string.append("\t");
        sam_string.append("AS:i:");
        sam_string.append(std::to_string(aln_score));
    } else {
        sam_string.append(qual);
    }
    append_tail();
}

void Sam::add_pair(
    const alignment &sam_aln1,
    const alignment &sam_aln2,
    const KSeq& record1,
    const KSeq& record2,
    const std::string &read1_rc,
    const std::string &read2_rc,
    int mapq1,
    int mapq2,
    bool is_proper,
    bool is_primary
) {
    int f1 = PAIRED | READ1;
    int f2 = PAIRED | READ2;
    if (!is_primary) {
        f1 |= SECONDARY;
        f2 |= SECONDARY;
    }

    int template_len1 = 0;
    bool both_aligned = !sam_aln1.is_unaligned && !sam_aln2.is_unaligned;
    if (both_aligned && sam_aln1.ref_id == sam_aln2.ref_id) {
        const int dist = sam_aln2.ref_start - sam_aln1.ref_start;
        if (dist > 0) {
            template_len1 = dist + sam_aln2.aln_length;
        }
        else {
            template_len1 = dist - sam_aln1.aln_length;
        }
    }
    if (is_proper) {
        f1 |= PROPER_PAIR;
        f2 |= PROPER_PAIR;
    }

    std::string mate_name1;
    std::string ref1;
    int ref_start1 = sam_aln1.ref_start;
    int ed1 = sam_aln1.ed;
    if (sam_aln1.is_unaligned) {
        f1 |= UNMAP;
        f2 |= MUNMAP;
        ref1 = "*";
        ref_start1 = -1;
        mate_name1 = "*";
    } else {
        if (sam_aln1.is_rc) {
            f1 |= REVERSE;
            f2 |= MREVERSE;
        }
        mate_name1 = references.names[sam_aln1.ref_id];
        ref1 = references.names[sam_aln1.ref_id];
    }

    std::string mate_name2;
    std::string ref2;
    int ref_start2 = sam_aln2.ref_start;
    int ed2 = sam_aln2.ed;
    if (sam_aln2.is_unaligned) {
        f2 |= UNMAP;
        f1 |= MUNMAP;
        ref_start2 = -1;
        ref2 = "*";
        mate_name2 = "*";
    } else {
        if (sam_aln2.is_rc) {
            f1 |= MREVERSE;
            f2 |= REVERSE;
        }
        mate_name2 = references.names[sam_aln2.ref_id];
        ref2 = references.names[sam_aln2.ref_id];
    }
    if (!sam_aln1.is_unaligned && !sam_aln2.is_unaligned && sam_aln1.ref_id == sam_aln2.ref_id) {
        mate_name1 = "=";
        mate_name2 = "=";
    }

    if (sam_aln1.is_unaligned) {
        add_unmapped_mate(record1, f1, mate_name2, ref_start2);
    } else {
        add_record(record1.name, f1, ref1, sam_aln1.ref_start, mapq1, sam_aln1.cigar, mate_name2, ref_start2, template_len1, record1.seq, read1_rc, record1.qual, ed1, sam_aln1.aln_score);
    }
    if (sam_aln2.is_unaligned) {
        add_unmapped_mate(record2, f2, mate_name1, ref_start1);
    } else {
        add_record(record2.name, f2, ref2, sam_aln2.ref_start, mapq2, sam_aln2.cigar, mate_name1, ref_start1, -template_len1, record2.seq, read2_rc, record2.qual, ed2, sam_aln2.aln_score);
    }
}

bool is_proper_pair(const alignment& sam_aln1, const alignment& sam_aln2, float mu, float sigma) {
    const int dist = sam_aln2.ref_start - sam_aln1.ref_start;
    const bool same_reference = sam_aln1.ref_id == sam_aln2.ref_id;
    const bool both_aligned = same_reference && !sam_aln1.is_unaligned && !sam_aln2.is_unaligned;
    const bool r1_r2 = !sam_aln1.is_rc && sam_aln2.is_rc && dist >= 0; // r1 ---> <---- r2
    const bool r2_r1 = !sam_aln2.is_rc && sam_aln1.is_rc && dist <= 0; // r2 ---> <---- r1
    const bool rel_orientation_good = r1_r2 || r2_r1;
    const bool insert_good = std::abs(dist) <= mu + 6 * sigma;

    return both_aligned && insert_good && rel_orientation_good;
}
