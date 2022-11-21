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
void Sam::add_unmapped(const KSeq& record, int flags) {
    sam_string.append(record.name);
    sam_string.append("\t");
    sam_string.append(std::to_string(flags));
    sam_string.append("\t*\t0\t" SAM_UNMAPPED_MAPQ_STRING "\t*\t*\t0\t0\t");
    sam_string.append(record.seq);
    sam_string.append("\t");
    sam_string.append(record.qual);
    sam_string.append("\n");
}

void Sam::add_unmapped_pair(const KSeq& r1, const KSeq& r2) {
    add_unmapped(r1, PAIRED | UNMAP | MUNMAP | READ1);
    add_unmapped(r2, PAIRED | UNMAP | MUNMAP | READ2);
}

void Sam::add(
    const alignment& sam_aln,
    const KSeq& record,
    const std::string& sequence_rc,
    bool is_secondary
) {
    int flags = 0;
    const std::string* output_read;
    if (!sam_aln.is_unaligned && sam_aln.is_rc) {
        flags |= REVERSE;
        output_read = &sequence_rc;
    } else {
        output_read = &record.seq;
    }
    if (is_secondary) {
        flags |= SECONDARY;
    }
    sam_string.append(record.name);  // QNAME
    sam_string.append("\t");
    sam_string.append(std::to_string(flags));  // FLAG
    sam_string.append("\t");
    sam_string.append(references.names[sam_aln.ref_id]);  // RNAME
    sam_string.append("\t");
    sam_string.append(std::to_string(sam_aln.ref_start));  // POS
    sam_string.append("\t");
    sam_string.append(std::to_string(sam_aln.mapq));  // MAPQ
    sam_string.append("\t");
    sam_string.append(sam_aln.cigar);  // CIGAR
    sam_string.append("\t*\t0\t0\t");  // RNEXT, PNEXT, TLEN
    sam_string.append(*output_read);  // SEQ
    sam_string.append("\t");
    if (!sam_aln.is_unaligned) {
        if (sam_aln.is_rc){
            auto qual_rev = record.qual;
            std::reverse(qual_rev.begin(), qual_rev.end());
            sam_string.append(qual_rev);  // QUAL
        } else {
            sam_string.append(record.qual);
        }
        sam_string.append("\t");
        sam_string.append("NM:i:");
        sam_string.append(std::to_string(sam_aln.ed));
        sam_string.append("\t");
        sam_string.append("AS:i:");
        sam_string.append(std::to_string((int) sam_aln.aln_score));
    } else {
        sam_string.append(record.qual);  // QUAL
    }
    sam_string.append("\n");
}


// Used by add_pair() to add one individual record
void Sam::add_one(
    const KSeq& record,
    int flags,
    const std::string& ref_name,
    const alignment& sam_aln,
    int mapq,
    const std::string& mate_name,
    int mate_ref_start,
    int template_len,
    const std::string& output_read,
    int ed
) {
    sam_string.append(record.name);
    sam_string.append("\t");
    sam_string.append(std::to_string(flags));
    sam_string.append("\t");
    sam_string.append(ref_name);
    sam_string.append("\t");
    sam_string.append(std::to_string(sam_aln.ref_start));
    sam_string.append("\t");
    sam_string.append(std::to_string(mapq));
    sam_string.append("\t");
    sam_string.append(sam_aln.cigar);
    sam_string.append("\t");
    sam_string.append(mate_name);
    sam_string.append("\t");
    sam_string.append(std::to_string(mate_ref_start));
    sam_string.append("\t");
    sam_string.append(std::to_string(template_len));
    sam_string.append("\t");
    if (!sam_aln.is_unaligned) {
        sam_string.append(output_read);
        sam_string.append("\t");
        if (sam_aln.is_rc){
            auto qual_rev = record.qual;
            std::reverse(qual_rev.begin(), qual_rev.end());
            sam_string.append(qual_rev);
        } else {
            sam_string.append(record.qual);
        }
        sam_string.append("\t");
        sam_string.append("NM:i:");
        sam_string.append(std::to_string(ed));
        sam_string.append("\t");
        sam_string.append("AS:i:");
        sam_string.append(std::to_string((int) sam_aln.aln_score));
    } else {
        sam_string.append(record.seq);
        sam_string.append("\t");
        sam_string.append(record.qual);
    }
    sam_string.append("\n");
}

void Sam::add_pair(
    alignment &sam_aln1,
    alignment &sam_aln2,
    const KSeq& record1,
    const KSeq& record2,
    const std::string &read1_rc,
    const std::string &read2_rc,
    int mapq1,
    int mapq2,
    float mu,
    float sigma,
    bool is_primary
) {

    // Commented lines below because we do not longer mark a read as not proper just because of the non-matching hash
    // Proper or non proper reads are further below only decided based on the expected distance and relative orientation they align to
//    if (sam_aln1.ed < 5){ // Flag alignments previously deemed as 'not proper' (based on matching strobemer hash ) to proper because of small ed
//        sam_aln1.is_proper = true;
//    }
//    if (sam_aln2.ed < 5){ // Flag alignments previously deemed as 'not proper' (based on matching strobemer hash ) to proper because of small ed
//        sam_aln2.is_proper = true;
//    }

    const int dist = sam_aln2.ref_start - sam_aln1.ref_start;
    int template_len1;
    if (dist > 0) {
        template_len1 = dist + sam_aln2.aln_length;
    }
    else {
        template_len1 = dist - sam_aln1.aln_length;
    }

    bool both_aligned = !sam_aln1.is_unaligned && !sam_aln2.is_unaligned;
    bool r1_r2 = !sam_aln1.is_rc && sam_aln2.is_rc && dist >= 0; // r1 ---> <---- r2
    bool r2_r1 = !sam_aln2.is_rc && sam_aln1.is_rc && dist <= 0; // r2 ---> <---- r1
    bool rel_orientation_good = r1_r2 || r2_r1;
    bool insert_good = std::abs(dist) <= mu + 6 * sigma;
    if (both_aligned && insert_good && rel_orientation_good) {
        sam_aln1.is_proper = true;
        sam_aln2.is_proper = true;
    } else {
        sam_aln1.is_proper = false;
        sam_aln2.is_proper = false;
    }

    int f1 = PAIRED | READ1;
    int f2 = PAIRED | READ2;
    if (!is_primary) {
        f1 |= SECONDARY;
        f2 |= SECONDARY;
    }
    if (sam_aln1.is_proper && sam_aln2.is_proper) {
        f1 |= PROPER_PAIR;
        f2 |= PROPER_PAIR;
    }

    std::string output_read1;
    output_read1 = record1.seq;
    std::string output_read2;
    output_read2 = record2.seq;
    if (sam_aln1.is_rc) {
        f1 |= REVERSE;
        f2 |= MREVERSE;
        output_read1 = read1_rc;
    }
    if (sam_aln2.is_rc) {
        f1 |= MREVERSE;
        f2 |= REVERSE;
        output_read2 = read2_rc;
    }

    std::string mate_name1;
    std::string mate_name2;
    if (sam_aln1.ref_id == sam_aln2.ref_id){
        mate_name1 = "=";
        mate_name2 = "=";
    } else{
        mate_name1 = references.names[sam_aln1.ref_id];
        mate_name2 = references.names[sam_aln2.ref_id];
    }

//    if ( (sam_aln1.is_unaligned) && (sam_aln2.is_unaligned) ){
//        f1 = PAIRED | UNMAP | MUNMAP;
//        f2 = PAIRED | UNMAP | MUNMAP;
//        m1_chr = "*";
//        m2_chr = "*";
//        sam_aln1.cigar = "*";
//        sam_aln2.cigar = "*";
//    } else if (sam_aln1.is_unaligned){
//        f1 = PAIRED | UNMAP;
//        m1_chr = "*";
//        sam_aln1.cigar = "*";
//        f2 |= MUNMAP;
//        f2 -= 32;
//    } else if (sam_aln2.is_unaligned){
//        f2 = 5;
//        m2_chr = "*";
//        sam_aln2.cigar = "*";
//        f1 |= MUNMAP;
//        f1 -= 32;
//    }
    std::string ref1 = references.names[sam_aln1.ref_id];
    std::string ref2 = references.names[sam_aln2.ref_id];
    int ed1 = sam_aln1.ed;
    int ed2 = sam_aln2.ed;

    if (sam_aln1.is_unaligned && sam_aln2.is_unaligned){
        f1 |= UNMAP;
        f1 |= MUNMAP;
        f2 |= UNMAP;
        f2 |= MUNMAP;
        sam_aln1.ref_start = 0;
        sam_aln2.ref_start = 0;
        template_len1 = 0;
        ref1 = "*";
        ref2 = "*";
        f1 |= (0u << 4);  // REVERSE
        f1 |= (0u << 5);  // MREVERSE
        f2 |= (0u << 4);  // REVERSE
        f2 |= (0u << 5);  // MREVERSE
        ed1 = 0;
        ed2 = 0;
        mapq1 = SAM_UNMAPPED_MAPQ;
        mapq2 = SAM_UNMAPPED_MAPQ;
    } else if (sam_aln1.is_unaligned){
        f1 |= UNMAP;
        f1 |= (0u << 4);  // REVERSE
        f2 |= MUNMAP;
        sam_aln1.ref_start = sam_aln2.ref_start;
        template_len1 = 0;
        ed1 = 0;
        mapq1 = SAM_UNMAPPED_MAPQ;
    } else if (sam_aln2.is_unaligned){
        f2 |= UNMAP;
        f2 |= (0u << 4);  // REVERSE
        f1 |= MUNMAP;
        sam_aln2.ref_start = sam_aln1.ref_start;
        template_len1 = 0;
        ed2 = 0;
        mapq2 = SAM_UNMAPPED_MAPQ;
    }

    add_one(record1, f1, ref1, sam_aln1, mapq1, mate_name2, sam_aln2.ref_start, template_len1, output_read1, ed1);
    add_one(record2, f2, ref2, sam_aln2, mapq2, mate_name1, sam_aln1.ref_start, -template_len1, output_read2, ed2);
}
