#include <ostream>
#include "aln.hpp"

void Sam::add_unmapped(std::string& name, std::string& sequence, std::string& qual) {
    sam_string.append(name);
    sam_string.append("\t4\t*\t0\t255\t*\t*\t0\t0\t");
    sam_string.append(sequence);
    sam_string.append("\t*\n");
}

void Sam::add(
    const alignment& sam_aln,
    const std::string& sequence,
    const std::string& sequence_rc,
    const std::string& query_acc,
    const std::string& qual,
    bool is_secondary
) {
    int flags = 0;
    std::string output_read;
    if (sam_aln.is_rc) {
        flags |= REVERSE;
        output_read = sequence_rc;
    } else {
        output_read = sequence;
    }
    if (is_secondary) {
        flags |= SECONDARY;
    }
    sam_string.append(query_acc);
    sam_string.append("\t");
    sam_string.append(std::to_string(flags));
    sam_string.append("\t");
    sam_string.append(acc_map[sam_aln.ref_id]);
    sam_string.append("\t");
    sam_string.append(std::to_string(sam_aln.ref_start));
    sam_string.append("\t");
    sam_string.append(std::to_string(sam_aln.mapq));
    sam_string.append("\t");
    sam_string.append(sam_aln.cigar);
    sam_string.append("\t*\t0\t0\t");
    if (!sam_aln.is_unaligned) {
        sam_string.append(output_read);
        sam_string.append("\t");
        if (sam_aln.is_rc){
            auto qual_rev = qual;
            std::reverse(qual_rev.begin(), qual_rev.end()); // reverse
            sam_string.append(qual_rev);
        } else {
            sam_string.append(qual);
        }
        sam_string.append("\t");
        sam_string.append("NM:i:");
        sam_string.append(std::to_string(sam_aln.ed));
        sam_string.append("\t");
        sam_string.append("AS:i:");
        sam_string.append(std::to_string((int) sam_aln.aln_score));
    } else {
        sam_string.append(sequence);
        sam_string.append("\t");
        sam_string.append(qual);
    }
    sam_string.append("\n");
}

void Sam::add_unmapped_pair(KSeq& r1, KSeq& r2) {
    sam_string.append(r1.name);
    sam_string.append("\t77\t*\t0\t255\t*\t*\t0\t0\t");
    sam_string.append(r1.seq);
    sam_string.append("\t*\n");
    sam_string.append(r2.name);
    sam_string.append("\t141\t*\t0\t255\t*\t*\t0\t0\t");
    sam_string.append(r2.seq);
    sam_string.append("\t*\n");
}

void Sam::add_pair(alignment &sam_aln1, alignment &sam_aln2, const std::string &read1, const std::string &read2, const std::string &read1_rc, const std::string &read2_rc, const std::string &query_acc1, const std::string &query_acc2, int &mapq1, int &mapq2, float &mu, float &sigma, int read_len, bool is_primary, const std::string &qual1, const std::string &qual2) {
    int f1 = PAIRED;
    int f2 = PAIRED; // template having multiple segments in sequencing
//    std::cerr <<  query_acc1 << " HERE 1: " << sam_aln1.is_unaligned << " " << query_acc2 << " HERE 2: " << sam_aln2.is_unaligned  << std::endl;
    if (!is_primary){
        f1 |= SECONDARY;
        f2 |= SECONDARY;
    }

    // Commented lines below because we do not longer mark a read as not proper just because of the non-matching hash
    // Proper or non proper reads are further below only decided based on the expected distance and relative orientation they align to
//    if (sam_aln1.ed < 5){ // Flag alignments previously deemed as 'not proper' (based on matching strobemer hash ) to proper because of small ed
//        sam_aln1.not_proper = false;
//    }
//    if (sam_aln2.ed < 5){ // Flag alignments previously deemed as 'not proper' (based on matching strobemer hash ) to proper because of small ed
//        sam_aln2.not_proper = false;
//    }

    int d, template_len1, template_len2;
    if (sam_aln1.ref_start < sam_aln2.ref_start){
        d = sam_aln2.ref_start - sam_aln1.ref_start;
        template_len1 = - d - sam_aln2.aln_length;
        template_len2 = d + sam_aln2.aln_length;
//        std::cerr << "case1: " <<  query_acc1 << "HERE 1: " << sam_aln1.ref_start << " " << query_acc2 << " HERE 2: " << sam_aln2.ref_start << "  " << d <<  " sam_aln2.aln_length " << sam_aln2.aln_length << " sam_aln1.aln_length " << sam_aln1.aln_length << std::endl;
    }
    else{
        d = sam_aln1.ref_start - sam_aln2.ref_start;
        template_len1 = d + sam_aln1.aln_length;
        template_len2 = - d - sam_aln1.aln_length;
//        std::cerr << "case2: " <<   query_acc1 << "HERE 1: " << sam_aln1.ref_start << " " << query_acc2 << " HERE 2: " << sam_aln2.ref_start << "  " << d <<  " sam_aln2.aln_length " << sam_aln2.aln_length << " sam_aln1.aln_length " << sam_aln1.aln_length << std::endl;
    }
    //    int d = sam_aln1.ref_start < sam_aln2.ref_start ? sam_aln2.ref_start - sam_aln1.ref_start : sam_aln1.ref_start - sam_aln2.ref_start;

    bool both_aligned = !sam_aln1.is_unaligned && !sam_aln2.is_unaligned;
    int a = sam_aln2.ref_start - sam_aln1.ref_start;
    bool r1_r2 = !sam_aln1.is_rc && sam_aln2.is_rc && a >= 0; // r1 ---> <---- r2
    bool r2_r1 = !sam_aln2.is_rc && sam_aln1.is_rc && a <= 0 ; // r2 ---> <---- r1
    bool rel_orientation_good = r1_r2 || r2_r1;
    bool insert_good = d <= mu + 6*sigma;
    if ( both_aligned && insert_good && rel_orientation_good ){
        sam_aln1.not_proper = false;
        sam_aln2.not_proper = false;
    } else {
        sam_aln1.not_proper = true;
        sam_aln2.not_proper = true;
    }

//    if ( d > (mu + 6*sigma) ){ // Flag alignments as 'not proper' because of too large deviation in insert size
//            sam_aln1.not_proper = true;
//            sam_aln2.not_proper = true;
//        }

    if ( (!sam_aln1.not_proper) && (!sam_aln2.not_proper)){ // if both segements in pair are properly aligned
        f1 |= PROPER_PAIR;
        f2 |= PROPER_PAIR;
    }

    std::string output_read1;
    output_read1 = read1;
    std::string output_read2;
    output_read2 = read2;
    f1 |= READ1;
    f2 |= READ2;
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

    std::string m1_chr;
    std::string m2_chr;
    if (sam_aln1.ref_id == sam_aln2.ref_id){
        m1_chr = "=";
        m2_chr = "=";
    } else{
        m1_chr = acc_map[sam_aln1.ref_id];
        m2_chr = acc_map[sam_aln2.ref_id];
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
    std::string ref1 = acc_map[sam_aln1.ref_id];
    std::string ref2 = acc_map[sam_aln2.ref_id];
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
        template_len2 = 0;
        ref1 = "*";
        ref2 = "*";
        f1 |= (0u << 4);
        f1 |= (0u << 5);
        f2 |= (0u << 4);
        f2 |= (0u << 5);
        ed1 = 0;
        ed2 = 0;
        mapq1 = 255;
        mapq2 = 255;
    } else if (sam_aln1.is_unaligned){
        f1 |= (1u << 2);
        f1 |= (0u << 4);
        f2 |= (1u << 3);
        sam_aln1.ref_start = sam_aln2.ref_start;
        template_len1 = 0;
        template_len2 = 0;
        ed1 = 0;
        mapq1 = 255;
    } else if (sam_aln2.is_unaligned){
        f2 |= (1u << 2);
        f2 |= (0u << 4);
        f1 |= (1u << 3);
        sam_aln2.ref_start = sam_aln1.ref_start;
        template_len1 = 0;
        template_len2 = 0;
        ed2 = 0;
        mapq2 = 255;
    }

//    if ( (sam_aln1.ref_start == 0)){ // && !sam_aln1.is_unaligned  ){
//        std::cerr << "OMG1" << std::endl;
//        std::cerr << query_acc1 << std::endl;
//    }
//    if ( (sam_aln2.ref_start == 0)){ // && !sam_aln2.is_unaligned  ){
//        std::cerr << "OMG2" << std::endl;
//        std::cerr << query_acc2 << std::endl;
//    }

    sam_string.append(query_acc1);
    sam_string.append("\t");
    sam_string.append(std::to_string(f1));
    sam_string.append("\t");
    sam_string.append(ref1);
    sam_string.append("\t");
    sam_string.append(std::to_string(sam_aln1.ref_start));
    sam_string.append("\t");
    sam_string.append(std::to_string(mapq1));
    sam_string.append("\t");
    sam_string.append(sam_aln1.cigar);
    sam_string.append("\t");
    sam_string.append(m2_chr);
    sam_string.append("\t");
    sam_string.append(std::to_string(sam_aln2.ref_start));
    sam_string.append("\t");
    sam_string.append(std::to_string(template_len1));
    sam_string.append("\t");
//    sam_string.append("\t*\t0\t0\t");
    if (!sam_aln1.is_unaligned) {
        sam_string.append(output_read1);
        sam_string.append("\t");
        if (sam_aln1.is_rc){
            auto qual_rev = qual1;
            std::reverse(qual_rev.begin(), qual_rev.end());
            sam_string.append(qual_rev);
        } else {
            sam_string.append(qual1);
        }
        sam_string.append("\t");
        sam_string.append("NM:i:");
        sam_string.append(std::to_string(ed1));
        sam_string.append("\t");
        sam_string.append("AS:i:");
        sam_string.append(std::to_string((int) sam_aln1.aln_score));
    } else {
        sam_string.append(read1);
        sam_string.append("\t");
        sam_string.append(qual1);
    }
    sam_string.append("\n");

    sam_string.append(query_acc2);
    sam_string.append("\t");
    sam_string.append(std::to_string(f2));
    sam_string.append("\t");
    sam_string.append(ref2);
    sam_string.append("\t");
    sam_string.append(std::to_string(sam_aln2.ref_start));
    sam_string.append("\t");
    sam_string.append(std::to_string(mapq2));
    sam_string.append("\t");
    sam_string.append(sam_aln2.cigar);
    sam_string.append("\t");
    sam_string.append(m1_chr);
    sam_string.append("\t");
    sam_string.append(std::to_string(sam_aln1.ref_start));
    sam_string.append("\t");
    sam_string.append(std::to_string(template_len2));
    sam_string.append("\t");
//    sam_string.append("\t*\t0\t0\t");
    if (!sam_aln2.is_unaligned) {
        sam_string.append(output_read2);
        sam_string.append("\t");
        if (sam_aln2.is_rc){
            auto qual2_rev = qual2;
            std::reverse(qual2_rev.begin(), qual2_rev.end()); // reverse
            sam_string.append(qual2_rev);
        } else {
            sam_string.append(qual2);
        }
        sam_string.append("\t");
        sam_string.append("NM:i:");
        sam_string.append(std::to_string(ed2));
        sam_string.append("\t");
        sam_string.append("AS:i:");
        sam_string.append(std::to_string((int) sam_aln2.aln_score));
    } else {
        sam_string.append(read2);
        sam_string.append("\t");
        sam_string.append(qual2);
    }
    sam_string.append("\n");
}
