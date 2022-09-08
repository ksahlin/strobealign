#include <ostream>
#include "aln.hpp"

void Sam::unmapped(std::string& name, std::string& sequence, std::string& qual) {
    sam_string.append(name);
    sam_string.append("\t4\t*\t0\t255\t*\t*\t0\t0\t");
    sam_string.append(sequence);
    sam_string.append("\t*\n");
}


void Sam::add(const alignment& sam_aln, const std::string& sequence, const std::string& sequence_rc, const std::string& query_acc, const std::string& qual) {
        int flags = 0;
        std::string output_read;
        if (sam_aln.is_rc) {
            flags |= REVERSE;
            output_read = sequence_rc;
        } else {
            output_read = sequence;
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

//        sam_string.append("\t*\tNM:i:");
//        sam_string.append(std::to_string(sam_aln.ed));
//        sam_string.append("\n");



}
