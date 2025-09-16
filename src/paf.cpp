#include "paf.hpp"

/* PAF columns (see https://github.com/lh3/miniasm/blob/master/PAF.md):
 * 1 query name
 * 2 query length
 * 3 query start (0-based)
 * 4 query end
 * 5 relative strand (+ or -)
 * 6 target name
 * 7 target length
 * 8 target start
 * 9 target end
 * 10 no. of matches
 * 11 alignment block length
 * 12 mapping quality (0-255; 255 for missing)
 */
void output_hits_paf_PE(std::string &paf_output, const Chain &n, const std::string &query_name, const References& references, int read_len, uint8_t mapq) {
    paf_output.append(query_name);
    paf_output.append("\t");
    paf_output.append(std::to_string(read_len));
    paf_output.append("\t");
    paf_output.append(std::to_string(n.query_start));
    paf_output.append("\t");
    paf_output.append(std::to_string(n.query_end));
    paf_output.append("\t");
    paf_output.append(n.is_revcomp ? "-" : "+");
    paf_output.append("\t");
    paf_output.append(references.names[n.ref_id]);
    paf_output.append("\t");
    paf_output.append(std::to_string(references.lengths[n.ref_id]));
    paf_output.append("\t");
    paf_output.append(std::to_string(n.ref_start));
    paf_output.append("\t");
    paf_output.append(std::to_string(n.ref_end));
    paf_output.append("\t");
    paf_output.append(std::to_string(n.anchors.size()));
    paf_output.append("\t");
    paf_output.append(std::to_string(n.ref_end - n.ref_start));
    paf_output.append("\t");
    paf_output.append(std::to_string(mapq));
    paf_output.append("\n");
}


void output_hits_paf(std::string &paf_output, const std::vector<Chain> &chains, const std::string& query_name, const References& references, int read_len, uint8_t mapq) {
    // Output results
    if (chains.empty()) {
        return;
    }
    // Only output single best hit based on: number of randstrobe-matches times span of the merged match.
    Chain n = chains[0];
    output_hits_paf_PE(paf_output, n, query_name, references, read_len, mapq);
}
