use std::fmt::Display;

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

#[derive(Default, Debug)]
pub struct PafRecord {
    pub query_name: String,
    pub query_length: u64,
    pub query_start: u64,
    pub query_end: u64,
    pub is_revcomp: bool,
    pub target_name: String,
    pub target_length: u64,
    pub target_start: u64,
    pub target_end: u64,
    pub n_matches: u64,
    pub alignment_length: u64,
    pub mapping_quality: Option<u8>,
}

impl Display for PafRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_name,
            self.query_length,
            self.query_start,
            self.query_end,
            if self.is_revcomp { "-" } else { "+" },
            self.target_name,
            self.target_length,
            self.target_start,
            self.target_end,
            self.n_matches,
            self.alignment_length,
            self.mapping_quality.unwrap_or(255),
        )
    }
}
