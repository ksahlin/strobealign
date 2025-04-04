use std::cell::Cell;
use crate::cigar::{Cigar, CigarOperation};
use crate::ssw::SswAligner;

#[derive(Debug, Clone)]
pub struct Scores {
    // match is a score, the others are penalties
    pub match_: u8,
    pub mismatch: u8,
    pub gap_open: u8,
    pub gap_extend: u8,
    pub end_bonus: u32,
}

impl Default for Scores {
    fn default() -> Scores {
        Scores {
            match_: 2,
            mismatch: 8,
            gap_open: 12,
            gap_extend: 1,
            end_bonus: 10,
        }
    }
}

#[derive(Debug)]
pub struct AlignmentInfo {
    pub cigar: Cigar,
    pub edit_distance: usize,
    pub ref_start: usize,
    pub ref_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub score: u32,
}

impl AlignmentInfo {
    pub fn ref_span(&self) -> usize {
        self.ref_end - self.ref_start
    }
}

#[derive(Clone)]
pub struct Aligner {
    pub scores: Scores, // TODO should not be pub?
    call_count: Cell<usize>,
    ssw_aligner: SswAligner,
}

impl Aligner {
    pub fn new(scores: Scores) -> Self {
        let ssw_aligner = SswAligner::new(scores.match_, scores.mismatch, scores.gap_open, scores.gap_extend);
        Aligner {
            scores,
            ssw_aligner,
            call_count: Cell::new(0usize),
        }
    }

    pub fn call_count(&self) -> usize {
        self.call_count.get()
    }

    /// The returned CIGAR string does not include soft clipped bases
    /// Use query_start and query_end to compute the number of clipped bases
    pub fn align(&self, query: &[u8], refseq: &[u8]) -> Option<AlignmentInfo> {
        self.call_count.set(self.call_count.get() + 1);

        if refseq.len() > 2000 {
            // TODO
            return None;
        }

        let mut alignment = self.ssw_aligner.align(query, refseq)?;

        // Try to extend to beginning of the query to get an end bonus
        let mut qstart = alignment.query_start;
        let mut rstart = alignment.ref_start;
        let mut score = alignment.score as i32;
        let mut edits = alignment.edit_distance;
        let mut front_cigar = Cigar::new();
        while qstart > 0 && rstart > 0 {
            qstart -= 1;
            rstart -= 1;
            if query[qstart] == refseq[rstart] {
                score += self.scores.match_ as i32;
                front_cigar.push(CigarOperation::Eq, 1);
            } else {
                score -= self.scores.mismatch as i32;
                front_cigar.push(CigarOperation::X, 1);
                edits += 1;
            }
        }
        let score_with_end_bonus = score + self.scores.end_bonus as i32;
        if qstart == 0 && score_with_end_bonus > alignment.score as i32 {
            if alignment.query_start > 0 {
                let mut cigar = front_cigar.reversed();
                cigar.extend(&alignment.cigar);
                alignment.cigar = cigar;
            }
            alignment.query_start = 0;
            alignment.ref_start = rstart;
            alignment.score = score_with_end_bonus as u32;
            alignment.edit_distance = edits;
        }

        // Try to extend to end of query to get an end bonus
        let mut qend = alignment.query_end;
        let mut rend = alignment.ref_end;
        score = alignment.score as i32;
        edits = alignment.edit_distance;
        let mut back_cigar = Cigar::new();
        while qend < query.len() && rend < refseq.len() {
            if query[qend] == refseq[rend] {
                score += self.scores.match_ as i32;
                back_cigar.push(CigarOperation::Eq, 1);
            } else {
                score -= self.scores.mismatch as i32;
                back_cigar.push(CigarOperation::X, 1);
                edits += 1;
            }
            qend += 1;
            rend += 1;
        }
        let score_with_end_bonus = score + self.scores.end_bonus as i32;
        if qend == query.len() && score_with_end_bonus > alignment.score as i32 {
            if alignment.query_end < query.len() {
                alignment.cigar.extend(&back_cigar);
            }
            alignment.query_end = query.len();
            alignment.ref_end = rend;
            alignment.score = score_with_end_bonus as u32;
            alignment.edit_distance = edits;
        }

        Some(alignment)
    }
}

pub fn hamming_distance(s: &[u8], t: &[u8]) -> Option<u32> {
    if s.len() != t.len() {
        return None;
    }

    let mismatches =
        s.iter().zip(t).map(
            |(x, y)| u32::from(x != y)
        ).sum();

    Some(mismatches)
}

/// Find highest-scoring segment between reference and query assuming only matches
/// and mismatches are allowed.
///
/// The end_bonus is added to the score if the segment extends until the end
/// of the query, once for each end.
fn highest_scoring_segment(query: &[u8], refseq: &[u8], match_: u8, mismatch: u8, end_bonus: u32) -> (usize, usize, u32) {
    let n = query.len();

    let mut start = 0; // start of the current segment
    let mut score = end_bonus as i32; // accumulated score so far in the current segment

    let mut best_start = 0;
    let mut best_end = 0;
    let mut best_score = 0;
    for i in 0..n {
        if query[i] == refseq[i] {
            score += match_ as i32;
        } else {
            score -= mismatch as i32;
        }
        if score < 0 {
            start = i + 1;
            score = 0;
        }
        if score > best_score {
            best_start = start;
            best_score = score;
            best_end = i + 1;
        }
    }
    if score + end_bonus as i32 > best_score {
        best_score = score + end_bonus as i32;
        best_end = query.len();
        best_start = start;
    }

    (best_start, best_end, best_score as u32)
}

/// The returned CIGAR string does not include soft clipped bases
/// Use query_start and query_end to compute the number of clipped bases
pub fn hamming_align(query: &[u8], refseq: &[u8], match_: u8, mismatch: u8, end_bonus: u32) -> Option<AlignmentInfo> {
    if query.len() != refseq.len() {
        return None;
    }
    let (segment_start, segment_end, score) = highest_scoring_segment(query, refseq, match_, mismatch, end_bonus);

    // Create CIGAR string and count mismatches
    let mut cigar = Cigar::new();
    let mut mismatches = 0;
    let query_segment = &query[segment_start..segment_end];
    let ref_segment = &refseq[segment_start..segment_end];
    for (q_c, r_c) in query_segment.iter().zip(ref_segment.iter()) {
        if q_c != r_c {
            mismatches += 1;
            cigar.push(CigarOperation::X, 1);
        } else {
            cigar.push(CigarOperation::Eq, 1);
        }
    }

    Some(AlignmentInfo {
        cigar,
        score,
        edit_distance: mismatches,
        ref_start: segment_start,
        ref_end: segment_end,
        query_start: segment_start,
        query_end: segment_end,
    })
}

#[cfg(test)]
mod test {
    use crate::aligner::{hamming_align, highest_scoring_segment, Aligner, Scores};

    #[test]
    fn test_ssw_align_no_result() {
        let aligner = Aligner::new(Scores::default());
        let query = b"TCTCTCCCTCTCTCTCTCTCCCTCCCTCTCTCTCCCTCTCTCTCTCTCTCTCCCTCCCTT";
        let refseq = b"GAGGGAGAGAGAGAGAGGGAGAGAGAGAGAGAG";
        let info = aligner.align(query, refseq);
        assert!(info.is_none());
    }

    #[test]
    fn test_call_count() {
        let aligner = Aligner::new(Scores::default());
        assert_eq!(aligner.call_count(), 0);
        let refseq = b"AAAAAACCCCCGGGGG";
        let query = b"CCCCC";
        let info = aligner.align(query, refseq);
        assert!(info.is_some());
        assert_eq!(aligner.call_count(), 1);
        aligner.align(query, refseq);
        assert_eq!(aligner.call_count(), 2);
    }

    #[test]
    fn test_highest_scoring_segment() {
        let (start, end, score) = highest_scoring_segment(b"", b"", 5, 7, 0);
        assert_eq!(start, 0);
        assert_eq!(end, 0);
        
        let (start, end, score) = highest_scoring_segment(b"AAAAAAAAAA", b"AAAAAATTTT", 5, 7, 0);
        assert_eq!(start, 0);
        assert_eq!(end, 6);
        let (start, end, score) = highest_scoring_segment(b"AAAAAAAAAA", b"TTTTAAAAAA", 5, 7, 0);
        assert_eq!(start, 4);
        assert_eq!(end, 10);
        
        assert_eq!(highest_scoring_segment(b"AAAAAAAAAA", b"AAAAAATTTT", 5, 7, 0), (0, 6, 30));
        assert_eq!(highest_scoring_segment(b"AAAAAAAAAA", b"TTTTAAAAAA", 5, 7, 0), (4, 10, 30));
        assert_eq!(highest_scoring_segment(b"AAAAAAAAAA", b"TTAAAAAATT", 5, 7, 0), (2, 8, 30));
        assert_eq!(highest_scoring_segment(b"AAAAAAAAAAAAAAA", b"TAAAAAATTTAAAAT", 5, 7, 0), (1, 7, 30));
    }

    #[test]
    fn test_highest_scoring_segment_with_soft_clipping() {
        assert_eq!(highest_scoring_segment(b"", b"", 2, 4, 5), (0, 0, 10));
        assert_eq!(highest_scoring_segment(b"TAAT", b"TAAA", 2, 4, 5), (0, 4, 3 * 2 - 4 + 10));
        assert_eq!(highest_scoring_segment(b"AAA", b"AAA", 2, 4, 5), (0, 3, 3 * 2 + 10));
        assert_eq!(highest_scoring_segment(b"TAAT", b"AAAA", 2, 4, 5), (0, 4, 10 + 2 * 2 - 2 * 4));
        assert_eq!(highest_scoring_segment(b"ATAATA", b"AAAAAA", 2, 4, 5), (0, 6, 4 * 2 - 2 * 4 + 10));
        assert_eq!(highest_scoring_segment(b"TTAATA", b"AAAAAA", 2, 4, 5), (2, 6, 3 * 2 - 1 * 4 + 5));
    }
    
    #[test]
    fn test_hamming_align_empty_sequences() {
        let info = hamming_align(b"", b"", 7, 5, 0).unwrap();
        assert!(info.cigar.is_empty());
        assert_eq!(info.edit_distance, 0);
        assert_eq!(info.score, 0);
        assert_eq!(info.ref_start, 0);
        assert_eq!(info.ref_span(), 0);
        assert_eq!(info.query_start, 0);
        assert_eq!(info.query_end, 0);
    }

    #[test]
    fn test_hamming_align_one_mismatch() {
        let info = hamming_align(
            b"AAXGGG",
            b"AAYGGG",
            1, 1, 0
        ).unwrap();
        assert_eq!(info.cigar.to_string(), "2=1X3=");
        assert_eq!(info.edit_distance, 1);
        assert_eq!(info.score, 4);
        assert_eq!(info.ref_start, 0);
        assert_eq!(info.ref_span(), 6);
        assert_eq!(info.query_start, 0);
        assert_eq!(info.query_end, 6);
    }

    #[test]
    fn test_hamming_align_soft_clipped() {
        let info = hamming_align(
            b"AXGGG",
            b"AYGGG",
            1, 4, 0
        ).unwrap();
        assert_eq!(info.cigar.to_string(), "3=");
        assert_eq!(info.edit_distance, 0);
        assert_eq!(info.score, 3);
        assert_eq!(info.ref_start, 2);
        assert_eq!(info.ref_span(), 3);
        assert_eq!(info.query_start, 2);
        assert_eq!(info.query_end, 5);
    }

    #[test]
    fn test_hamming_align_wildcard() {
        let info = hamming_align(
            b"NAACCG",
            b"TAACCG",
            3, 7, 0
        ).unwrap();
        assert_eq!(info.cigar.to_string(), "5=");
        assert_eq!(info.edit_distance, 0);
        assert_eq!(info.score, 5 * 3);
        assert_eq!(info.ref_start, 1);
        assert_eq!(info.ref_span(), 5);
        assert_eq!(info.query_start, 1);
        assert_eq!(info.query_end, 6);
    }

    #[test]
    fn test_hamming_align_wildcard_at_end() {
        let info = hamming_align(
            b"AACCGN",
            b"AACCGT",
            3, 7, 0
        ).unwrap();
        assert_eq!(info.cigar.to_string(), "5=");
        assert_eq!(info.edit_distance, 0);
        assert_eq!(info.score, 5 * 3);
        assert_eq!(info.ref_start, 0);
        assert_eq!(info.ref_span(), 5);
        assert_eq!(info.query_start, 0);
        assert_eq!(info.query_end, 5);
    }

    #[test]
    fn test_hamming_align_both_ends_soft_clipped() {
        // negative total score, soft clipping on both ends
        let info = hamming_align(
            b"NAAAAAAAAAAAAAA",
            b"TAAAATTTTTTTTTT",
            3, 7, 0
        ).unwrap();
        assert_eq!(info.cigar.to_string(), "4=");
        assert_eq!(info.edit_distance, 0);
        assert_eq!(info.score, 4 * 3);
        assert_eq!(info.ref_start, 1);
        assert_eq!(info.ref_span(), 4);
        assert_eq!(info.query_start, 1);
        assert_eq!(info.query_end, 5);
    }

    #[test]
    fn test_hamming_align_long_soft_clipping() {
        let info = hamming_align(
            b"NAAAAAAAAAAAAAA",
            b"TAAAATTTAAAAAAT",
            3, 7, 0
        ).unwrap();
        assert_eq!(info.cigar.to_string(), "6=");
        assert_eq!(info.edit_distance, 0);
        assert_eq!(info.score, 6 * 3);
        assert_eq!(info.ref_start, 8);
        assert_eq!(info.ref_span(), 6);
        assert_eq!(info.query_start, 8);
        assert_eq!(info.query_end, 14);
    }

    #[test]
    fn test_hamming_align_short_soft_clipping() {
        let info = hamming_align(
            b"AAAAAAAAAAAAAAA",
            b"TAAAAAATTTAAAAT",
            3, 7, 0
        ).unwrap();
        assert_eq!(info.cigar.to_string(), "6=");
        assert_eq!(info.edit_distance, 0);
        assert_eq!(info.score, 6 * 3);
        assert_eq!(info.ref_start, 1);
        assert_eq!(info.ref_span(), 6);
        assert_eq!(info.query_start, 1);
        assert_eq!(info.query_end, 7);
    }
}
