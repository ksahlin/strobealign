use std::cell::Cell;
use crate::cigar::{Cigar, CigarOperation};
use crate::ssw::SswAligner;

pub struct Scores {
    // match is a score, the others are penalties
    pub match_: u8,
    pub mismatch: u8,
    pub gap_open: u8,
    pub gap_extend: u8,
    pub end_bonus: u8,
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
    pub edit_distance: u32,
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

    fn call_count(&self) -> usize {
        self.call_count.get()
    }

    pub fn align(&self, query: &[u8], refseq: &[u8]) -> Option<AlignmentInfo> {
        self.call_count.set(self.call_count.get() + 1);

        if refseq.len() > 2000 {
            // TODO
            return None;
        }

        let mut alignment: AlignmentInfo = self.ssw_aligner.align(query, refseq)?.into();

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
fn highest_scoring_segment(query: &[u8], refseq: &[u8], match_: u8, mismatch: u8, end_bonus: u8) -> (usize, usize, u32) {
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

pub fn hamming_align(query: &[u8], refseq: &[u8], match_: u8, mismatch: u8, end_bonus: u8) -> Option<AlignmentInfo> {
    if query.len() != refseq.len() {
        return None;
    }

    let (segment_start, segment_end, score) = highest_scoring_segment(query, refseq, match_, mismatch, end_bonus);

    let mut cigar = Cigar::new();
    if segment_start > 0 {
        cigar.push(CigarOperation::Softclip, segment_start);
    }

    // Create CIGAR string and count mismatches
    let mut counter = 0;
    let mut prev_is_match = false;
    let mut mismatches = 0;
    let mut first = true;
    for i in segment_start..segment_end {
        let is_match = query[i] == refseq[i];
        mismatches += u32::from(!is_match);
        if !first && is_match != prev_is_match {
            cigar.push(if prev_is_match { CigarOperation::Eq } else { CigarOperation::X }, counter);
            counter = 0;
        }
        counter += 1;
        prev_is_match = is_match;
        first = false;
    }
    if !first {
        cigar.push(if prev_is_match { CigarOperation::Eq } else { CigarOperation::X }, counter);
    }

    let soft_right = query.len() - segment_end;
    if soft_right > 0 {
        cigar.push(CigarOperation::Softclip, soft_right);
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
