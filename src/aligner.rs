//! Low-level alignment functions

use std::cell::Cell;
use crate::cigar::{Cigar, CigarOperation};

pub struct Scores {
    // match is a score, the others are penalties (all are nonnegative)
    pub match_: u32,
    pub mismatch: u32,
    pub gap_open: u32,
    pub gap_extend: u32,
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
    // ssw_aligner: StripedSmithWaterman::Aligner,
    // ssw_filter: StripedSmithWaterman::Filter,
}

impl Aligner {
    /*Aligner(Scores parameters)
        : parameters(parameters)
        , ssw_aligner(StripedSmithWaterman::Aligner(parameters.match_, parameters.mismatch, parameters.gap_open, parameters.gap_extend))
    { }*/

    pub fn new(scores: Scores) -> Self {
        Aligner {
            scores,
            call_count: Cell::new(0usize),
        }
    }

    fn call_count(&self) -> usize {
        self.call_count.get()
    }

    pub fn align(&self, query: &[u8], refseq: &[u8]) -> Option<AlignmentInfo> {
        self.call_count.set(self.call_count.get() + 1);
        /*
        // AlignmentInfo aln;
        // int32_t maskLen = query.len() / 2;
        maskLen = std::max(maskLen, 15);
        if refseq.len() > 2000 {
            //        std::cerr << "ALIGNMENT TO REF LONGER THAN 2000bp - REPORT TO DEVELOPER. Happened for read: " <<  query << " ref len:" << ref.length() << std::endl;
            return None;
        }

        let alignment_ssw: StripedSmithWaterman::Alignment;

        // query must be NULL-terminated
        let flag = ssw_aligner.Align(query.c_str(), refseq.c_str(), refseq.len(), self.ssw_filter, &alignment_ssw, maskLen);
        if flag != 0 {
            return None;
        }

        let mut aln = AlignmentInfo { };
        aln.edit_distance = alignment_ssw.mismatches;
        aln.cigar = Cigar::parse(alignment_ssw.cigar);
        aln.score = alignment_ssw.sw_score;
        aln.ref_start = alignment_ssw.ref_begin;
        // end positions are off by 1 in SSW
        aln.ref_end = alignment_ssw.ref_end + 1;
        aln.query_start = alignment_ssw.query_begin;
        aln.query_end = alignment_ssw.query_end + 1;

        // Try to extend to beginning of the query to get an end bonus
        let mut qstart = aln.query_start;
        let mut rstart = aln.ref_start;
        let mut score = aln.score;
        let mut edits = aln.edit_distance;
        let mut front_cigar = Cigar::new();
        while qstart > 0 && rstart > 0 {
            qstart -= 1;
            rstart -= 1;
            if query[qstart] == refseq[rstart] {
                score += self.scores.match_;
                front_cigar.push(CigarOperation::Eq, 1);
            } else {
                score -= self.scores.mismatch;
                front_cigar.push(CigarOperation::X, 1);
                edits += 1;
            }
        }
        if qstart == 0 && score + self.scores.end_bonus > aln.score {
            if aln.query_start > 0 {
                // TODO assert_eq!(aln.cigar.m_ops[0] & 0xF, CIGAR_SOFTCLIP);
                aln.cigar.m_ops.erase(aln.cigar.m_ops.begin());  // remove soft clipping
                front_cigar.reverse();
                front_cigar += aln.cigar;
                aln.cigar = front_cigar;
            }
            aln.query_start = 0;
            aln.ref_start = rstart;
            aln.score = score + self.scores.end_bonus;
            aln.edit_distance = edits;
        }

        // Try to extend to end of query to get an end bonus
        let mut qend = aln.query_end;
        let mut rend = aln.ref_end;
        score = aln.score;
        edits = aln.edit_distance;
        let mut back_cigar = Cigar::new();
        while qend < query.len() && rend < refseq.len() {
            if query[qend] == refseq[rend] {
                score += self.scores.match_;
                back_cigar.push(CigarOperation::Eq, 1);
            } else {
                score -= self.scores.mismatch;
                back_cigar.push(CigarOperation::X, 1);
                edits += 1;
            }
            qend += 1;
            rend += 1;
        }
        if qend == query.len() && score + self.scores.end_bonus > aln.score {
            if aln.query_end < query.len() {
                // TODO assert((aln.cigar.m_ops[aln.cigar.m_ops.size() - 1] & 0xf) == CIGAR_SOFTCLIP);
                aln.cigar.m_ops.pop_back();
                aln.cigar += back_cigar;
            }
            aln.query_end = query.len();
            aln.ref_end = rend;
            aln.score = score + self.scores.end_bonus;
            aln.edit_distance = edits;
        }

        Some(aln)
        */
        None
    }
}

pub fn hamming_distance(s: &[u8], t: &[u8]) -> Option<u32> {
    if s.len() != t.len() {
        return None;
    }

    let mismatches =
        s.iter().zip(t).map(
            |(x, y)| if x != y { 1 } else { 0 }
        ).sum();

    Some(mismatches)
}

/// Find highest-scoring segment between reference and query assuming only matches
/// and mismatches are allowed.
///
/// The end_bonus is added to the score if the segment extends until the end
/// of the query, once for each end.
fn highest_scoring_segment(query: &[u8], refseq: &[u8], match_: u32, mismatch: u32, end_bonus: u32) -> (usize, usize, u32) {
    let n = query.len();

    let mut start = 0; // start of the current segment
    let mut score = end_bonus; // accumulated score so far in the current segment

    let mut best_start = 0;
    let mut best_end = 0;
    let mut best_score = 0;
    for i in 0..n {
        if query[i] == refseq[i] {
            score += match_;
        } else {
            score -= mismatch;
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
    if score + end_bonus > best_score {
        best_score = score + end_bonus;
        best_end = query.len();
        best_start = start;
    }

    (best_start, best_end, best_score)
}

pub fn hamming_align(query: &[u8], refseq: &[u8], match_: u32, mismatch: u32, end_bonus: u32) -> Option<AlignmentInfo> {
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
        mismatches += if is_match { 0 } else { 1 };
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
