use block_aligner::{
    cigar::{Cigar as BlockCigar, Operation},
    scan_block::{Block, PaddedBytes},
    scores::{AAMatrix, AAProfile, Gaps, NucMatrix, Profile},
};

use crate::{
    aligner::{hamming_align_global, AlignmentInfo, Scores},
    chainer::Anchor,
    cigar::{Cigar, CigarOperation},
};

#[derive(Clone)]
pub struct PiecewiseAligner {
    scores: Scores,
    k: usize,
    gaps: Gaps,
    xdrop: i32,
}

pub struct AlignmentResult {
    pub score: i32,
    pub query_start: usize,
    pub query_end: usize,
    pub ref_start: usize,
    pub ref_end: usize,
    pub cigar: Cigar,
}

impl PiecewiseAligner {
    pub fn new(scores: Scores, k: usize, xdrop: i32) -> Self {
        let gaps = Gaps {
            open: -(scores.gap_open as i8),
            extend: -(scores.gap_extend as i8),
        };
        PiecewiseAligner {
            scores,
            k,
            gaps,
            xdrop,
        }
    }

    // Returns a full global alignment using blockaligner
    fn global_alignment(&self, query: &[u8], reference: &[u8]) -> AlignmentResult {
        if query.is_empty() || reference.is_empty() {
            return AlignmentResult {
                score: 0,
                query_start: 0,
                query_end: 0,
                ref_start: 0,
                ref_end: 0,
                cigar: Cigar::new(),
            };
        }

        // Blockaligner block ranges, set to maximum values for full global alignment
        // Must be powers of 2
        let max_len = query.len().max(reference.len());
        let min_size = max_len.max(32).next_power_of_two();
        let max_size = max_len.max(128).next_power_of_two();

        // Create padded sequences
        let q_padded = PaddedBytes::from_bytes::<NucMatrix>(query, max_size);
        let r_padded = PaddedBytes::from_bytes::<NucMatrix>(reference, max_size);

        // Create simple scoring matrix (match/mismatch)
        let matrix = NucMatrix::new_simple(self.scores.match_ as i8, -(self.scores.mismatch as i8));

        // Make a global alignment call with traceback
        let mut block = Block::<true, false>::new(query.len(), reference.len(), max_size);
        block.align(
            &q_padded,
            &r_padded,
            &matrix,
            self.gaps,
            min_size..=max_size,
            0, // 0 x-drop threshold for global alignment
        );
        let res = block.res();

        // Retrieve CIGAR
        let mut cigar = BlockCigar::new(res.query_idx, res.reference_idx);
        block.trace().cigar_eq(
            &q_padded,
            &r_padded,
            res.query_idx,
            res.reference_idx,
            &mut cigar,
        );

        AlignmentResult {
            score: res.score,
            query_start: 0,
            query_end: res.query_idx,
            ref_start: 0,
            ref_end: res.reference_idx,
            cigar: build_cigar(&cigar),
        }
    }

    // Returns an x-drop alignment using blockaligner
    fn xdrop_alignment(&self, query: &[u8], reference: &[u8], reverse: bool) -> AlignmentResult {
        if query.is_empty() || reference.is_empty() {
            return AlignmentResult {
                score: 0,
                query_start: 0,
                query_end: 0,
                ref_start: 0,
                ref_end: 0,
                cigar: Cigar::new(),
            };
        }

        // Blockaligner block ranges, set to 20% and 50% of longest sequence
        // Must be powers of 2
        let max_len = query.len().max(reference.len());
        let min_size = (max_len / 5).max(32).next_power_of_two();
        let max_size = (max_len / 2).max(128).next_power_of_two();

        // In case we x-drop in reverse (aligning before the first anchor),
        // we exchange the query for the reversed reference and vice-versa
        let (query_used, ref_used) = if reverse {
            let mut query_rev = query.to_vec();
            query_rev.reverse();
            let mut ref_rev = reference.to_vec();
            ref_rev.reverse();
            (query_rev, ref_rev)
        } else {
            (query.to_vec(), reference.to_vec())
        };

        // Create padded sequences
        // blockaligner only has a AA profile, no nucleotide profile, so we have to use AA matrix
        let q_padded = PaddedBytes::from_bytes::<AAMatrix>(&query_used, max_size);
        let r_padded = PaddedBytes::from_bytes::<AAMatrix>(&ref_used, max_size);

        // Create profile to align the reference on the query with end bonus
        // We need to build a position-specific scoring matrix
        // blockaligner only has a AA profile, no nucleotide profile
        let mut profile =
            AAProfile::new(query_used.len(), max_size, -(self.scores.gap_extend as i8));

        // Set scores for each nucleotide combination
        for i in 1..=query_used.len() {
            let query_char = query_used[i - 1];
            for c in [b'A', b'C', b'G', b'T', b'N'] {
                if c == query_char {
                    profile.set(i, c, self.scores.match_ as i8);
                } else {
                    profile.set(i, c, -(self.scores.mismatch as i8));
                };
            }
        }

        // Set gap costs
        profile.set_all_gap_open_C(-(self.scores.gap_open as i8) - -(self.scores.gap_extend as i8));
        profile.set_all_gap_close_C(0);
        profile.set_all_gap_open_R(-(self.scores.gap_open as i8) - -(self.scores.gap_extend as i8));

        // Give a bonus score to the end
        let bonus_pos = query_used.len();
        for c in [b'A', b'C', b'G', b'T', b'N'] {
            let current_score = profile.get(bonus_pos, c);
            profile.set(bonus_pos, c, current_score + self.scores.end_bonus as i8);
        }

        // Make an x-drop alignment call with traceback
        let mut block = Block::<true, true>::new(ref_used.len(), query_used.len(), max_size);
        block.align_profile(&r_padded, &profile, min_size..=max_size, self.xdrop);
        let res = block.res();

        // Retrieve CIGAR
        let mut cigar = BlockCigar::new(res.query_idx, res.reference_idx);
        block.trace().cigar_eq(
            &r_padded,
            &q_padded,
            res.query_idx,
            res.reference_idx,
            &mut cigar,
        );

        if reverse {
            // Return the cigar reversed with insertions and deletions swapped
            AlignmentResult {
                score: res.score,
                query_start: query.len() - res.reference_idx,
                query_end: query.len(),
                ref_start: reference.len() - res.query_idx,
                ref_end: reference.len(),
                cigar: build_cigar_reverse_swap_id(&cigar),
            }
        } else {
            // Return the cigar with insertions and deletions swapped
            AlignmentResult {
                score: res.score,
                query_start: 0,
                query_end: res.reference_idx,
                ref_start: 0,
                ref_end: res.query_idx,
                cigar: build_cigar_swap_id(&cigar),
            }
        }
    }

    fn align_before_first_anchor(
        &self,
        query: &[u8],
        refseq: &[u8],
        first_anchor: &Anchor,
        padding: usize,
        result: &mut AlignmentInfo,
    ) {
        if first_anchor.query_start > 0 && first_anchor.ref_start > 0 {
            let query_part = &query[..first_anchor.query_start];
            let ref_start = 0
                .max(first_anchor.ref_start as isize - query_part.len() as isize - padding as isize)
                as usize;
            let ref_part = &refseq[ref_start..first_anchor.ref_start];

            let pre_align = self.xdrop_alignment(query_part, ref_part, true);

            if pre_align.score == 0 {
                result.query_start = first_anchor.query_start;
                result.ref_start = first_anchor.ref_start;
                result
                    .cigar
                    .push(CigarOperation::Softclip, result.query_start);
            } else {
                result.score += pre_align.score as u32;
                result.query_start = pre_align.query_start;
                result.ref_start = ref_start + pre_align.ref_start;
                if result.query_start > 0 {
                    result
                        .cigar
                        .push(CigarOperation::Softclip, result.query_start);
                }
                result.cigar.extend(&pre_align.cigar);
            }
        } else {
            result.query_start = first_anchor.query_start;
            result.ref_start = first_anchor.ref_start;
            if result.query_start == 0 {
                result.score += self.scores.end_bonus;
            } else {
                result
                    .cigar
                    .push(CigarOperation::Softclip, result.query_start);
            }
        }
    }

    fn align_after_last_anchor(
        &self,
        query: &[u8],
        refseq: &[u8],
        last_anchor: &Anchor,
        padding: usize,
        result: &mut AlignmentInfo,
    ) {
        let last_anchor_query_end = last_anchor.query_start + self.k;
        let last_anchor_ref_end = last_anchor.ref_start + self.k;

        if last_anchor_query_end < query.len() && last_anchor_ref_end < refseq.len() {
            let query_part = &query[last_anchor_query_end..];
            let ref_end = refseq
                .len()
                .min(last_anchor_ref_end + query_part.len() + padding);
            let ref_part = &refseq[last_anchor_ref_end..ref_end];

            let post_align = self.xdrop_alignment(query_part, ref_part, false);

            if post_align.score == 0 {
                result.query_end = last_anchor_query_end;
                result.ref_end = last_anchor_ref_end;
                result.cigar.push(
                    CigarOperation::Softclip,
                    query.len() - last_anchor_query_end,
                );
            } else {
                result.score += post_align.score as u32;
                result.query_end = last_anchor_query_end + post_align.query_end;
                result.ref_end = last_anchor_ref_end + post_align.ref_end;
                if result.query_end < query.len() {
                    result
                        .cigar
                        .push(CigarOperation::Softclip, query.len() - result.query_end);
                }
                result.cigar.extend(&post_align.cigar);
            }
        } else {
            result.query_end = last_anchor_query_end;
            result.ref_end = last_anchor_ref_end;
            if result.query_end == query.len() {
                result.score = self.scores.end_bonus;
            } else {
                result
                    .cigar
                    .push(CigarOperation::Softclip, query.len() - result.query_end);
            }
        }
    }

    pub fn piecewise_extension(
        &self,
        query: &[u8],
        refseq: &[u8],
        anchors: &[Anchor],
        padding: usize,
    ) -> AlignmentInfo {
        let mut result = AlignmentInfo {
            cigar: Cigar::new(),
            edit_distance: 0,
            ref_start: 0,
            ref_end: 0,
            query_start: 0,
            query_end: 0,
            score: 0,
        };

        // Maybe use the original vector instead of copying it ?
        let mut anchors: Vec<Anchor> = anchors.to_vec();
        remove_spurious_anchors(&mut anchors);

        self.align_before_first_anchor(query, refseq, &anchors[0], padding, &mut result);

        result.score += self.k as u32 * self.scores.match_ as u32;
        result.cigar.push(CigarOperation::Eq, self.k);

        for i in 1..anchors.len() {
            let anchor = &anchors[i];
            let prev_anchor = &anchors[i - 1];

            let curr_start_query = anchor.query_start;
            let curr_start_ref = anchor.ref_start;
            let prev_end_query = prev_anchor.query_start + self.k;
            let prev_end_ref = prev_anchor.ref_start + self.k;

            let query_diff = curr_start_query as isize - prev_end_query as isize;
            let ref_diff = curr_start_ref as isize - prev_end_ref as isize;

            // Early return for piecewise global alignment heuristic:
            // Because of spurious strobemer matches on the end of the query, chains tend to have a cluster of spurious
            // anchors on the ends of the query, this heuristic stop trusting the anchors if they create a suspiciously
            // large deletion on the end of the query
            if query.len() - curr_start_query <= 200
                && ref_diff - query_diff
                    >= query.len() as isize - prev_end_query as isize / 2 - self.k as isize
            {
                self.align_after_last_anchor(
                    query,
                    refseq,
                    prev_anchor,
                    (padding as isize + ref_diff - query_diff) as usize,
                    &mut result,
                );
                result.edit_distance = result.cigar.edit_distance();
                return result;
            }

            if ref_diff > 0 && query_diff > 0 {
                let query_part = &query[prev_end_query..query_diff as usize];
                let ref_part = &refseq[prev_end_ref..ref_diff as usize];

                println!("{:?} and {:?}", query_part, ref_part);

                if ref_diff == query_diff {
                    let hamming_aligned = hamming_align_global(
                        query_part,
                        ref_part,
                        self.scores.match_,
                        self.scores.mismatch,
                    );

                    if hamming_aligned.score
                        >= self.scores.match_ as u32
                            * ((query_part.len() as f32) * 0.85).ceil() as u32
                    {
                        result.score += hamming_aligned.score;
                        result.cigar.extend(&hamming_aligned.cigar);

                        result.score += self.k as u32 * self.scores.match_ as u32;
                        result.cigar.push(CigarOperation::Eq, self.k);
                    }
                }

                let aligned = self.global_alignment(query_part, ref_part);

                result.score = (aligned.score + result.score as i32) as u32;
                result.cigar.extend(&aligned.cigar);

                result.score += self.k as u32 * self.scores.match_ as u32;
                result.cigar.push(CigarOperation::Eq, self.k);
            } else {
                {
                    // Overlap between anchors, no need to align
                    if ref_diff < query_diff {
                        let inserted_part = (-ref_diff + query_diff) as usize;
                        result.score += (-(self.scores.gap_open as isize)
                            + (inserted_part as isize - 1) * -(self.scores.gap_extend as isize))
                            as u32;
                        result.cigar.push(CigarOperation::Insertion, inserted_part);

                        let matching_part = (self.k as isize + ref_diff) as usize;
                        result.score += (matching_part * self.scores.match_ as usize) as u32;
                        result.cigar.push(CigarOperation::Eq, matching_part);
                    } else if ref_diff > query_diff {
                        let deleted_part = (-query_diff + ref_diff) as usize;
                        result.score += (-(self.scores.gap_open as isize)
                            + (deleted_part as isize - 1) * -(self.scores.gap_extend as isize))
                            as u32;
                        result.cigar.push(CigarOperation::Deletion, deleted_part);

                        let matching_part = (self.k as isize + query_diff) as usize;
                        result.score += (matching_part * self.scores.match_ as usize) as u32;
                        result.cigar.push(CigarOperation::Eq, matching_part);
                    } else {
                        let matching_part = (self.k as isize + query_diff) as usize;
                        result.score += (matching_part * self.scores.match_ as usize) as u32;
                        result.cigar.push(CigarOperation::Eq, matching_part);
                    }
                }
            }
        }

        self.align_after_last_anchor(query, refseq, anchors.last().unwrap(), padding, &mut result);

        result.edit_distance = result.cigar.edit_distance();
        result
    }
}

// Converts a blockaligner Cigar to our Cigar format
fn build_cigar(block_cigar: &BlockCigar) -> Cigar {
    let mut result = Cigar::new();

    for i in 0..block_cigar.len() {
        let oplen = block_cigar.get(i);
        let op_code = match oplen.op {
            Operation::M => CigarOperation::Match,
            Operation::Eq => CigarOperation::Eq,
            Operation::X => CigarOperation::X,
            Operation::I => CigarOperation::Insertion,
            Operation::D => CigarOperation::Deletion,
            _ => continue,
        };
        result.push(op_code, oplen.len);
    }

    result
}

// Converts a blockaligner Cigar to our Cigar format with I/D swapped
fn build_cigar_swap_id(block_cigar: &BlockCigar) -> Cigar {
    let mut result = Cigar::new();

    for i in 0..block_cigar.len() {
        let oplen = block_cigar.get(i);
        let op_code = match oplen.op {
            Operation::M => CigarOperation::Match,
            Operation::Eq => CigarOperation::Eq,
            Operation::X => CigarOperation::X,
            Operation::I => CigarOperation::Deletion, // Swapped
            Operation::D => CigarOperation::Insertion, // Swapped
            _ => continue,
        };
        result.push(op_code, oplen.len);
    }

    result
}

// Converts a blockaligner Cigar to our Cigar format in reverse order with I/D swapped
fn build_cigar_reverse_swap_id(block_cigar: &BlockCigar) -> Cigar {
    let mut result = Cigar::new();

    for i in (0..block_cigar.len()).rev() {
        let oplen = block_cigar.get(i);
        let op_code = match oplen.op {
            Operation::M => CigarOperation::Match,
            Operation::Eq => CigarOperation::Eq,
            Operation::X => CigarOperation::X,
            Operation::I => CigarOperation::Deletion, // Swapped
            Operation::D => CigarOperation::Insertion, // Swapped
            _ => continue,
        };
        result.push(op_code, oplen.len);
    }

    result
}

fn remove_spurious_anchors(anchors: &mut Vec<Anchor>) {
    if anchors.len() < 2 {
        return;
    }

    // First pruning:
    // We remove all cluster of anchors creating 2 canceling indels within the tolerance range.
    // This will only happen inside the chain, so even if they were part of the best scoring path
    // they will be retrieved in the global alignment.

    let diag_tolerance = 5;

    let mut tracked_indel = 0;
    let mut deviation_start = None;

    let mut i = 1;
    while i < anchors.len() {
        let query_diff = (anchors[i].query_start as isize) - (anchors[i - 1].query_start as isize);
        let ref_diff = (anchors[i].ref_start as isize) - (anchors[i - 1].ref_start as isize);

        let indel = query_diff - ref_diff;

        // 2 anchors are consider to deviate from the diagonal if
        // they create an indel bigger than the diagonal tolerance
        if indel.abs() > diag_tolerance {
            if let Some(index) = deviation_start {
                if (tracked_indel + indel).abs() <= diag_tolerance {
                    // if we find an indel canceling the tracked indel, we
                    // remove all anchors covered between the 2 indel positions
                    anchors.drain(index..i);
                    i -= i - index;
                    deviation_start = None;
                    tracked_indel = 0;
                }
            } else {
                deviation_start = Some(i);
                tracked_indel = indel;
            }
        }
        i += 1;
    }

    // Second pruning:
    // We remove anchors of the ends of the chain if they create any indel of min_indel_size.
    // If they were part of the best scoring path, they should be retireved by the x-drop alignment.
    // the idea is that spurious anchors are much more difficult to detect on the ends of the chain,
    // so we look at a small ratio of anchors and remove any anchors creating deviations.

    let edge_prune_ratio = 0.1;
    let max_prune_count = ((anchors.len() as f32) * edge_prune_ratio).ceil() as usize;

    for i in 1..anchors.len().max(max_prune_count) {
        let query_diff = (anchors[i].query_start as isize) - (anchors[i - 1].query_start as isize);
        let ref_diff = (anchors[i].ref_start as isize) - (anchors[i - 1].ref_start as isize);

        let indel = query_diff - ref_diff;

        if indel.abs() != 0 {
            anchors.drain(..i);
            break;
        }
    }

    let max_prune_count = ((anchors.len() as f32) * (1.0 - edge_prune_ratio)).ceil() as usize;

    for i in (max_prune_count..anchors.len()).rev() {
        let query_diff = (anchors[i].query_start as isize) - (anchors[i - 1].query_start as isize);
        let ref_diff = (anchors[i].ref_start as isize) - (anchors[i - 1].ref_start as isize);

        let indel = query_diff - ref_diff;

        if indel.abs() != 0 {
            anchors.drain(i..);
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_global_perfect_match() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            0,
        );
        let result = aligner.global_alignment(b"AAATTT", b"AAATTT");
        assert_eq!(result.score, 6 * Scores::default().match_ as i32);
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 6);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 6);
        assert_eq!(result.cigar.to_string(), "6=");
    }

    #[test]
    fn test_global_complete_mismatch() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            0,
        );
        let result = aligner.global_alignment(b"AAA", b"TTT");
        assert_eq!(result.score, 3 * -(Scores::default().mismatch as i32));
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 3);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 3);
        assert_eq!(result.cigar.to_string(), "3X");
    }

    #[test]
    fn test_global_single_mismatch() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            0,
        );
        let result = aligner.global_alignment(b"AAATAA", b"AAAAAA");
        assert_eq!(
            result.score,
            5 * Scores::default().match_ as i32 - Scores::default().mismatch as i32
        );
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 6);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 6);
        assert_eq!(result.cigar.to_string(), "3=1X2=");
    }

    #[test]
    fn test_global_gap_in_query() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            0,
        );
        let result = aligner.global_alignment(b"AAATTT", b"AAAATTTT");
        assert_eq!(
            result.score,
            6 * Scores::default().match_ as i32
                - Scores::default().gap_open as i32
                - Scores::default().gap_extend as i32
        );
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 6);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 8);
        assert_eq!(result.cigar.to_string(), "3=2D3=");
    }

    #[test]
    fn test_global_gap_in_reference() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            0,
        );
        let result = aligner.global_alignment(b"AAAATTTT", b"AAATTT");
        assert_eq!(
            result.score,
            6 * Scores::default().match_ as i32
                - Scores::default().gap_open as i32
                - Scores::default().gap_extend as i32,
        );
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 8);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 6);
        assert_eq!(result.cigar.to_string(), "3=2I3=");
    }

    #[test]
    fn test_global_gap_at_query_start() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            0,
        );
        let result = aligner.global_alignment(b"AAATTT", b"TTAAATTT");
        assert_eq!(
            result.score,
            6 * Scores::default().match_ as i32
                - Scores::default().gap_open as i32
                - Scores::default().gap_extend as i32
        );
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 6);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 8);
        assert_eq!(result.cigar.to_string(), "2D6=");
    }

    #[test]
    fn test_global_gap_at_query_end() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            0,
        );
        let result = aligner.global_alignment(b"AAATTT", b"AAATTTAA");

        assert_eq!(
            result.score,
            6 * Scores::default().match_ as i32
                - Scores::default().gap_open as i32
                - Scores::default().gap_extend as i32
        );
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 6);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 8);
        assert_eq!(result.cigar.to_string(), "6=2D");
    }

    #[test]
    fn test_global_gap_at_reference_start() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            0,
        );
        let result = aligner.global_alignment(b"TTAAATTT", b"AAATTT");

        assert_eq!(
            result.score,
            6 * Scores::default().match_ as i32
                - Scores::default().gap_open as i32
                - Scores::default().gap_extend as i32
        );
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 8);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 6);
        assert_eq!(result.cigar.to_string(), "2I6=");
    }

    #[test]
    fn test_global_gap_at_reference_end() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            0,
        );
        let result = aligner.global_alignment(b"AAATTTAA", b"AAATTT");

        assert_eq!(
            result.score,
            6 * Scores::default().match_ as i32
                - Scores::default().gap_open as i32
                - Scores::default().gap_extend as i32
        );
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 8);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 6);
        assert_eq!(result.cigar.to_string(), "6=2I");
    }

    #[test]
    fn test_global_multiple_mismatches() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            0,
        );
        let result = aligner.global_alignment(b"ATATATATAT", b"AAAAAAAAAA");

        assert_eq!(
            result.score,
            5 * Scores::default().match_ as i32 - 5 * Scores::default().mismatch as i32
        );
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 10);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 10);
        assert_eq!(result.cigar.to_string(), "1=1X1=1X1=1X1=1X1=1X");
    }

    #[test]
    fn test_global_complex_alignment() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            0,
        );
        let result = aligner.global_alignment(b"AAACTTAAACCTT", b"AAAATTGAAATT");
        assert_eq!(
            result.score,
            10 * Scores::default().match_ as i32
                - Scores::default().mismatch as i32
                - 2 * Scores::default().gap_open as i32
                - Scores::default().gap_extend as i32
        );
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 13);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 12);
        assert_eq!(result.cigar.to_string(), "3=1X2=1D3=2I2=");
    }

    #[test]
    fn test_xdrop_forward_perfect_match() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            100,
        );
        let result = aligner.xdrop_alignment(b"AAATTT", b"AAATTT", false);
        assert_eq!(result.score, 6 * 2 + 10);
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 6);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 6);
        assert_eq!(result.cigar.to_string(), "6=");
    }

    #[test]
    fn test_xdrop_forward_with_mismatch() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            100,
        );
        let result = aligner.xdrop_alignment(b"AAATAA", b"AAAAAA", false);
        assert_eq!(result.score, 5 * 2 - 8 + 10);
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 6);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 6);
        assert_eq!(result.cigar.to_string(), "3=1X2=");
    }

    #[test]
    fn test_xdrop_forward_early_termination() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            100,
        );
        let result = aligner.xdrop_alignment(b"AAAAAAATTTTTTT", b"AAAAAAACCCCCCC", false);
        assert_eq!(result.score, 7 * 2);
        assert_eq!(result.query_start, 0);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.query_end, 7);
        assert_eq!(result.ref_end, 7);
        assert_eq!(result.cigar.to_string(), "7=");
    }

    #[test]
    fn test_xdrop_forward_end_bonus_extends_alignment() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 50,
            },
            0,
            100,
        );
        let result = aligner.xdrop_alignment(b"AAAAAAAATT", b"AAAAAAAAAA", false);
        assert_eq!(result.score, 8 * 2 - 8 * 2 + 50);
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 10);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 10);
        assert_eq!(result.cigar.to_string(), "8=2X");
    }

    #[test]
    fn test_xdrop_forward_with_gap() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            100,
        );
        let result = aligner.xdrop_alignment(b"AAAAACCTTT", b"AAAAATTT", false);
        assert_eq!(result.score, 8 * 2 - 12 - 1 + 10);
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 10);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 8);
        assert_eq!(result.cigar.to_string(), "5=2I3=");
    }

    #[test]
    fn test_xdrop_forward_gap_in_reference() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            100,
        );
        let result = aligner.xdrop_alignment(b"AAAAATTT", b"AAAAACCTTT", false);
        assert_eq!(result.score, 8 * 2 - 12 - 1 + 10);
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 8);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 10);
        assert_eq!(result.cigar.to_string(), "5=2D3=");
    }

    #[test]
    fn test_xdrop_reverse_perfect_match() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            100,
        );
        let result = aligner.xdrop_alignment(b"AAATTT", b"AAATTT", true);
        assert_eq!(result.score, 6 * 2 + 10);
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 6);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 6);
        assert_eq!(result.cigar.to_string(), "6=");
    }

    #[test]
    fn test_xdrop_reverse_with_mismatch() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            100,
        );
        let result = aligner.xdrop_alignment(b"AATAAA", b"AAAAAA", true);
        assert_eq!(result.score, 5 * 2 - 8 + 10);
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 6);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 6);
        assert_eq!(result.cigar.to_string(), "2=1X3=");
    }

    #[test]
    fn test_xdrop_reverse_early_termination() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            100,
        );
        let result = aligner.xdrop_alignment(b"TTTTTTTAAAAAAA", b"CCCCCCCAAAAAAA", true);
        assert_eq!(result.score, 7 * 2);
        assert_eq!(result.query_start, 7);
        assert_eq!(result.query_end, 14);
        assert_eq!(result.ref_start, 7);
        assert_eq!(result.ref_end, 14);
        assert_eq!(result.cigar.to_string(), "7=");
    }

    #[test]
    fn test_xdrop_reverse_end_bonus_extends_alignment() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 50,
            },
            0,
            100,
        );
        let result = aligner.xdrop_alignment(b"TTAAAAAAAA", b"AAAAAAAAAA", true);
        assert_eq!(result.score, 8 * 2 - 8 * 2 + 50);
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 10);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 10);
        assert_eq!(result.cigar.to_string(), "2X8=");
    }

    #[test]
    fn test_xdrop_reverse_with_gap() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            100,
        );
        let result = aligner.xdrop_alignment(b"TTTCCAAAAA", b"TTTAAAAA", true);
        assert_eq!(result.score, 8 * 2 - 12 - 1 + 10);
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 10);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 8);
        assert_eq!(result.cigar.to_string(), "3=2I5=");
    }

    #[test]
    fn test_xdrop_reverse_gap_in_reference() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            0,
            100,
        );
        let result = aligner.xdrop_alignment(b"TTTAAAAA", b"TTTCCAAAAA", true);
        assert_eq!(result.score, 8 * 2 - 12 - 1 + 10);
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 8);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 10);
        assert_eq!(result.cigar.to_string(), "3=2D5=");
    }
}
