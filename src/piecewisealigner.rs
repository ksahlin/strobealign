use crate::{
    aligner::{AlignmentInfo, Scores, hamming_align_global},
    chainer::Anchor,
    cigar::{Cigar, CigarOperation},
};
use block_aligner::{
    cigar::{Cigar as BlockCigar, Operation},
    scan_block::{Block, PaddedBytes},
    scores::{AAMatrix, AAProfile, Gaps, NucMatrix, Profile},
};

// Maximum value for blockaligner's block sizes, higher values might cause a crash
const MAXIMUM_BLOCK_SIZE: usize = 8192;

#[derive(Clone, Copy)]
enum XdropMode {
    GlobalStartLocalEnd,
    LocalStartGlobalEnd,
}

#[derive(Clone)]
pub struct PiecewiseAligner {
    scores: Scores,
    k: usize,
    gaps: Gaps,
    matrix: NucMatrix,
    xdrop: i32,
}

struct AlignmentResult {
    pub score: i32,
    pub query_start: usize,
    pub query_end: usize,
    pub ref_start: usize,
    pub ref_end: usize,
    pub cigar: Cigar,
}

impl Default for AlignmentResult {
    fn default() -> Self {
        AlignmentResult {
            score: 0,
            query_start: 0,
            query_end: 0,
            ref_start: 0,
            ref_end: 0,
            cigar: Cigar::new(),
        }
    }
}

impl PiecewiseAligner {
    pub fn new(scores: Scores, k: usize, xdrop: i32) -> Self {
        let gaps = Gaps {
            open: -(scores.gap_open as i8),
            extend: -(scores.gap_extend as i8),
        };
        let matrix = NucMatrix::new_simple(scores.match_ as i8, -(scores.mismatch as i8));
        PiecewiseAligner {
            scores,
            k,
            gaps,
            matrix,
            xdrop,
        }
    }

    /// Performs a global alignment between two sequences.
    ///
    /// The alignment uses the blockaligner crate with block sizes
    /// determined based on the largest length of the 2 sequences.
    ///
    /// # Parameters
    ///
    /// * `query` - The query sequence as a byte slice
    /// * `reference` - The reference sequence as a byte slice
    ///
    /// # Returns
    ///
    /// An [`AlignmentResult`] containing:
    /// - The alignment score
    /// - Start and end positions in both sequences
    /// - A CIGAR string representing the alignment operations
    fn global_alignment(&self, query: &[u8], reference: &[u8]) -> AlignmentResult {
        if query.is_empty() || reference.is_empty() {
            return AlignmentResult::default();
        }

        // Blockaligner block ranges, set to maximum values for full global alignment
        // Must be powers of 2
        let max_len = query.len().max(reference.len());
        let block_size = max_len.next_power_of_two().clamp(32, MAXIMUM_BLOCK_SIZE);

        // Create padded sequences
        let mut q_padded = PaddedBytes::new::<NucMatrix>(query.len(), block_size);
        q_padded.set_bytes::<NucMatrix>(query, block_size);
        let mut r_padded = PaddedBytes::new::<NucMatrix>(reference.len(), block_size);
        r_padded.set_bytes::<NucMatrix>(reference, block_size);

        // Make a global alignment call with traceback
        let mut block = Block::<true, false>::new(query.len(), reference.len(), block_size);
        block.align(
            &q_padded,
            &r_padded,
            &self.matrix,
            self.gaps,
            block_size..=block_size,
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

    /// Performs an alignment with one end fixed and the other allowed to terminate
    /// early based on the x-drop mode.
    ///
    /// This method computes a semi-global alignment where one end of the alignment is
    /// anchored while the other end can terminate early if the alignment score drops
    /// too far below the maximum observed score (determined by the x-drop value).
    ///
    /// The alignment uses the blockaligner crate with block sizes chosen dynamically
    /// as 20% of the maximum length of the 2 sequences and 50% of the maximum length
    /// of the 2 sequences.
    ///
    /// # Parameters
    ///
    /// * `query` - The query sequence as a byte slice
    /// * `reference` - The reference sequence as a byte slice
    /// * `mode` - Controls which end is fixed:
    ///   - [`XdropMode::GlobalStartLocalEnd`]: Alignment must start at position 0 of
    ///     both sequences but may end before the last position
    ///   - [`XdropMode::LocalStartGlobalEnd`]: Alignment must end at the last position
    ///     of both sequences but may start after position 0. This is implemented by
    ///     reversing both sequences, aligning them, and transforming the result back
    ///     to forward coordinates.
    ///
    /// # Returns
    ///
    /// An [`AlignmentResult`] containing:
    /// - The alignment score
    /// - Start and end positions in both sequences
    /// - A CIGAR string representing the alignment operations
    fn xdrop_alignment(&self, query: &[u8], reference: &[u8], mode: XdropMode) -> AlignmentResult {
        if query.is_empty() || reference.is_empty() {
            return AlignmentResult::default();
        }

        // Blockaligner block ranges, set to 20% and 50% of longest sequence
        // Must be powers of 2
        let max_len = query.len().max(reference.len());
        let min_size = (max_len / 5)
            .next_power_of_two()
            .clamp(32, MAXIMUM_BLOCK_SIZE);
        let max_size = (max_len / 2)
            .next_power_of_two()
            .clamp(128, MAXIMUM_BLOCK_SIZE);

        // Create padded sequences
        // blockaligner only has a AA profile, no nucleotide profile, so we have to use AA matrix
        let mut q_padded = PaddedBytes::new::<AAMatrix>(query.len(), max_size);
        let mut r_padded = PaddedBytes::new::<AAMatrix>(reference.len(), max_size);

        match mode {
            XdropMode::GlobalStartLocalEnd => {
                q_padded.set_bytes::<AAMatrix>(query, max_size);
                r_padded.set_bytes::<AAMatrix>(reference, max_size);
            }
            XdropMode::LocalStartGlobalEnd => {
                q_padded.set_bytes_rev::<AAMatrix>(query, max_size);
                r_padded.set_bytes_rev::<AAMatrix>(reference, max_size);
            }
        }

        // Create profile to align the reference on the query with end bonus
        // We need to build a position-specific scoring matrix
        // blockaligner only has a AA profile, no nucleotide profile
        let profile = make_aa_profile(query, &self.scores, max_size, mode);

        // Make an x-drop alignment call with traceback
        let mut block = Block::<true, true>::new(reference.len(), query.len(), max_size);
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
        match mode {
            XdropMode::GlobalStartLocalEnd => AlignmentResult {
                score: res.score,
                query_start: 0,
                query_end: res.reference_idx,
                ref_start: 0,
                ref_end: res.query_idx,
                cigar: build_cigar_swap_indel(&cigar),
            },
            XdropMode::LocalStartGlobalEnd => AlignmentResult {
                score: res.score,
                query_start: query.len() - res.reference_idx,
                query_end: query.len(),
                ref_start: reference.len() - res.query_idx,
                ref_end: reference.len(),
                cigar: build_cigar_reverse_swap_indel(&cigar),
            },
        }
    }

    /// Aligns the query sequence region before the first anchor.
    ///
    /// This method handles the prefix alignment by taking the query sequence up to the
    /// first anchor's start position and aligning it against a corresponding region in
    /// the reference sequence. The reference region is extended backward from the anchor
    /// by the query prefix length plus padding to allow for potential indels.
    ///
    /// # Parameters
    ///
    /// * `query` - The full query sequence as a byte slice
    /// * `refseq` - The full reference sequence as a byte slice
    /// * `first_anchor` - The first anchor
    /// * `padding` - Additional bases to include in the reference region to account for
    ///   potential insertions/deletions
    ///
    /// # Returns
    ///
    /// An [`AlignmentResult`] containing:
    /// - The alignment score (including end bonus if alignment starts at the beginning of
    ///   the query)
    /// - Start and end positions in both sequences
    /// - A CIGAR string representing the alignment operations
    ///
    /// If there is no sequence before the first anchor or the alignment score is 0,
    /// returns a result indicating the alignment starts at the first anchor position.
    fn align_before_first_anchor(
        &self,
        query: &[u8],
        refseq: &[u8],
        first_anchor: &Anchor,
        padding: usize,
    ) -> AlignmentResult {
        if first_anchor.query_start > 0 && first_anchor.ref_start > 0 {
            let query_part = &query[..first_anchor.query_start];
            let ref_start = first_anchor
                .ref_start
                .saturating_sub(query_part.len() + padding);
            let ref_part = &refseq[ref_start..first_anchor.ref_start];

            let mut pre_align =
                self.xdrop_alignment(query_part, ref_part, XdropMode::LocalStartGlobalEnd);

            if pre_align.score == 0 {
                AlignmentResult {
                    query_start: first_anchor.query_start,
                    ref_start: first_anchor.ref_start,
                    ..Default::default()
                }
            } else {
                pre_align.ref_start += ref_start;
                pre_align
            }
        } else {
            AlignmentResult {
                query_start: first_anchor.query_start,
                ref_start: first_anchor.ref_start,
                score: if first_anchor.query_start == 0 {
                    self.scores.end_bonus as i32
                } else {
                    0
                },
                ..Default::default()
            }
        }
    }

    /// Aligns the query sequence region after the last anchor.
    ///
    /// This method handles the suffix alignment by taking the query sequence after the
    /// last anchor's end position and aligning it against a corresponding region in
    /// the reference sequence. The reference region is extended forward from the anchor
    /// by the query suffix length plus padding to allow for potential indels.
    ///
    /// # Parameters
    ///
    /// * `query` - The full query sequence as a byte slice
    /// * `refseq` - The full reference sequence as a byte slice
    /// * `last_anchor` - The last anchor
    /// * `padding` - Additional bases to include in the reference region to account for
    ///   potential insertions/deletions
    ///
    /// # Returns
    ///
    /// An [`AlignmentResult`] containing:
    /// - The alignment score (including end bonus if alignment reaches the end of the query)
    /// - Start and end positions in both sequences
    /// - A CIGAR string representing the alignment operations
    ///
    /// If there is no sequence after the last anchor or the alignment score is 0,
    /// returns a result indicating the alignment ends at the last anchor position.
    fn align_after_last_anchor(
        &self,
        query: &[u8],
        refseq: &[u8],
        last_anchor: &Anchor,
        padding: usize,
    ) -> AlignmentResult {
        let last_anchor_query_end = last_anchor.query_start + self.k;
        let last_anchor_ref_end = last_anchor.ref_start + self.k;

        if last_anchor_query_end < query.len() && last_anchor_ref_end < refseq.len() {
            let query_part = &query[last_anchor_query_end..];
            let ref_end = refseq
                .len()
                .min(last_anchor_ref_end + query_part.len() + padding);
            let ref_part = &refseq[last_anchor_ref_end..ref_end];

            let mut post_align =
                self.xdrop_alignment(query_part, ref_part, XdropMode::GlobalStartLocalEnd);

            if post_align.score == 0 {
                AlignmentResult {
                    query_end: last_anchor_query_end,
                    ref_end: last_anchor_ref_end,
                    ..Default::default()
                }
            } else {
                post_align.query_end += last_anchor_query_end;
                post_align.ref_end += last_anchor_ref_end;
                post_align
            }
        } else {
            AlignmentResult {
                query_end: last_anchor_query_end,
                ref_end: last_anchor_ref_end,
                score: if last_anchor_query_end == query.len() {
                    self.scores.end_bonus as i32
                } else {
                    0
                },
                ..Default::default()
            }
        }
    }

    /// Performs piecewise alignment between query and reference sequences using anchor.
    ///
    /// This method implements a piecewise alignment extension strategy that uses previously
    /// calculated anchors from chaining (exact k-mer matches) to guide the alignment.
    /// The alignment is constructed by:
    /// -  Aligning the sequence before the first anchor
    /// -  Adding each anchor as a perfect match of length k
    /// -  Aligning regions between consecutive anchors
    /// -  Aligning the sequence after the last anchor
    ///
    /// Between anchors, the method uses different strategies based on the distance:
    /// - If both sequences have positive gaps of equal length, attempts Hamming distance
    ///   alignment first (faster), falling back to full global alignment if quality is poor
    /// - If gaps are unequal or negative (overlapping anchors), handles insertions/deletions
    ///   directly without alignment
    ///
    /// # Parameters
    ///
    /// * `query` - The query sequence as a byte slice
    /// * `refseq` - The reference sequence as a byte slice
    /// * `anchors` - A slice of anchor sorted from last to first (reverse order)
    /// * `padding` - Additional bases to include when aligning prefix/suffix regions to
    ///   account for potential insertions/deletions
    ///
    /// # Returns
    ///
    /// An [`Option<AlignmentInfo>`] containing:
    /// - `Some(AlignmentInfo)` with the alignment score, CIGAR string, edit distance,
    ///   and start/end positions in both sequences if the alignment score is positive
    /// - `None` if the alignment score is <= 0
    pub fn extend_piecewise(
        &self,
        query: &[u8],
        refseq: &[u8],
        anchors: &[Anchor],
        padding: usize,
    ) -> Option<AlignmentInfo> {
        let AlignmentResult {
            mut cigar,
            ref_start,
            query_start,
            mut score,
            ..
        } = self.align_before_first_anchor(query, refseq, anchors.last()?, padding);

        score += self.k as i32 * self.scores.match_ as i32;
        cigar.push(CigarOperation::Eq, self.k);

        for i in (1..anchors.len()).rev() {
            let curr_query_start = anchors[i - 1].query_start;
            let curr_ref_start = anchors[i - 1].ref_start;
            let prev_query_end = anchors[i].query_start + self.k;
            let prev_ref_end = anchors[i].ref_start + self.k;

            let query_diff = curr_query_start as isize - prev_query_end as isize;
            let ref_diff = curr_ref_start as isize - prev_ref_end as isize;

            if ref_diff > 0 && query_diff > 0 {
                let query_part = &query[prev_query_end..prev_query_end + query_diff as usize];
                let ref_part = &refseq[prev_ref_end..prev_ref_end + ref_diff as usize];

                if ref_diff == query_diff {
                    let hamming_aligned = hamming_align_global(
                        query_part,
                        ref_part,
                        self.scores.match_,
                        self.scores.mismatch,
                    );

                    if hamming_aligned.score > 0
                        && hamming_aligned.score
                            > self.scores.match_ as u32
                                * ((query_part.len() as f32) * 0.85).ceil() as u32
                    {
                        score += hamming_aligned.score as i32;
                        cigar.extend(&hamming_aligned.cigar);

                        score += self.k as i32 * self.scores.match_ as i32;
                        cigar.push(CigarOperation::Eq, self.k);
                        continue;
                    }
                }

                let aligned = self.global_alignment(query_part, ref_part);

                score += aligned.score;
                cigar.extend(&aligned.cigar);

                score += self.k as i32 * self.scores.match_ as i32;
                cigar.push(CigarOperation::Eq, self.k);
            } else {
                // Overlap between anchors, no need to align
                if ref_diff < query_diff {
                    let inserted_part = (query_diff - ref_diff) as usize;

                    score += -(self.scores.gap_open as i32)
                        + (inserted_part as i32 - 1) * -(self.scores.gap_extend as i32);
                    cigar.push(CigarOperation::Insertion, inserted_part);

                    let matching_part = (self.k as isize + ref_diff) as usize;
                    score += matching_part as i32 * self.scores.match_ as i32;
                    cigar.push(CigarOperation::Eq, matching_part);
                } else if ref_diff > query_diff {
                    let deleted_part = (ref_diff - query_diff) as usize;
                    score += -(self.scores.gap_open as i32)
                        + (deleted_part as i32 - 1) * -(self.scores.gap_extend as i32);
                    cigar.push(CigarOperation::Deletion, deleted_part);

                    let matching_part = (self.k as isize + query_diff) as usize;
                    score += matching_part as i32 * self.scores.match_ as i32;
                    cigar.push(CigarOperation::Eq, matching_part);
                } else {
                    let matching_part = (self.k as isize + query_diff) as usize;
                    score += matching_part as i32 * self.scores.match_ as i32;
                    cigar.push(CigarOperation::Eq, matching_part);
                }
            }
        }

        let AlignmentResult {
            cigar: end_cigar,
            ref_end,
            query_end,
            score: end_score,
            ..
        } = self.align_after_last_anchor(query, refseq, anchors.first()?, padding);

        cigar.extend(&end_cigar);
        score += end_score;
        let edit_distance = cigar.edit_distance();

        if score > 0 {
            Some(AlignmentInfo {
                cigar,
                edit_distance,
                ref_start,
                ref_end,
                query_start,
                query_end,
                score: score as u32,
            })
        } else {
            None
        }
    }
}

/// Converts a blockaligner CIGAR string to our internal CIGAR format.
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

/// Converts a blockaligner CIGAR string to our internal CIGAR format with
/// insertions and deletions swapped.
fn build_cigar_swap_indel(block_cigar: &BlockCigar) -> Cigar {
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

/// Converts a blockaligner CIGAR string to our internal CIGAR format in reverse
/// order with insertions and deletions swapped.
fn build_cigar_reverse_swap_indel(block_cigar: &BlockCigar) -> Cigar {
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

/// Creates a position-specific scoring matrix (profile) for aligning the reference
/// against the query with end bonuses.
///
/// # Parameters
///
/// * `query` - The query sequence to build the profile from
/// * `scores` - Scoring parameters
/// * `max_size` - Maximum block size for alignment
/// * `mode` - Alignment mode determining sequence traversal direction
///
/// # Returns
///
/// An [`AAProfile`] configured with position-specific scores for each nucleotide
/// pairing and appropriate end bonuses based on the alignment mode.
fn make_aa_profile(query: &[u8], scores: &Scores, max_size: usize, mode: XdropMode) -> AAProfile {
    let mut profile = AAProfile::new(query.len(), max_size, -(scores.gap_extend as i8));

    // Set scores for each nucleotide combination
    for i in 1..=query.len() {
        // Select the query character at this position depending on the alignment mode.
        // `GlobalStartLocalEnd`: traverse the query from start to end (forward).
        // `LocalStartGlobalEnd`: traverse the query from end to start (reversed).
        let query_char = match mode {
            XdropMode::GlobalStartLocalEnd => query[i - 1],
            XdropMode::LocalStartGlobalEnd => query[query.len() - i],
        };

        for &c in b"ACGTN" {
            if c == query_char {
                profile.set(i, c, scores.match_ as i8);
            } else {
                profile.set(i, c, -(scores.mismatch as i8));
            };
        }
    }

    // Set gap costs
    profile.set_all_gap_open_C(-(scores.gap_open as i8) - -(scores.gap_extend as i8));
    profile.set_all_gap_close_C(0);
    profile.set_all_gap_open_R(-(scores.gap_open as i8) - -(scores.gap_extend as i8));

    // Give a bonus score to the end
    let bonus_pos = query.len();
    for &c in b"ACGTN" {
        let current_score = profile.get(bonus_pos, c);
        profile.set(bonus_pos, c, current_score + scores.end_bonus as i8);
    }

    profile
}

/// Removes spurious anchors from a chain of anchor points.
///
/// This method applies two pruning strategies to filter out unreliable anchors that
/// may interfere with the optimal alignment. Spurious anchors can arise from repetitive
/// sequences or random k-mer matches and can disrupt the alignment path.
///
/// # Pruning Strategies
///
/// ## First Pruning: Canceling Indels
/// Removes clusters of anchors that create two cancelling indels within the diagonal
/// tolerance (5 bases).
///
/// ## Second Pruning: Edge Anchors
/// Removes anchors at the beginning and end of the chain (up to 20% of the chain
/// length) if they create any indel.
///
/// # Parameters
///
/// * `anchors` - A mutable vector of anchor points to be pruned in place
pub fn remove_spurious_anchors(anchors: &mut Vec<Anchor>) {
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
        let query_diff = (anchors[i - 1].query_start as isize) - (anchors[i].query_start as isize);
        let ref_diff = (anchors[i - 1].ref_start as isize) - (anchors[i].ref_start as isize);

        let indel = query_diff - ref_diff;

        // 2 anchors are consider to deviate from the diagonal if
        // they create an indel bigger than the diagonal tolerance
        if indel.abs() > diag_tolerance {
            if let Some(index) = deviation_start
                && (tracked_indel + indel).abs() <= diag_tolerance
            {
                // if we find an indel canceling the tracked indel, we
                // remove all anchors covered between the 2 indel positions
                anchors.drain(index..i);
                i = index;
                deviation_start = None;
                tracked_indel = 0;
                continue;
            }
            deviation_start = Some(i);
            tracked_indel = indel;
        }
        i += 1;
    }

    // Second pruning:
    // We remove anchors of the ends of the chain if they create any indel.
    // If they were part of the best scoring path, they should be retireved by the x-drop alignment.
    // the idea is that spurious anchors are much more difficult to detect on the ends of the chain,
    // so we look at a small ratio of anchors and remove any anchors creating deviations.

    let edge_prune_ratio = 0.20;
    let max_prune_count = ((anchors.len() as f32) * edge_prune_ratio).ceil() as usize;

    for i in 1..anchors.len().min(max_prune_count) {
        let query_diff = (anchors[i].query_start as isize) - (anchors[i - 1].query_start as isize);
        let ref_diff = (anchors[i].ref_start as isize) - (anchors[i - 1].ref_start as isize);

        let indel = query_diff - ref_diff;

        if indel.abs() != 0 {
            anchors.drain(..i);
            break;
        }
    }

    for i in (anchors.len().saturating_sub(max_prune_count) + 1..anchors.len()).rev() {
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
        assert_eq!(result.score, 6 * 2);
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
        assert_eq!(result.score, 3 * -8);
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
        assert_eq!(result.score, 5 * 2 - 8);
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
        assert_eq!(result.score, 6 * 2 - 12 - 1);
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
        assert_eq!(result.score, 6 * 2 - 12 - 1,);
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
        assert_eq!(result.score, 6 * 2 - 12 - 1);
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

        assert_eq!(result.score, 6 * 2 - 12 - 1);
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

        assert_eq!(result.score, 6 * 2 - 12 - 1);
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

        assert_eq!(result.score, 6 * 2 - 12 - 1);
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

        assert_eq!(result.score, 5 * 2 - 5 * 8);
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
        assert_eq!(result.score, 10 * 2 - 8 - 2 * 12 - 1);
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
        let result = aligner.xdrop_alignment(b"AAATTT", b"AAATTT", XdropMode::GlobalStartLocalEnd);
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
        let result = aligner.xdrop_alignment(b"AAATAA", b"AAAAAA", XdropMode::GlobalStartLocalEnd);
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
        let result = aligner.xdrop_alignment(
            b"AAAAAAATTTTTTT",
            b"AAAAAAACCCCCCC",
            XdropMode::GlobalStartLocalEnd,
        );
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
        let result =
            aligner.xdrop_alignment(b"AAAAAAAATT", b"AAAAAAAAAA", XdropMode::GlobalStartLocalEnd);
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
        let result =
            aligner.xdrop_alignment(b"AAAAACCTTT", b"AAAAATTT", XdropMode::GlobalStartLocalEnd);
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
        let result =
            aligner.xdrop_alignment(b"AAAAATTT", b"AAAAACCTTT", XdropMode::GlobalStartLocalEnd);
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
        let result = aligner.xdrop_alignment(b"AAATTT", b"AAATTT", XdropMode::LocalStartGlobalEnd);
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
        let result = aligner.xdrop_alignment(b"AATAAA", b"AAAAAA", XdropMode::LocalStartGlobalEnd);
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
        let result = aligner.xdrop_alignment(
            b"TTTTTTTAAAAAAA",
            b"CCCCCCCAAAAAAA",
            XdropMode::LocalStartGlobalEnd,
        );
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
        let result =
            aligner.xdrop_alignment(b"TTAAAAAAAA", b"AAAAAAAAAA", XdropMode::LocalStartGlobalEnd);
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
        let result =
            aligner.xdrop_alignment(b"TTTCCAAAAA", b"TTTAAAAA", XdropMode::LocalStartGlobalEnd);
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
        let result =
            aligner.xdrop_alignment(b"TTTAAAAA", b"TTTCCAAAAA", XdropMode::LocalStartGlobalEnd);
        assert_eq!(result.score, 8 * 2 - 12 - 1 + 10);
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 8);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 10);
        assert_eq!(result.cigar.to_string(), "3=2D5=");
    }

    #[test]
    fn test_remove_spurious_anchors_control() {
        let expected = vec![
            Anchor {
                ref_id: 0,
                ref_start: 0,
                query_start: 0,
            },
            Anchor {
                ref_id: 0,
                ref_start: 10,
                query_start: 10,
            },
            Anchor {
                ref_id: 0,
                ref_start: 20,
                query_start: 20,
            },
            Anchor {
                ref_id: 0,
                ref_start: 30,
                query_start: 30,
            },
            Anchor {
                ref_id: 0,
                ref_start: 40,
                query_start: 40,
            },
            Anchor {
                ref_id: 0,
                ref_start: 50,
                query_start: 50,
            },
            Anchor {
                ref_id: 0,
                ref_start: 60,
                query_start: 60,
            },
            Anchor {
                ref_id: 0,
                ref_start: 70,
                query_start: 70,
            },
            Anchor {
                ref_id: 0,
                ref_start: 80,
                query_start: 80,
            },
            Anchor {
                ref_id: 0,
                ref_start: 90,
                query_start: 90,
            },
            Anchor {
                ref_id: 0,
                ref_start: 100,
                query_start: 100,
            },
        ];
        let mut result = expected.clone();
        remove_spurious_anchors(&mut result);
        assert_eq!(expected, result);
    }

    #[test]
    fn test_remove_spurious_anchors_inside_trim() {
        let mut result = vec![
            Anchor {
                ref_id: 0,
                ref_start: 0,
                query_start: 0,
            },
            Anchor {
                ref_id: 0,
                ref_start: 10,
                query_start: 10,
            },
            Anchor {
                ref_id: 0,
                ref_start: 20,
                query_start: 20,
            },
            Anchor {
                ref_id: 0,
                ref_start: 130,
                query_start: 30,
            },
            Anchor {
                ref_id: 0,
                ref_start: 140,
                query_start: 40,
            },
            Anchor {
                ref_id: 0,
                ref_start: 150,
                query_start: 50,
            },
            Anchor {
                ref_id: 0,
                ref_start: 60,
                query_start: 60,
            },
            Anchor {
                ref_id: 0,
                ref_start: 70,
                query_start: 50,
            },
            Anchor {
                ref_id: 0,
                ref_start: 80,
                query_start: 80,
            },
            Anchor {
                ref_id: 0,
                ref_start: 90,
                query_start: 90,
            },
            Anchor {
                ref_id: 0,
                ref_start: 100,
                query_start: 100,
            },
        ];
        remove_spurious_anchors(&mut result);
        let expected = vec![
            Anchor {
                ref_id: 0,
                ref_start: 0,
                query_start: 0,
            },
            Anchor {
                ref_id: 0,
                ref_start: 10,
                query_start: 10,
            },
            Anchor {
                ref_id: 0,
                ref_start: 20,
                query_start: 20,
            },
            Anchor {
                ref_id: 0,
                ref_start: 60,
                query_start: 60,
            },
            Anchor {
                ref_id: 0,
                ref_start: 80,
                query_start: 80,
            },
            Anchor {
                ref_id: 0,
                ref_start: 90,
                query_start: 90,
            },
            Anchor {
                ref_id: 0,
                ref_start: 100,
                query_start: 100,
            },
        ];
        assert_eq!(expected, result);
    }

    #[test]
    fn test_remove_spurious_anchors_ends_trim() {
        let mut result = vec![
            Anchor {
                ref_id: 0,
                ref_start: 5,
                query_start: 0,
            },
            Anchor {
                ref_id: 0,
                ref_start: 15,
                query_start: 10,
            },
            Anchor {
                ref_id: 0,
                ref_start: 20,
                query_start: 20,
            },
            Anchor {
                ref_id: 0,
                ref_start: 30,
                query_start: 30,
            },
            Anchor {
                ref_id: 0,
                ref_start: 40,
                query_start: 40,
            },
            Anchor {
                ref_id: 0,
                ref_start: 50,
                query_start: 50,
            },
            Anchor {
                ref_id: 0,
                ref_start: 60,
                query_start: 60,
            },
            Anchor {
                ref_id: 0,
                ref_start: 70,
                query_start: 70,
            },
            Anchor {
                ref_id: 0,
                ref_start: 80,
                query_start: 80,
            },
            Anchor {
                ref_id: 0,
                ref_start: 90,
                query_start: 95,
            },
            Anchor {
                ref_id: 0,
                ref_start: 100,
                query_start: 105,
            },
        ];
        remove_spurious_anchors(&mut result);
        let expected = vec![
            Anchor {
                ref_id: 0,
                ref_start: 20,
                query_start: 20,
            },
            Anchor {
                ref_id: 0,
                ref_start: 30,
                query_start: 30,
            },
            Anchor {
                ref_id: 0,
                ref_start: 40,
                query_start: 40,
            },
            Anchor {
                ref_id: 0,
                ref_start: 50,
                query_start: 50,
            },
            Anchor {
                ref_id: 0,
                ref_start: 60,
                query_start: 60,
            },
            Anchor {
                ref_id: 0,
                ref_start: 70,
                query_start: 70,
            },
            Anchor {
                ref_id: 0,
                ref_start: 80,
                query_start: 80,
            },
        ];
        assert_eq!(expected, result);
    }

    #[test]
    fn test_extend_piecewise_matching() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            5,
            100,
        );
        let query = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let refseq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let mut anchors = vec![
            Anchor {
                ref_id: 0,
                ref_start: 40,
                query_start: 30,
            },
            Anchor {
                ref_id: 0,
                ref_start: 30,
                query_start: 20,
            },
            Anchor {
                ref_id: 0,
                ref_start: 20,
                query_start: 10,
            },
            Anchor {
                ref_id: 0,
                ref_start: 10,
                query_start: 0,
            },
        ];
        let result = aligner
            .extend_piecewise(query, refseq, &mut anchors, 5)
            .unwrap();
        assert_eq!(result.score, 50 * 2 + 10 * 2);
        assert_eq!(result.edit_distance, 0);
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 50);
        assert_eq!(result.ref_start, 10);
        assert_eq!(result.ref_end, 60);
        assert_eq!(result.cigar.to_string(), "50=");
    }

    #[test]
    fn test_extend_piecewise_unmappable() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            5,
            100,
        );
        let query = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let refseq = b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        let mut anchors = vec![
            Anchor {
                ref_id: 0,
                ref_start: 30,
                query_start: 20,
            },
            Anchor {
                ref_id: 0,
                ref_start: 20,
                query_start: 20,
            },
            Anchor {
                ref_id: 0,
                ref_start: 10,
                query_start: 10,
            },
        ];
        let result = aligner.extend_piecewise(query, refseq, &mut anchors, 5);
        assert!(result.is_none());
    }

    #[test]
    fn test_extend_piecewise_overlapping_anchors() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            5,
            100,
        );
        let query = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let refseq =
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let mut anchors = vec![
            Anchor {
                ref_id: 0,
                ref_start: 40,
                query_start: 30,
            },
            Anchor {
                ref_id: 0,
                ref_start: 37,
                query_start: 27,
            },
            Anchor {
                ref_id: 0,
                ref_start: 35,
                query_start: 25,
            },
            Anchor {
                ref_id: 0,
                ref_start: 30,
                query_start: 20,
            },
            Anchor {
                ref_id: 0,
                ref_start: 25,
                query_start: 20,
            },
            Anchor {
                ref_id: 0,
                ref_start: 20,
                query_start: 15,
            },
            Anchor {
                ref_id: 0,
                ref_start: 10,
                query_start: 5,
            },
        ];
        let result = aligner
            .extend_piecewise(query, refseq, &mut anchors, 5)
            .unwrap();
        assert_eq!(result.score, 50 * 2 + 10 * 2 - 12 - 4 * 1);
        assert_eq!(result.edit_distance, 5);
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 50);
        assert_eq!(result.ref_start, 5);
        assert_eq!(result.ref_end, 60);
        assert_eq!(result.cigar.to_string(), "25=5D25=");
    }

    #[test]
    fn test_extend_piecewise_complex() {
        let aligner = PiecewiseAligner::new(
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: 12,
                gap_extend: 1,
                end_bonus: 10,
            },
            5,
            100,
        );
        let query = b"CTTTTAAAAATTTTAAAAATGGTTTCAAAAATTCCTAAAAATTTTTCCCCC";
        let refseq = b"TTTTTAAAAATTTTTAAAAATTTTTAAAAATTTTTAAAAATTTTTAAAAA";
        let mut anchors = vec![
            Anchor {
                ref_id: 0,
                ref_start: 35,
                query_start: 36,
            },
            Anchor {
                ref_id: 0,
                ref_start: 15,
                query_start: 14,
            },
            Anchor {
                ref_id: 0,
                ref_start: 5,
                query_start: 5,
            },
        ];
        let result = aligner
            .extend_piecewise(query, refseq, &mut anchors, 5)
            .unwrap();
        assert_eq!(
            result.score,
            (13 + 6 + 3 + 7 + 11) * 2 + 10 - 12 * 2 - 1 - 8 * 4
        );
        assert_eq!(result.edit_distance, 7);
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 46);
        assert_eq!(result.ref_start, 0);
        assert_eq!(result.ref_end, 45);
        assert_eq!(result.cigar.to_string(), "1X13=1D6=2I3=1X7=2X11=");
    }
}
