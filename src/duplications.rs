use std::cmp::Reverse;

use crate::aligner::{Aligner, Scores};
use crate::cigar::{Cigar, CigarOperation, OpLen};
use crate::details::AlignmentDetails;
use crate::modes::align::{Alignment, AlignmentWithMapq, MappingParameters};
use crate::packed_seq::PackedSeqSlice;
use crate::read::Read;
use crate::refseq::RefSequence;
use crate::simdaligner::SplitReferenceAlignment;

/// Splits each alignment that contains a tandem duplication into one alignment per copy. Each
/// copy keeps the mapping quality of the alignment it came from.
///
/// The returned alignments keep the order of the input alignments: the copies of the first
/// input alignment come first, then those of the second, and so on. The copies of one input
/// alignment are sorted by score, best first. With the primary alignment passed first, the
/// first returned alignment is therefore the best copy of the primary.
pub fn split_duplications(
    aligner: &Aligner,
    refseq: &RefSequence,
    read: &Read,
    alignments: Vec<AlignmentWithMapq>,
    parameters: &MappingParameters,
    details: &mut AlignmentDetails,
) -> Vec<AlignmentWithMapq> {
    let mut resolved = Vec::new();
    for AlignmentWithMapq { alignment, mapq } in alignments {
        let mut split = reinterpret_insertions(
            aligner,
            alignment,
            refseq,
            read,
            parameters.min_duplication_length,
            parameters.min_duplication_matches,
            parameters.duplication_band,
            details,
        );

        // Sort only within one split, so that the copies of the primary alignment stay ahead of
        // those of the supplementary alignment
        split.sort_by_key(|alignment| Reverse(alignment.score));
        resolved.extend(
            split
                .into_iter()
                .map(|alignment| AlignmentWithMapq { alignment, mapq }),
        );
    }

    resolved
}

/// Cuts the alignment into one alignment per copy wherever an insertion of at least `min_length`
/// bases turns out to be a tandem duplication. A `min_length` of zero turns this off. Always
/// returns at least the alignment it was given.
pub fn reinterpret_insertions(
    aligner: &Aligner,
    alignment: Alignment,
    refseq: &RefSequence,
    read: &Read,
    min_length: usize,
    min_matches: f32,
    band: usize,
    details: &mut AlignmentDetails,
) -> Vec<Alignment> {
    if min_length == 0 {
        return vec![alignment];
    }
    let query = if alignment.is_revcomp {
        read.rc()
    } else {
        read.seq()
    };
    let contig = refseq.contig(alignment.contig_id);

    let mut alignments = Vec::new();
    let mut current_cigar = Cigar::new();
    let mut current_start = alignment.contig_start;
    let mut current_clip_left = alignment.soft_clip_left;
    let mut ref_pos = alignment.contig_start;
    let mut query_pos = alignment.soft_clip_left;

    for OpLen { op, len } in alignment.cigar.iter() {
        if op == CigarOperation::Insertion
            && len >= min_length
            && let Some((copy_start, split)) = find_duplication(
                aligner,
                &contig,
                &query[query_pos..query_pos + len],
                ref_pos,
                min_length,
                min_matches,
                band,
                details,
            )
        {
            // The head arm carries the alignment being built on past the insertion point, and
            // the tail arm opens the next one at the second copy.
            let mut upstream = current_cigar.clone();
            upstream.extend(&split.head.cigar);

            if let Some(head) = build_alignment(
                &alignment,
                upstream,
                current_start,
                current_clip_left,
                query.len(),
                &aligner.scores,
            ) {
                alignments.push(head);
                current_cigar = split.tail.cigar;
                current_start = copy_start;
                current_clip_left = query_pos + split.tail.query_start;
                query_pos += len;
                continue;
            }
        }

        current_cigar.push(op, len);
        let (query_len, ref_len) = query_ref_len(op, len);
        query_pos += query_len;
        ref_pos += ref_len;
    }

    if alignments.is_empty() {
        return vec![alignment];
    }

    let Some(last) = build_alignment(
        &alignment,
        current_cigar,
        current_start,
        current_clip_left,
        query.len(),
        &aligner.scores,
    ) else {
        return vec![alignment];
    };
    alignments.push(last);
    details.duplications += alignments.len() - 1;

    alignments
}

/// Aligns the bases inserted at `ref_pos` against both references surrounding the
/// insertion region and returns a duplication pattern: 2 alignments bridging the
/// query region by jumping from one reference to the other. None when the
/// duplication pattern is not good enough to report.
fn find_duplication(
    aligner: &Aligner,
    contig: &PackedSeqSlice,
    insertion: &[u8],
    ref_pos: usize,
    min_length: usize,
    min_matches: f32,
    band: usize,
    details: &mut AlignmentDetails,
) -> Option<(usize, SplitReferenceAlignment)> {
    details.duplication_candidates += 1;
    let len = insertion.len();
    let padding = len + band;
    let left_start = ref_pos.saturating_sub(len + padding);
    let right_end = (ref_pos + len + padding).min(contig.len());
    let left_reference = contig.decode(left_start, ref_pos);
    let right_reference = contig.decode(ref_pos, right_end);

    let split = aligner.realign_insertion(insertion, &left_reference, &right_reference, band)?;

    let matched = split.head.cigar.matches() + split.tail.cigar.matches();
    let copy_len =
        (split.head.ref_end - split.head.ref_start) + (split.tail.ref_end - split.tail.ref_start);

    // The duplication must match most of the inserted sequence: matched bases over the
    // insertion length is a cheap proxy for the identity of the two copies.
    if (matched as f32) < min_matches * len as f32 || copy_len < min_length {
        return None;
    }

    Some((left_start + split.tail.ref_start, split))
}

/// Makes an alignment for one copy of a duplication split. It needs its own
/// cigar, reference start, soft clips, score and edit distance, the rest comes
/// from the original alignment. Returns None if it scores zero or less.
fn build_alignment(
    source: &Alignment,
    cigar: Cigar,
    ref_start: usize,
    soft_clip_left: usize,
    read_len: usize,
    scores: &Scores,
) -> Option<Alignment> {
    let (query_span, ref_span) = query_ref_span(&cigar);
    let soft_clip_right = read_len - (soft_clip_left + query_span);
    let unclipped_ends = (soft_clip_left == 0) as u32 + (soft_clip_right == 0) as u32;
    let score = scores.score(&cigar) + (unclipped_ends * scores.end_bonus) as i32;
    if score <= 0 {
        return None;
    }

    Some(Alignment {
        contig_id: source.contig_id,
        is_revcomp: source.is_revcomp,
        contig_start: ref_start,
        length: ref_span,
        soft_clip_left,
        soft_clip_right,
        edit_distance: cigar.edit_distance(),
        score: score as u32,
        cigar,
        gapped: source.gapped,
        anchors: source.anchors,
        collisions: source.collisions,
    })
}

fn query_ref_span(cigar: &Cigar) -> (usize, usize) {
    let mut query = 0;
    let mut reference = 0;

    for OpLen { op, len } in cigar.iter() {
        let (q, r) = query_ref_len(op, len);
        query += q;
        reference += r;
    }

    (query, reference)
}

fn query_ref_len(op: CigarOperation, len: usize) -> (usize, usize) {
    match op {
        CigarOperation::Match | CigarOperation::Eq | CigarOperation::X => (len, len),
        CigarOperation::Insertion | CigarOperation::Softclip => (len, 0),
        CigarOperation::Deletion | CigarOperation::Skip => (0, len),
        CigarOperation::Hardclip | CigarOperation::Pad => (0, 0),
    }
}
