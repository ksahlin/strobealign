use crate::aligner::{Aligner, Scores};
use crate::cigar::{Cigar, CigarOperation};
use crate::modes::align::Alignment;
use crate::packed_seq::PackedSeqSlice;
use crate::read::Read;
use crate::refseq::RefSequence;
use crate::simdaligner::SplitReferenceAlignment;

/// Cuts the alignment into one segment per copy wherever an insertion of at least `min_length`
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

    let mut segments = Vec::new();
    let mut segment_cigar = Cigar::new();
    let mut segment_start = alignment.contig_start;
    let mut segment_clip_left = alignment.soft_clip_left;
    let mut ref_pos = alignment.contig_start;
    let mut query_pos = alignment.soft_clip_left;

    for (op, len) in alignment.cigar.iter() {
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
            )
        {
            // The right arm carries the segment being built on past the insertion point, and
            // the left arm opens the next one at the second copy.
            let mut upstream = segment_cigar.clone();
            upstream.extend(&split.right.cigar);

            if let Some(segment) = clip_segment(
                &alignment,
                upstream,
                segment_start,
                segment_clip_left,
                query.len(),
                &aligner.scores,
            ) {
                segments.push(segment);
                segment_cigar = split.left.cigar;
                segment_start = copy_start;
                segment_clip_left = query_pos + split.left.query_start;
                query_pos += len;
                continue;
            }
        }

        segment_cigar.push(op, len);
        let (query_len, ref_len) = query_ref_len(op, len);
        query_pos += query_len;
        ref_pos += ref_len;
    }

    if segments.is_empty() {
        return vec![alignment];
    }

    let Some(last) = clip_segment(
        &alignment,
        segment_cigar,
        segment_start,
        segment_clip_left,
        query.len(),
        &aligner.scores,
    ) else {
        return vec![alignment];
    };
    segments.push(last);

    segments
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
) -> Option<(usize, SplitReferenceAlignment)> {
    let len = insertion.len();
    let padding = len + band;
    let left_start = ref_pos.saturating_sub(len + padding);
    let right_end = (ref_pos + len + padding).min(contig.len());
    let left_reference = contig.decode(left_start, ref_pos);
    let right_reference = contig.decode(ref_pos, right_end);

    let split = aligner.realign_insertion(insertion, &left_reference, &right_reference, band)?;

    let matched = split.left.cigar.matches() + split.right.cigar.matches();
    let copy_len =
        (split.right.ref_end - split.right.ref_start) + (split.left.ref_end - split.left.ref_start);

    // The duplication must match most of the inserted sequence.
    if (matched as f32) < min_matches * len as f32 || copy_len < min_length {
        return None;
    }

    Some((left_start + split.left.ref_start, split))
}

/// Builds the alignment of one segment from the cigar it covers, soft-clipping the read on
/// both sides and scoring the segment on its own. None if it scores zero or less.
fn clip_segment(
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

    for (op, len) in cigar.iter() {
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
