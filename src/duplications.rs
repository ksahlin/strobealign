use crate::aligner::{Aligner, Scores};
use crate::cigar::{Cigar, CigarOperation, OpLen};
use crate::details::AlignmentDetails;
use crate::modes::align::Alignment;
use crate::packed_seq::PackedSeqSlice;
use crate::read::Read;
use crate::refseq::RefSequence;

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

    let mut segments = Vec::new();
    let mut segment_cigar = Cigar::new();
    let mut segment_start = alignment.contig_start;
    let mut segment_clip_left = alignment.soft_clip_left;
    let mut ref_pos = alignment.contig_start;
    let mut query_pos = alignment.soft_clip_left;

    for OpLen { op, len } in alignment.cigar.iter() {
        if op == CigarOperation::Insertion
            && len >= min_length
            && let Some(duplication) = find_duplication(
                aligner,
                &contig,
                &query[query_pos..query_pos + len],
                ref_pos,
                min_length,
                min_matches,
                band,
                details,
            )
            && let Some(inner_segments) = build_inner_segments(
                &alignment,
                &duplication,
                query_pos,
                query.len(),
                &aligner.scores,
            )
        {
            // The head arm carries the segment being built on past the insertion point, and
            // the tail arm opens the next one at the last copy.
            let mut upstream = segment_cigar.clone();
            upstream.extend(&duplication.head.cigar);

            if let Some(segment) = build_segment(
                &alignment,
                upstream,
                segment_start,
                segment_clip_left,
                query.len(),
                &aligner.scores,
            ) {
                segments.push(segment);
                segments.extend(inner_segments);
                segment_start = duplication.tail.contig_start;
                segment_clip_left = query_pos + duplication.tail.query_start;
                segment_cigar = duplication.tail.cigar;
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

    let Some(last) = build_segment(
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
    details.duplications += segments.len() - 1;

    segments
}

struct Duplication {
    head: DuplicatedCopy,
    inner: Vec<DuplicatedCopy>,
    tail: DuplicatedCopy,
}

struct DuplicatedCopy {
    contig_start: usize,
    query_start: usize,
    cigar: Cigar,
}

/// None when the inserted bases are not a tandem duplication. A tiling that is not good enough
/// falls through to the split rather than ending the search.
fn find_duplication(
    aligner: &Aligner,
    contig: &PackedSeqSlice,
    insertion: &[u8],
    ref_pos: usize,
    min_length: usize,
    min_matches: f32,
    band: usize,
    details: &mut AlignmentDetails,
) -> Option<Duplication> {
    details.duplication_candidates += 1;
    let len = insertion.len();
    let padding = len + band;
    let left_start = ref_pos.saturating_sub(len + padding);
    let right_end = (ref_pos + len + padding).min(contig.len());
    let left_reference = contig.decode(left_start, ref_pos);
    let right_reference = contig.decode(ref_pos, right_end);

    if let Some(duplication) = find_multi_copy(
        aligner,
        insertion,
        &left_reference,
        &right_reference,
        ref_pos,
        min_length,
        min_matches,
        band,
    ) {
        return Some(duplication);
    }

    find_single_copy(
        aligner,
        insertion,
        &left_reference,
        &right_reference,
        ref_pos,
        min_length,
        min_matches,
        band,
    )
}

fn find_single_copy(
    aligner: &Aligner,
    insertion: &[u8],
    left_reference: &[u8],
    right_reference: &[u8],
    ref_pos: usize,
    min_length: usize,
    min_matches: f32,
    band: usize,
) -> Option<Duplication> {
    let split = aligner.realign_insertion(insertion, left_reference, right_reference, band)?;
    let copy_len =
        (split.head.ref_end - split.head.ref_start) + (split.tail.ref_end - split.tail.ref_start);
    if copy_len < min_length {
        return None;
    }

    // Matched bases over the insertion length is a cheap proxy for the identity of the copies.
    let matched = split.head.cigar.matches() + split.tail.cigar.matches();
    if (matched as f32) < min_matches * insertion.len() as f32 {
        return None;
    }

    Some(Duplication {
        head: DuplicatedCopy {
            contig_start: ref_pos,
            query_start: 0,
            cigar: split.head.cigar,
        },
        inner: vec![],
        tail: DuplicatedCopy {
            contig_start: ref_pos - left_reference.len() + split.tail.ref_start,
            query_start: split.tail.query_start,
            cigar: split.tail.cigar,
        },
    })
}

fn find_multi_copy(
    aligner: &Aligner,
    insertion: &[u8],
    left_reference: &[u8],
    right_reference: &[u8],
    ref_pos: usize,
    min_length: usize,
    min_matches: f32,
    band: usize,
) -> Option<Duplication> {
    // Not enough room for multiple copies
    if 2 * insertion.len() < 3 * min_length {
        return None;
    }

    // Find approximate unit length (original copy size)
    let (backward, forward) =
        aligner.insertion_reach(insertion, left_reference, right_reference, band)?;
    let unit_len = forward + backward;
    if unit_len < min_length {
        return None;
    }

    // Estimate copy numbers
    let copy_count = (insertion.len() as f64 / unit_len as f64).round() as usize;
    if copy_count < 2 {
        return None;
    }

    // Creating the tiled_reference: (prefix + suffix) * copy_count.
    let unit_prefix = &right_reference[..forward];
    let unit_suffix = &left_reference[left_reference.len() - backward..];
    let mut tiled_reference = Vec::with_capacity(unit_len * copy_count);
    for _ in 0..copy_count {
        tiled_reference.extend_from_slice(unit_prefix);
        tiled_reference.extend_from_slice(unit_suffix);
    }

    // Aligning the insertion against the tiled reference
    let tiled_alignment = aligner.realign_copies(insertion, &tiled_reference, band)?;

    // Cut the cigar of the tiled_alignment to get the cigar of each copy
    let mut copy_bounds = Vec::with_capacity(copy_count);
    for index in 0..copy_count {
        copy_bounds.push(forward + index * unit_len);
    }
    let copy_cigars = cut_copy_cigars(tiled_alignment.cigar, &copy_bounds);

    // Separate first, last and inner copies
    let [first, inner_cigars @ .., last] = &copy_cigars[..] else {
        return None;
    };

    // Start of the unit on the reference
    let unit_start = ref_pos - backward;

    // Duplication head, deletions trimmed at it's end
    let mut cigar = first.clone();
    if let Some(last_op) = cigar.last()
        && last_op.op == CigarOperation::Deletion
    {
        cigar.pop();
    }
    let (mut query_start, _) = query_ref_span(&cigar);
    let head = DuplicatedCopy {
        contig_start: ref_pos,
        query_start: 0,
        cigar,
    };

    // Inner duplications, deletions trimmed on both ends
    let mut inner = Vec::with_capacity(inner_cigars.len());
    for copy_cigar in inner_cigars {
        let mut cigar = copy_cigar.clone();
        if let Some(last_op) = cigar.last()
            && last_op.op == CigarOperation::Deletion
        {
            cigar.pop();
        }

        let mut ref_skip = 0;
        if let Some(first_op) = cigar.first()
            && first_op.op == CigarOperation::Deletion
        {
            cigar.pop_front();
            ref_skip = first_op.len;
        }
        let (query_span, _) = query_ref_span(&cigar);

        if !cigar.is_empty() {
            inner.push(DuplicatedCopy {
                contig_start: unit_start + ref_skip,
                query_start,
                cigar,
            });
        }

        query_start += query_span;
    }

    // Duplication tail, deletions trimmed at it's front
    let mut cigar = last.clone();
    let mut ref_skip = 0;
    if let Some(first_op) = cigar.first()
        && first_op.op == CigarOperation::Deletion
    {
        cigar.pop_front();
        ref_skip = first_op.len;
    }
    let tail = DuplicatedCopy {
        contig_start: unit_start + ref_skip,
        query_start,
        cigar,
    };

    // Matched bases over the insertion length is a cheap proxy for the identity of the copies.
    let mut matched = head.cigar.matches() + tail.cigar.matches();
    for copy in &inner {
        matched += copy.cigar.matches();
    }
    if (matched as f32) < min_matches * insertion.len() as f32 {
        return None;
    }

    Some(Duplication { head, inner, tail })
}

fn cut_copy_cigars(cigar: Cigar, copy_bounds: &[usize]) -> Vec<Cigar> {
    debug_assert!(copy_bounds.is_sorted());
    let mut copy_cigars = Vec::with_capacity(copy_bounds.len() + 1);
    let mut copy_start = 0;
    let mut rest = cigar;

    for &ref_bound in copy_bounds {
        let (copy_cigar, remainder) = split_cigar_at(&rest, ref_bound - copy_start);
        copy_cigars.push(copy_cigar);
        copy_start = ref_bound;
        rest = remainder;
    }
    copy_cigars.push(rest);

    copy_cigars
}

fn split_cigar_at(cigar: &Cigar, ref_bound: usize) -> (Cigar, Cigar) {
    let mut before = Cigar::new();
    let mut after = Cigar::new();
    let mut consumed = 0;

    for OpLen { op, len } in cigar.iter() {
        let (_, ref_len) = query_ref_len(op, len);
        // Insertions do no advance on the reference
        if ref_len == 0 {
            if after.is_empty() {
                before.push(op, len);
            } else {
                after.push(op, len);
            }
            continue;
        }
        // split the cigar between what goes before and after
        let taken = len.min(ref_bound.saturating_sub(consumed));
        before.push(op, taken);
        after.push(op, len - taken);
        consumed += len;
    }

    (before, after)
}

/// Builds a segment for every copy between the head and the tail of a multi
/// copy duplication. None is one of them has a score <= 0.
fn build_inner_segments(
    source: &Alignment,
    duplication: &Duplication,
    query_pos: usize,
    read_len: usize,
    scores: &Scores,
) -> Option<Vec<Alignment>> {
    let mut segments = Vec::with_capacity(duplication.inner.len());
    for copy in &duplication.inner {
        segments.push(build_segment(
            source,
            copy.cigar.clone(),
            copy.contig_start,
            query_pos + copy.query_start,
            read_len,
            scores,
        )?);
    }

    Some(segments)
}

/// Makes an alignment for a segment of a duplication split. It needs its own
/// cigar, reference start, soft clips, score and edit distance, the rest comes
/// from the original alignment. Returns None if it aligns nothing or scores zero or less.
fn build_segment(
    source: &Alignment,
    cigar: Cigar,
    ref_start: usize,
    soft_clip_left: usize,
    read_len: usize,
    scores: &Scores,
) -> Option<Alignment> {
    if cigar.is_empty() {
        return None;
    }
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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn cut_copy_cigars_bounds() {
        let cigars = cut_copy_cigars("20=".parse().unwrap(), &[5, 15]);
        assert_eq!(cigars.len(), 3);
        assert_eq!(cigars[0].to_string(), "5=");
        assert_eq!(cigars[1].to_string(), "10=");
        assert_eq!(cigars[2].to_string(), "5=");

        let cigars = cut_copy_cigars("10=5I10=".parse().unwrap(), &[10]);
        assert_eq!(cigars[0].to_string(), "10=5I");
        assert_eq!(cigars[1].to_string(), "10=");

        let cigars = cut_copy_cigars("10=".parse().unwrap(), &[10, 20]);
        assert_eq!(cigars.len(), 3);
        assert!(cigars[1].is_empty());
        assert!(cigars[2].is_empty());
    }
}
