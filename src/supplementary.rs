use std::ops::Range;

use crate::aligner::Aligner;
use crate::chain::Chain;
use crate::details::AlignmentDetails;
use crate::modes::align::{
    Alignment, AlignmentWithMapq, MappingParameters, extend_seed_in_interval,
};
use crate::read::Read;
use crate::refseq::RefSequence;

/// Extend the chains that may continue the primary alignment over the part of
/// the read it does not cover and return the best one. Returns `None` if no
/// chain is a candidate, or if the best extension scores below
/// `min_supplementary_score`.
pub fn get_supplementary_alignment(
    aligner: &Aligner,
    refseq: &RefSequence,
    read: &Read,
    k: usize,
    primary: &Alignment,
    chains: &mut [Chain],
    parameters: &MappingParameters,
    details: &mut AlignmentDetails,
) -> Option<AlignmentWithMapq> {
    let mut candidates: Vec<&mut Chain> = chains
        .iter_mut()
        .filter(|chain| {
            is_candidate(
                primary,
                refseq,
                chain,
                read.len(),
                parameters.max_supplementary_overlap,
            )
        })
        .collect();
    details.supplementary_candidates += candidates.len();

    let dropoff = candidates.first()?.score * parameters.dropoff_threshold;
    let max_gap = parameters.max_supplementary_gap as i64;
    let mut best_alignment = None;
    let mut best_score = 0;
    let mut second_best_score = 0;

    let primary_span = alignment_read_span(primary, read.len());
    for chain in candidates.iter_mut().take(parameters.max_tries) {
        if chain.score < dropoff {
            break;
        }

        let interval = if chain_read_span(chain, read.len()).start >= primary_span.start {
            primary_span.end..read.len()
        } else {
            0..primary_span.start
        };
        if interval.is_empty() {
            continue;
        }
        // The interval is extended on the chain's strand
        let interval = if chain.is_revcomp {
            read.len() - interval.end..read.len() - interval.start
        } else {
            interval
        };
        let Some(alignment) = extend_seed_in_interval(
            aligner,
            chain,
            refseq,
            read,
            parameters.use_ssw,
            k,
            interval,
        ) else {
            continue;
        };
        details.supplementary_extended += 1;
        details.add_alignment(&alignment);

        let (read_gap, ref_gap) = alignment_gaps(primary, &alignment, read.len());
        if read_gap.abs() > max_gap && ref_gap.abs() > max_gap {
            continue;
        }

        if alignment.score > best_score {
            second_best_score = best_score;
            best_score = alignment.score;
            best_alignment = Some(alignment);
        } else if alignment.score > second_best_score {
            second_best_score = alignment.score;
        }
    }

    if best_score < parameters.min_supplementary_score {
        return None;
    }

    let alignment = best_alignment?;
    let mapq = (60 * best_score.saturating_sub(second_best_score)).div_ceil(best_score) as u8;
    details.supplementary += 1;

    Some(AlignmentWithMapq { alignment, mapq })
}

/// Whether a chain is worth extending into a supplementary alignment: it is on
/// the same contig as the primary, and reaches back into the primary on the
/// read by at most `overlap_fraction` of the primary's read span. A chain on
/// the other strand may be an inversion.
fn is_candidate(
    primary: &Alignment,
    refseq: &RefSequence,
    chain: &Chain,
    read_len: usize,
    overlap_fraction: f32,
) -> bool {
    let (contig_id, _) = refseq.unflatten(chain.ref_start);
    if primary.contig_id != contig_id {
        return false;
    }

    let primary_span = alignment_read_span(primary, read_len);
    let chain_span = chain_read_span(chain, read_len);
    let read_gap = if chain_span.start >= primary_span.start {
        chain_span.start as i64 - primary_span.end as i64
    } else {
        primary_span.start as i64 - chain_span.end as i64
    };

    let overlap_tolerance = (primary_span.len() as f32 * overlap_fraction) as i64;

    read_gap >= -overlap_tolerance
}

/// The part of the original read a chain covers.
fn chain_read_span(chain: &Chain, read_len: usize) -> Range<usize> {
    if chain.is_revcomp {
        read_len - chain.query_end..read_len - chain.query_start
    } else {
        chain.query_start..chain.query_end
    }
}

/// The part of the original read an alignment covers.
fn alignment_read_span(alignment: &Alignment, read_len: usize) -> Range<usize> {
    if alignment.is_revcomp {
        alignment.soft_clip_right..read_len - alignment.soft_clip_left
    } else {
        alignment.soft_clip_left..read_len - alignment.soft_clip_right
    }
}

/// Returns the gaps between two alignments on the read and on the reference,
/// as (read_gap, ref_gap). A negative gap means the alignments overlap.
///
/// If the alignments are on different strands, an inversion may join them at
/// either of their ends on the reference, so the reference gap is the smallest
/// distance between an end of one and an end of the other.
fn alignment_gaps(primary: &Alignment, candidate: &Alignment, read_len: usize) -> (i64, i64) {
    if candidate.is_revcomp != primary.is_revcomp {
        let primary_span = alignment_read_span(primary, read_len);
        let candidate_span = alignment_read_span(candidate, read_len);
        let read_gap = if candidate_span.start >= primary_span.start {
            candidate_span.start as i64 - primary_span.end as i64
        } else {
            primary_span.start as i64 - candidate_span.end as i64
        };
        let primary_end = primary.contig_start + primary.length;
        let candidate_end = candidate.contig_start + candidate.length;
        let ref_gap = primary
            .contig_start
            .abs_diff(candidate.contig_start)
            .min(primary.contig_start.abs_diff(candidate_end))
            .min(primary_end.abs_diff(candidate.contig_start))
            .min(primary_end.abs_diff(candidate_end));
        return (read_gap, ref_gap as i64);
    }

    if candidate.soft_clip_left >= primary.soft_clip_left {
        (
            candidate.soft_clip_left as i64 - (read_len - primary.soft_clip_right) as i64,
            candidate.contig_start as i64 - (primary.contig_start + primary.length) as i64,
        )
    } else {
        (
            primary.soft_clip_left as i64 - (read_len - candidate.soft_clip_right) as i64,
            primary.contig_start as i64 - (candidate.contig_start + candidate.length) as i64,
        )
    }
}
