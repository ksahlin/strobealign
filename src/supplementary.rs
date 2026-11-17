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

    let primary_start = primary.soft_clip_left;
    let primary_end = read.len() - primary.soft_clip_right;
    for chain in candidates.iter_mut().take(parameters.max_tries) {
        if chain.score < dropoff {
            break;
        }

        let interval = if chain.query_start >= primary_start {
            primary_end..read.len()
        } else {
            0..primary_start
        };
        if interval.is_empty() {
            continue;
        }
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
/// the same contig and strand as the primary, and reaches back into the
/// primary on the read by at most `overlap_fraction` of the primary's read
/// span.
fn is_candidate(
    primary: &Alignment,
    refseq: &RefSequence,
    chain: &Chain,
    read_len: usize,
    overlap_fraction: f32,
) -> bool {
    let (contig_id, _) = refseq.unflatten(chain.ref_start);
    if primary.is_revcomp != chain.is_revcomp || primary.contig_id != contig_id {
        return false;
    }

    let read_gap = if chain.query_start >= primary.soft_clip_left {
        chain.query_start as i64 - (read_len - primary.soft_clip_right) as i64
    } else {
        primary.soft_clip_left as i64 - chain.query_end as i64
    };

    let primary_span = (read_len - primary.soft_clip_left - primary.soft_clip_right) as i64;
    let overlap_tolerance = (primary_span as f32 * overlap_fraction) as i64;

    read_gap >= -overlap_tolerance
}

/// Returns the gaps between two alignments on the read and on the reference,
/// as (read_gap, ref_gap). A negative gap means the alignments overlap.
fn alignment_gaps(primary: &Alignment, candidate: &Alignment, read_len: usize) -> (i64, i64) {
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
