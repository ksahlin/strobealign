use crate::aligner::Aligner;
use crate::chain::Chain;
use crate::details::AlignmentDetails;
use crate::modes::align::{Alignment, MappingParameters, Segment, extend_seed};
use crate::read::Read;
use crate::refseq::RefSequence;

/// Extend the chains that fit the pattern of a large variation with the
/// primary and return the best one. Returns `None` if no chain fits or clears
/// the score floor.
pub fn get_supplementary_segment(
    aligner: &Aligner,
    refseq: &RefSequence,
    read: &Read,
    k: usize,
    primary: &Alignment,
    chains: &mut [Chain],
    parameters: &MappingParameters,
    details: &mut AlignmentDetails,
) -> Option<Segment> {
    let mut candidates: Vec<&mut Chain> = chains
        .iter_mut()
        .filter(|chain| {
            fits_variant_pattern(
                primary,
                refseq,
                chain,
                read.len(),
                parameters.max_chain_distance as i64,
                parameters.max_gap_distance as i64,
            )
        })
        .collect();
    details.supplementary_candidates += candidates.len();

    let dropoff = candidates.first()?.score * parameters.dropoff_threshold;
    let max_gap = parameters.max_segment_gap as i64;
    let mut best_alignment = None;
    let mut best_score = 0;
    let mut second_best_score = 0;

    for chain in candidates.iter_mut().take(parameters.max_tries) {
        if chain.score < dropoff {
            break;
        }
        let Some(alignment) = extend_seed(aligner, chain, refseq, read, parameters.use_ssw, k)
        else {
            continue;
        };
        details.supplementary_extended += 1;
        details.add_alignment(&alignment);

        let (read_gap, ref_gap) = segment_gaps(primary, &alignment, read.len());
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

    Some(Segment { alignment, mapq })
}

/// Whether a chain and the primary are the two sides of a deletion, insertion
/// or duplication.
fn fits_variant_pattern(
    primary: &Alignment,
    refseq: &RefSequence,
    chain: &Chain,
    read_len: usize,
    max_chain_distance: i64,
    max_gap_distance: i64,
) -> bool {
    let (contig_id, contig_start) = refseq.unflatten(chain.ref_start);
    let contig_start = contig_start.0 as usize;
    if primary.is_revcomp != chain.is_revcomp || primary.contig_id != contig_id {
        return false;
    }

    let (read_gap, ref_gap) = if chain.query_start >= primary.soft_clip_left {
        (
            chain.query_start as i64 - (read_len - primary.soft_clip_right) as i64,
            contig_start as i64 - (primary.contig_start + primary.length) as i64,
        )
    } else {
        (
            primary.soft_clip_left as i64 - chain.query_end as i64,
            primary.contig_start as i64 - (contig_start + chain.ref_span()) as i64,
        )
    };

    // A candidate has to bring read bases the primary does not have, which is also what
    // keeps the chain the primary was built from out of its own search
    if read_gap < 0 {
        return false;
    }

    // An insertion: the contig carries straight on while the read jumps the
    // inserted bases
    if (0..=max_chain_distance).contains(&ref_gap) {
        return true;
    }

    // A deletion/duplication: the read carries straight on while the contig
    // jumps bases
    if read_gap <= max_chain_distance && (-max_gap_distance..max_gap_distance).contains(&ref_gap) {
        return true;
    }

    false
}

/// How far apart two alignments end up on the read and on the reference.
fn segment_gaps(primary: &Alignment, candidate: &Alignment, read_len: usize) -> (i64, i64) {
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
