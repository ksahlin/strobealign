use std::cmp::Reverse;

use crate::aligner::Aligner;
use crate::chainer::Chainer;
use crate::details::{Details, NamDetails};
use crate::index::StrobemerIndex;
use crate::insertsize::InsertSizeDistribution;
use crate::io::fasta::RefSequence;
use crate::io::paf::PafRecord;
use crate::io::record::{End, SequenceRecord};
use crate::mapper::{Alignment, extend_seed, mapping_quality, rescue_align};
use crate::math::normal_pdf;
use crate::mcsstrategy::McsStrategy;
use crate::nam::{Nam, get_nams_by_chaining, reverse_nam_if_needed, sort_nams};
use crate::read::Read;
use crate::shuffle::shuffle_best;

use fastrand::Rng;

/// Map a single-end read to the reference and return PAF records
///
/// This implements "mapping-only" mode in which no base-level alignments are computed
pub fn map_single_end_read(
    record: &SequenceRecord,
    index: &StrobemerIndex,
    references: &[RefSequence],
    rescue_distance: usize,
    mcs_strategy: McsStrategy,
    chainer: &Chainer,
    rng: &mut Rng,
) -> (Vec<PafRecord>, Details) {
    let (mut nam_details, mut nams) = get_nams_by_chaining(
        &record.sequence,
        index,
        chainer,
        rescue_distance,
        mcs_strategy,
    );
    nam_details.time_sort_nams = sort_nams(&mut nams, rng);

    if nams.is_empty() {
        (vec![], nam_details.into())
    } else {
        let mapq = mapping_quality(&nams);
        (
            vec![paf_record_from_nam(
                &nams[0],
                &record.name,
                references,
                record.sequence.len(),
                Some(mapq),
                End::None,
            )],
            nam_details.into(),
        )
    }
}

/// Map a single-end read to the reference and estimate abundances
///
/// This implements abundance estimation mode (`--aemb`)
pub fn abundances_single_end_read(
    record: &SequenceRecord,
    index: &StrobemerIndex,
    abundances: &mut [f64],
    rescue_distance: usize,
    mcs_strategy: McsStrategy,
    chainer: &Chainer,
    rng: &mut Rng,
) {
    let (_, mut nams) = get_nams_by_chaining(
        &record.sequence,
        index,
        chainer,
        rescue_distance,
        mcs_strategy,
    );
    sort_nams(&mut nams, rng);

    let n_best = nams
        .iter()
        .take_while(|nam| nam.score == nams[0].score)
        .count();
    let weight = record.sequence.len() as f64 / n_best as f64;
    for nam in &nams[0..n_best] {
        abundances[nam.ref_id] += weight;
    }
}

/// Convert Nam into PAF record
fn paf_record_from_nam(
    nam: &Nam,
    name: &str,
    references: &[RefSequence],
    query_length: usize,
    mapq: Option<u8>,
    end: End,
) -> PafRecord {
    PafRecord {
        query_name: name.into(),
        end,
        query_length: query_length as u64,
        query_start: nam.query_start as u64,
        query_end: nam.query_end as u64,
        is_revcomp: nam.is_revcomp,
        target_name: references[nam.ref_id].name.clone(),
        target_length: references[nam.ref_id].sequence.len() as u64,
        target_start: nam.ref_start as u64,
        target_end: nam.ref_end as u64,
        matching_bases: nam.matching_bases as u64,
        alignment_length: (nam.ref_end - nam.ref_start) as u64,
        mapping_quality: mapq,
    }
}

/// Map a paired-end read pair to the reference and return PAF records
///
/// This implements "mapping-only" mode in which no base-level alignments are computed
pub fn map_paired_end_read(
    r1: &SequenceRecord,
    r2: &SequenceRecord,
    index: &StrobemerIndex,
    references: &[RefSequence],
    rescue_distance: usize,
    insert_size_distribution: &mut InsertSizeDistribution,
    mcs_strategy: McsStrategy,
    chainer: &Chainer,
    rng: &mut Rng,
) -> (Vec<PafRecord>, Details) {
    let (mut nam_details1, mut nams1) =
        get_nams_by_chaining(&r1.sequence, index, chainer, rescue_distance, mcs_strategy);
    let (mut nam_details2, mut nams2) =
        get_nams_by_chaining(&r2.sequence, index, chainer, rescue_distance, mcs_strategy);

    if nams1.is_empty() && nams2.is_empty() {
        nam_details1 += nam_details2;
        return (vec![], nam_details1.into());
    }

    let mut nam_pairs = get_nam_pairs(
        &mut nams1,
        &mut nams2,
        insert_size_distribution.mu,
        insert_size_distribution.sigma,
        &nam_details1,
        &nam_details2,
    );
    shuffle_best(&mut nam_pairs, |p| p.score, rng);

    nam_details1.time_sort_nams = sort_nams(&mut nams1, rng);
    nam_details2.time_sort_nams = sort_nams(&mut nams2, rng);

    let mut records = vec![];
    match get_best_paired_mapping_location(&nam_pairs, &nams1, &nams2, insert_size_distribution) {
        MappedNams::Individual(best1, best2) => {
            if let Some(nam1) = best1 {
                records.push(paf_record_from_nam(
                    nam1,
                    &r1.name,
                    references,
                    r1.sequence.len(),
                    None,
                    End::One,
                ));
            }
            if let Some(nam2) = best2 {
                records.push(paf_record_from_nam(
                    nam2,
                    &r2.name,
                    references,
                    r2.sequence.len(),
                    None,
                    End::Two,
                ));
            }
        }
        MappedNams::Pair(nam1, nam2, _) => {
            records.push(paf_record_from_nam(
                nam1,
                &r1.name,
                references,
                r1.sequence.len(),
                None,
                End::One,
            ));
            records.push(paf_record_from_nam(
                nam2,
                &r2.name,
                references,
                r2.sequence.len(),
                None,
                End::Two,
            ));
        }
    }

    nam_details1 += nam_details2;
    (records, nam_details1.into())
}

/// Map a paired-end read pair to the reference and estimate abundances
///
/// This implements abundance estimation mode (`--aemb`)
pub fn abundances_paired_end_read(
    r1: &SequenceRecord,
    r2: &SequenceRecord,
    index: &StrobemerIndex,
    abundances: &mut [f64],
    rescue_distance: usize,
    insert_size_distribution: &mut InsertSizeDistribution,
    mcs_strategy: McsStrategy,
    chainer: &Chainer,
    rng: &mut Rng,
) {
    let (nam_details1, mut nams1) =
        get_nams_by_chaining(&r1.sequence, index, chainer, rescue_distance, mcs_strategy);
    let (nam_details2, mut nams2) =
        get_nams_by_chaining(&r2.sequence, index, chainer, rescue_distance, mcs_strategy);

    if nams1.is_empty() && nams2.is_empty() {
        return;
    }

    let mut nam_pairs = get_nam_pairs(
        &mut nams1,
        &mut nams2,
        insert_size_distribution.mu,
        insert_size_distribution.sigma,
        &nam_details1,
        &nam_details2,
    );
    shuffle_best(&mut nam_pairs, |p| p.score, rng);

    sort_nams(&mut nams1, rng);
    sort_nams(&mut nams2, rng);

    match get_best_paired_mapping_location(&nam_pairs, &nams1, &nams2, insert_size_distribution) {
        MappedNams::Individual(_, _) => {
            for (nams, read_len) in [(&nams1, r1.sequence.len()), (&nams2, r2.sequence.len())] {
                let n_best = nams
                    .iter()
                    .take_while(|nam| nam.score == nams[0].score)
                    .count();
                let weight = read_len as f64 / n_best as f64;
                for nam in &nams[0..n_best] {
                    abundances[nam.ref_id] += weight;
                }
            }
        }
        MappedNams::Pair(_, _, joint_score) => {
            let n_best = nam_pairs
                .iter()
                .take_while(|nam_pair| nam_pair.score == joint_score)
                .count();
            let weight_r1 = r1.sequence.len() as f64 / n_best as f64;
            let weight_r2 = r2.sequence.len() as f64 / n_best as f64;
            for NamPair {
                nam1,
                nam2,
                score: _,
            } in &nam_pairs[..n_best]
            {
                abundances[nam1.ref_id] += weight_r1;
                abundances[nam2.ref_id] += weight_r2;
            }
        }
    }
}

enum MappedNams<'a> {
    /// Two independent best NAMs (one per read)
    Individual(Option<&'a Nam>, Option<&'a Nam>),
    /// A proper paired NAMs (nam1, nam2, pairing score)
    Pair(&'a Nam, &'a Nam, f64),
}

/// Choose between:
/// - the best individual NAMs
/// - the best proper pair of NAMs
///
/// Also updates the insert size distribution using confident pairs.
///
/// For paired-end mapping and abundance estimation modes only
fn get_best_paired_mapping_location<'a>(
    nam_pairs: &'a [NamPair],
    nams1: &'a [Nam],
    nams2: &'a [Nam],
    insert_size_distribution: &mut InsertSizeDistribution,
) -> MappedNams<'a> {
    let best_nam1 = nams1.first();
    let best_nam2 = nams2.first();

    // Score if reads are treated independently.
    let individual_score = best_nam1.map_or(0.0, |nam| nam.score as f64)
        + best_nam2.map_or(0.0, |nam| nam.score as f64);

    // Prefer a proper pair only if it beats a penalized individual mapping.
    // Divisor 2 is penalty for being mapped individually
    if let Some(NamPair { nam1, nam2, score }) = nam_pairs.first()
        && *score >= individual_score / 2.0
    {
        // Update insert size using confident proper pairs.
        if insert_size_distribution.sample_size < 400 {
            insert_size_distribution.update(nam1.ref_start.abs_diff(nam2.ref_start));
        }

        MappedNams::Pair(nam1, nam2, *score)
    } else {
        MappedNams::Individual(best_nam1, best_nam2)
    }
}

/// Properly paired nams
#[derive(Debug)]
pub struct NamPair {
    pub nam1: Nam,
    pub nam2: Nam,
    pub score: f64,
}

/// Build all plausible forward/revcomp mapping pairings
pub fn get_nam_pairs(
    nams1: &mut [Nam],
    nams2: &mut [Nam],
    mu: f32,
    sigma: f32,
    details1: &NamDetails,
    details2: &NamDetails,
) -> Vec<NamPair> {
    let mut nam_pairs = vec![];
    if nams1.is_empty() || nams2.is_empty() {
        return nam_pairs;
    }

    let (fwd1, rev1): (&mut [Nam], &mut [Nam]) =
        split_nams_by_orientation_checked(nams1, details1.both_orientations);
    let (fwd2, rev2): (&mut [Nam], &mut [Nam]) =
        split_nams_by_orientation_checked(nams2, details2.both_orientations);

    if !fwd1.is_empty() && !rev2.is_empty() {
        fwd1.sort_unstable_by_key(|nam| (nam.ref_id, nam.projected_ref_start()));
        rev2.sort_unstable_by_key(|nam| (nam.ref_id, nam.projected_ref_start()));
        nam_pairs.extend(find_pairs(fwd1, rev2, mu, sigma, false));
    }

    if !fwd2.is_empty() && !rev1.is_empty() {
        fwd2.sort_unstable_by_key(|nam| (nam.ref_id, nam.projected_ref_start()));
        rev1.sort_unstable_by_key(|nam| (nam.ref_id, nam.projected_ref_start()));
        nam_pairs.extend(find_pairs(fwd2, rev1, mu, sigma, true));
    }

    nam_pairs.sort_unstable_by(|a, b| b.score.total_cmp(&a.score));
    nam_pairs
}

/// Split nams into (forward, revcomp),
/// if only 1 orientation exists, returns it and a empty slice
fn split_nams_by_orientation_checked(nams: &mut [Nam], both: bool) -> (&mut [Nam], &mut [Nam]) {
    if both {
        split_nams_by_orientation(nams)
    } else if nams[0].is_revcomp {
        (&mut [], nams)
    } else {
        (nams, &mut [])
    }
}

/// In-place partition of NAMs by orientation:
/// forward on the left, revcomp on the right.
/// Returns two slices separating (forward, revcomp)
fn split_nams_by_orientation(nams: &mut [Nam]) -> (&mut [Nam], &mut [Nam]) {
    let mut left = 0;
    let mut right = nams.len();

    while left < right {
        if nams[left].is_revcomp {
            right -= 1;
            nams.swap(left, right);
        } else {
            left += 1;
        }
    }

    nams.split_at_mut(left)
}

/// Find most forward/revcomp pairs using a two-pointer scan.
/// Assumes both slices are sorted by (ref_id, projected_ref_start).
fn find_pairs(fwd: &[Nam], rev: &[Nam], mu: f32, sigma: f32, swap_order: bool) -> Vec<NamPair> {
    let mut out = Vec::new();
    let max_dist = (mu + 10.0 * sigma).ceil() as usize; // distance cutoff from insert size distribution
    let mut rev_ptr = 0;
    let mut last_paired = None;

    for f in fwd {
        // Advance revcomp pointer to the first possible candidate
        while rev_ptr < rev.len()
            && (rev[rev_ptr].ref_id < f.ref_id
                || rev[rev_ptr].ref_id == f.ref_id
                    && rev[rev_ptr].projected_ref_start() < f.projected_ref_start())
        {
            rev_ptr += 1;
        }
        if rev_ptr == rev.len() {
            break;
        }
        if rev[rev_ptr].ref_id > f.ref_id {
            continue;
        }

        // Scan window of revcomp nams within distance limit.
        let mut best = None;
        let mut i = rev_ptr;
        while i < rev.len()
            && rev[i].ref_id == f.ref_id
            && (rev[i].projected_ref_start() - f.projected_ref_start()) <= max_dist
        {
            let r = &rev[i];
            // The pairing score gets a bonus based on the reference distance of the two chosen nam
            // paired and from our current knowledge of the reference distance distribution
            let x = f.ref_start.abs_diff(r.ref_start);
            let score = f.score as f64
                + r.score as f64
                + 0.001f64.max((normal_pdf(x as f32, mu, sigma) + 1.0).ln() as f64);

            if best.is_none_or(|(_, highest_score)| score > highest_score) {
                best = Some((i, score));
            }
            i += 1;
        }

        // Highest scoring candidate
        let Some((best_id, score)) = best else {
            continue;
        };
        let r = &rev[best_id];

        // If the same revcomp nam was paired previously, keep only the better scoring pair.
        if let Some((last_id, prev_score)) = last_paired
            && last_id == best_id
        {
            if score <= prev_score {
                continue; // keep the previous better pair
            }
            out.pop(); // replace it
        }

        out.push(if swap_order {
            NamPair {
                nam1: r.clone(),
                nam2: f.clone(),
                score,
            }
        } else {
            NamPair {
                nam1: f.clone(),
                nam2: r.clone(),
                score,
            }
        });

        last_paired = Some((best_id, score));
        rev_ptr = best_id;
    }

    out
}

/// A scored alignment pair
#[derive(Debug, Clone)]
pub struct PairedAlignments {
    pub score: f64,
    pub alignment1: Alignment,
    pub alignment2: Alignment,
}

/// Align both reads and collect all plausible paired and individual alignments.
/// First tries NAMs that could form proper pairs, then falls back to
/// unpaired chains with mate rescue if the max_tries is not yet exhausted.
pub fn get_paired_alignment(
    aligner: &Aligner,
    mut nams1: Vec<Nam>,
    mut nams2: Vec<Nam>,
    max_tries: usize,
    dropoff: f32,
    references: &[RefSequence],
    read1: &Read,
    read2: &Read,
    details1: &mut Details,
    details2: &mut Details,
    mu: f32,
    sigma: f32,
    k: usize,
    rng: &mut Rng,
) -> (Vec<PairedAlignments>, Vec<Alignment>, Vec<Alignment>) {
    let mut nam_pairs = get_nam_pairs(
        &mut nams1,
        &mut nams2,
        mu,
        sigma,
        &details1.nam,
        &details2.nam,
    );
    let max_score = nam_pairs.first().map_or(0.0, |p| p.score as f32);

    let mut paired_alignments = vec![];
    let mut alignments1 = vec![];
    let mut alignments2 = vec![];

    for p in nam_pairs.iter_mut() {
        if p.score as f32 / max_score < dropoff || paired_alignments.len() == max_tries {
            break;
        }
        let a1 = align_nam(aligner, &mut p.nam1, references, read1, k, details1);
        let a2 = align_nam(aligner, &mut p.nam2, references, read2, k, details2);
        if let (Some(a1), Some(a2)) = (&a1, &a2) {
            paired_alignments.push(PairedAlignments {
                score: compute_combined_score(a1, a2, mu, sigma),
                alignment1: a1.clone(),
                alignment2: a2.clone(),
            });
        }
        alignments1.extend(a1);
        alignments2.extend(a2);
    }

    // Fallback to unpaired NAMs only if we didn't have enough alignment pairs
    if paired_alignments.len() < max_tries {
        let mut unpaired_nams = get_unpaired_nams(nams1, nams2, &nam_pairs);
        let max_score = max_score.max(unpaired_nams.first().map_or(0.0, |u| u.nam.score));
        for u in unpaired_nams.iter_mut() {
            if u.nam.score / max_score < dropoff || paired_alignments.len() == max_tries {
                break;
            }
            let (a1, a2) = if u.read1 {
                (
                    align_nam(aligner, &mut u.nam, references, read1, k, details1),
                    rescue_align(aligner, &u.nam, references, read2, mu, sigma, k, details2),
                )
            } else {
                (
                    rescue_align(aligner, &u.nam, references, read1, mu, sigma, k, details1),
                    align_nam(aligner, &mut u.nam, references, read2, k, details2),
                )
            };
            if let (Some(a1), Some(a2)) = (&a1, &a2) {
                paired_alignments.push(PairedAlignments {
                    score: compute_combined_score(a1, a2, mu, sigma),
                    alignment1: a1.clone(),
                    alignment2: a2.clone(),
                });
            }
            alignments1.extend(a1);
            alignments2.extend(a2);
        }
    }
    alignments1.sort_unstable_by_key(|a| Reverse(a.score));
    alignments2.sort_unstable_by_key(|a| Reverse(a.score));
    shuffle_best(&mut alignments1, |a| a.score, rng);
    shuffle_best(&mut alignments2, |a| a.score, rng);

    if let Some(a1) = alignments1.first()
        && let Some(a2) = alignments2.first()
    {
        paired_alignments.push(PairedAlignments {
            score: compute_combined_score(a1, a2, mu, sigma),
            alignment1: a1.clone(),
            alignment2: a2.clone(),
        });
    }

    paired_alignments.sort_unstable_by(|a, b| b.score.total_cmp(&a.score));
    deduplicate_scored_pairs(&mut paired_alignments);
    shuffle_best(&mut paired_alignments, |a| a.score, rng);

    (paired_alignments, alignments1, alignments2)
}

/// Compute the combined score for a read pair.
/// Properly oriented pairs within the expected insert size get a log-normal
/// distance bonus
fn compute_combined_score(a1: &Alignment, a2: &Alignment, mu: f32, sigma: f32) -> f64 {
    let r1_r2 = a2.is_revcomp
        && !a1.is_revcomp
        && a1.ref_start <= a2.ref_start
        && (a2.ref_start - a1.ref_start) < (mu + 10.0 * sigma) as usize;
    let r2_r1 = a1.is_revcomp
        && !a2.is_revcomp
        && a2.ref_start <= a1.ref_start
        && (a1.ref_start - a2.ref_start) < (mu + 10.0 * sigma) as usize;

    if r1_r2 || r2_r1 {
        let x = a1.ref_start.abs_diff(a2.ref_start);
        a1.score as f64
            + a2.score as f64
            + (-20.0f64 + 0.001).max(normal_pdf(x as f32, mu, sigma).ln() as f64)
    } else {
        a1.score as f64 + a2.score as f64 - 20.0
    }
}

/// Remove consecutive identical alignment pairs and leave only the first.
fn deduplicate_scored_pairs(pairs: &mut Vec<PairedAlignments>) {
    pairs.dedup_by(|a, b| {
        a.alignment1.ref_start == b.alignment1.ref_start
            && a.alignment1.reference_id == b.alignment1.reference_id
            && a.alignment2.ref_start == b.alignment2.ref_start
            && a.alignment2.reference_id == b.alignment2.reference_id
    });
}

/// Align a NAM against the reference,
fn align_nam(
    aligner: &Aligner,
    nam: &mut Nam,
    references: &[RefSequence],
    read: &Read,
    k: usize,
    details: &mut Details,
) -> Option<Alignment> {
    let consistent = reverse_nam_if_needed(nam, read, references, k);
    details.inconsistent_nams += !consistent as usize;
    // Can do piecewise here
    let aln = extend_seed(aligner, nam, references, read, consistent, true);
    details.tried_alignment += 1;
    details.gapped += aln.as_ref().map_or(0, |a| a.gapped as usize);
    aln
}

/// Nam without any pairing partner
#[derive(Debug)]
struct UnpairedNam {
    nam: Nam,
    read1: bool,
}

/// Returns the nams that did not get paired
fn get_unpaired_nams(nams1: Vec<Nam>, nams2: Vec<Nam>, nam_pairs: &[NamPair]) -> Vec<UnpairedNam> {
    let mut paired1 = vec![false; nams1.len()];
    let mut paired2 = vec![false; nams2.len()];

    for pair in nam_pairs {
        paired1[pair.nam1.nam_id] = true;
        paired2[pair.nam2.nam_id] = true;
    }

    let mut unpaired_nams =
        Vec::with_capacity((nams1.len() + nams2.len()).saturating_sub(nam_pairs.len()));

    for nam in nams1 {
        if !paired1[nam.nam_id] {
            unpaired_nams.push(UnpairedNam { nam, read1: true });
        }
    }
    for nam in nams2 {
        if !paired2[nam.nam_id] {
            unpaired_nams.push(UnpairedNam { nam, read1: false });
        }
    }

    unpaired_nams.sort_unstable_by(|a, b| b.nam.score.total_cmp(&a.nam.score));
    unpaired_nams
}

pub fn make_alignment(
    reference_id: usize,
    ref_start: usize,
    score: u32,
    is_revcomp: bool,
) -> Alignment {
    Alignment {
        reference_id,
        ref_start,
        score,
        is_revcomp,
        edit_distance: 0,
        soft_clip_left: 0,
        soft_clip_right: 0,
        length: 50,
        cigar: Default::default(),
        gapped: false,
        rescued: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{chainer::Anchor, nam::Nam};

    fn make_nam(
        nam_id: usize,
        ref_id: usize,
        ref_start: usize,
        ref_end: usize,
        score: f32,
        is_revcomp: bool,
    ) -> Nam {
        Nam {
            nam_id,
            ref_id,
            ref_start,
            ref_end,
            query_start: 0,
            query_end: ref_end - ref_start,
            anchors: vec![Anchor {
                ref_id: 0,
                ref_start: 0,
                query_start: 0,
            }],
            matching_bases: ref_end - ref_start,
            score,
            is_revcomp,
        }
    }

    #[test]
    fn find_pairs_swap_order() {
        let mu = 300.0_f32;
        let sigma = 50.0_f32;
        let fwd = vec![make_nam(1, 0, 100, 150, 20.0, false)];
        let rev = vec![make_nam(2, 0, 300, 350, 20.0, true)];
        let pairs = find_pairs(&fwd, &rev, mu, sigma, true);
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].nam1.nam_id, 2);
        assert_eq!(pairs[0].nam2.nam_id, 1);
    }

    #[test]
    fn find_pairs_selects_best_scoring_candidate() {
        let mu = 300.0_f32;
        let sigma = 50.0_f32;
        let fwd = vec![make_nam(0, 0, 100, 150, 20.0, false)];
        let rev = vec![
            make_nam(1, 0, 300, 350, 10.0, true),
            make_nam(2, 0, 400, 450, 30.0, true),
        ];
        let pairs = find_pairs(&fwd, &rev, mu, sigma, false);
        assert!(!pairs.is_empty());
        assert_eq!(pairs[0].nam2.nam_id, 2);
    }

    #[test]
    fn find_pairs_score() {
        let mu = 300.0_f32;
        let sigma = 50.0_f32;
        let fwd = vec![make_nam(0, 0, 100, 150, 20.0, false)];
        let rev = vec![make_nam(1, 0, 350, 400, 15.0, true)];
        let pairs = find_pairs(&fwd, &rev, mu, sigma, false);
        assert_eq!(pairs.len(), 1);
        assert!(pairs[0].score > 20.0 + 15.0);
    }

    #[test]
    fn get_nam_pairs_empty() {
        let mut nams1: Vec<Nam> = vec![];
        let mut nams2: Vec<Nam> = vec![];
        let d1 = NamDetails {
            both_orientations: false,
            ..Default::default()
        };
        let d2 = NamDetails {
            both_orientations: false,
            ..Default::default()
        };
        let pairs = get_nam_pairs(&mut nams1, &mut nams2, 300.0, 50.0, &d1, &d2);
        assert!(pairs.is_empty());
    }

    #[test]
    fn get_nam_pairs_fwd_orientation() {
        let mut nams1 = vec![make_nam(0, 0, 100, 150, 20.0, false)];
        let mut nams2 = vec![make_nam(1, 0, 350, 400, 20.0, true)];
        let d1 = NamDetails {
            both_orientations: false,
            ..Default::default()
        };
        let d2 = NamDetails {
            both_orientations: false,
            ..Default::default()
        };
        let pairs = get_nam_pairs(&mut nams1, &mut nams2, 300.0, 50.0, &d1, &d2);
        assert!(!pairs.is_empty());
        for w in pairs.windows(2) {
            assert!(w[0].score >= w[1].score);
        }
    }

    #[test]
    fn get_nam_pairs_rev_orientation() {
        let mut nams1 = vec![make_nam(0, 0, 350, 400, 20.0, true)];
        let mut nams2 = vec![make_nam(1, 0, 100, 150, 20.0, false)];
        let d1 = NamDetails {
            both_orientations: false,
            ..Default::default()
        };
        let d2 = NamDetails {
            both_orientations: false,
            ..Default::default()
        };
        let pairs = get_nam_pairs(&mut nams1, &mut nams2, 300.0, 50.0, &d1, &d2);
        assert!(!pairs.is_empty());
    }

    #[test]
    fn get_unpaired_nam() {
        let nam_a = make_nam(0, 0, 100, 150, 10.0, false);
        let nam_b = make_nam(1, 0, 200, 250, 5.0, false);
        let nam_c = make_nam(0, 0, 300, 350, 8.0, true);
        let nam_d = make_nam(1, 0, 400, 450, 20.0, true);
        let nams1 = vec![nam_a, nam_b];
        let nams2 = vec![nam_c, nam_d];
        let paired = vec![NamPair {
            nam1: make_nam(0, 0, 100, 150, 10.0, false),
            nam2: make_nam(1, 0, 400, 450, 20.0, true),
            score: 100.0,
        }];
        let unpaired = get_unpaired_nams(nams1, nams2, &paired);
        assert_eq!(unpaired.len(), 2);
        assert_eq!(unpaired[0].nam.nam_id, 0);
        assert_eq!(unpaired[1].nam.nam_id, 1);
    }

    #[test]
    fn get_unpaired_nams_all_paired() {
        let nams1 = vec![make_nam(0, 0, 100, 150, 10.0, false)];
        let nams2 = vec![make_nam(0, 0, 350, 400, 10.0, true)];
        let paired = vec![NamPair {
            nam1: make_nam(0, 0, 100, 150, 10.0, false),
            nam2: make_nam(0, 0, 350, 400, 10.0, true),
            score: 50.0,
        }];
        let unpaired = get_unpaired_nams(nams1, nams2, &paired);
        assert!(unpaired.is_empty());
    }

    #[test]
    fn get_unpaired_nams_none_paired() {
        let nams1 = vec![make_nam(0, 0, 100, 150, 10.0, false)];
        let nams2 = vec![make_nam(0, 0, 350, 400, 10.0, true)];
        let unpaired = get_unpaired_nams(nams1, nams2, &[]);
        assert_eq!(unpaired.len(), 2);
    }

    #[test]
    fn get_unpaired_nams_flag_correct() {
        let nams1 = vec![make_nam(0, 0, 100, 150, 10.0, false)];
        let nams2 = vec![make_nam(0, 0, 200, 250, 10.0, true)];
        let unpaired = get_unpaired_nams(nams1, nams2, &[]);
        let u0 = unpaired
            .iter()
            .find(|u| u.nam.nam_id == 0 && u.read1)
            .unwrap();
        let u1 = unpaired
            .iter()
            .find(|u| u.nam.nam_id == 0 && !u.read1)
            .unwrap();
        assert!(u0.read1);
        assert!(!u1.read1);
    }

    #[test]
    fn compute_combined_score_bonus() {
        let mu = 300.0_f32;
        let sigma = 50.0_f32;
        let a1 = make_alignment(0, 100, 40, false);
        let a2 = make_alignment(0, 380, 40, true);
        let score = compute_combined_score(&a1, &a2, mu, sigma);
        assert!(score > (a1.score + a2.score) as f64 - 20.0);
    }

    #[test]
    fn compute_combined_score_penalty() {
        let mu = 300.0_f32;
        let sigma = 50.0_f32;
        let a1 = make_alignment(0, 100, 40, false);
        let a2 = make_alignment(0, 380, 40, false);
        let score = compute_combined_score(&a1, &a2, mu, sigma);
        assert_eq!(score, (a1.score + a2.score) as f64 - 20.0);
    }

    #[test]
    fn deduplicate_removes_consecutive_pairs() {
        let aln = make_alignment(0, 100, 40, false);
        let mut pairs = vec![
            PairedAlignments {
                score: 100.0,
                alignment1: aln.clone(),
                alignment2: aln.clone(),
            },
            PairedAlignments {
                score: 90.0,
                alignment1: aln.clone(),
                alignment2: aln.clone(),
            },
        ];
        deduplicate_scored_pairs(&mut pairs);
        assert_eq!(pairs.len(), 1);
    }

    #[test]
    fn deduplicate_keeps_distinct_pairs() {
        let aln_a = make_alignment(0, 100, 40, false);
        let aln_b = make_alignment(0, 500, 40, true);
        let mut pairs = vec![
            PairedAlignments {
                score: 100.0,
                alignment1: aln_a.clone(),
                alignment2: aln_a.clone(),
            },
            PairedAlignments {
                score: 80.0,
                alignment1: aln_b.clone(),
                alignment2: aln_b.clone(),
            },
        ];
        deduplicate_scored_pairs(&mut pairs);
        assert_eq!(pairs.len(), 2);
    }

    #[test]
    fn deduplicate_single_pair_unchanged() {
        let aln = make_alignment(0, 100, 40, false);
        let mut pairs = vec![PairedAlignments {
            score: 100.0,
            alignment1: aln.clone(),
            alignment2: aln.clone(),
        }];
        deduplicate_scored_pairs(&mut pairs);
        assert_eq!(pairs.len(), 1);
    }

    #[test]
    fn split_orientation() {
        let mut nams = vec![
            make_nam(0, 0, 100, 150, 10.0, false),
            make_nam(1, 0, 200, 250, 10.0, true),
            make_nam(2, 0, 300, 350, 10.0, false),
            make_nam(3, 0, 400, 450, 10.0, true),
        ];
        let (fwd, rev) = split_nams_by_orientation(&mut nams);
        assert_eq!(fwd.len(), 2);
        assert_eq!(rev.len(), 2);
        assert!(fwd.iter().all(|n| !n.is_revcomp));
        assert!(rev.iter().all(|n| n.is_revcomp));
    }

    #[test]
    fn split_checked_all_fwd() {
        let mut nams = vec![
            make_nam(0, 0, 100, 150, 10.0, false),
            make_nam(1, 0, 200, 250, 10.0, false),
        ];
        let (fwd, rev) = split_nams_by_orientation_checked(&mut nams, false);
        assert_eq!(fwd.len(), 2);
        assert_eq!(rev.len(), 0);
    }

    #[test]
    fn split_checked_all_rev() {
        let mut nams = vec![
            make_nam(0, 0, 100, 150, 10.0, true),
            make_nam(1, 0, 200, 250, 10.0, true),
        ];
        let (fwd, rev) = split_nams_by_orientation_checked(&mut nams, false);
        assert_eq!(fwd.len(), 0);
        assert_eq!(rev.len(), 2);
    }

    #[test]
    fn find_pairs_within_distance() {
        let mu = 300.0_f32;
        let sigma = 50.0_f32;
        let fwd = vec![make_nam(0, 0, 100, 150, 20.0, false)];
        let rev = vec![make_nam(1, 0, 350, 400, 20.0, true)];
        let pairs = find_pairs(&fwd, &rev, mu, sigma, false);
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].nam1.nam_id, 0);
        assert_eq!(pairs[0].nam2.nam_id, 1);
    }

    #[test]
    fn find_pairs_beyond_distance() {
        let mu = 300.0_f32;
        let sigma = 50.0_f32;
        let fwd = vec![make_nam(0, 0, 100, 150, 20.0, false)];
        let rev = vec![make_nam(1, 0, 10_100, 10_150, 20.0, true)];
        let pairs = find_pairs(&fwd, &rev, mu, sigma, false);
        assert!(pairs.is_empty());
    }
}
