use fastrand::Rng;

use crate::chainer::Chainer;
use crate::details::{Details, NamDetails};
use crate::fasta::RefSequence;
use crate::fastq::{End, SequenceRecord};
use crate::index::StrobemerIndex;
use crate::insertsize::InsertSizeDistribution;
use crate::mapper::mapping_quality;
use crate::math::normal_pdf;
use crate::mcsstrategy::McsStrategy;
use crate::nam::{Nam, get_nams_by_chaining, sort_nams};
use crate::paf::PafRecord;

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
    sort_nams(&mut nams, &mut nam_details, rng);

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
    let (mut nam_details, mut nams) = get_nams_by_chaining(
        &record.sequence,
        index,
        chainer,
        rescue_distance,
        mcs_strategy,
    );
    sort_nams(&mut nams, &mut nam_details, rng);

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

    let nam_pairs = get_nam_pairs(
        &mut nams1,
        &mut nams2,
        insert_size_distribution.mu,
        insert_size_distribution.sigma,
        &nam_details1,
        &nam_details2,
    );

    sort_nams(&mut nams1, &mut nam_details1, rng);
    sort_nams(&mut nams2, &mut nam_details2, rng);

    let mut records = vec![];
    match get_best_paired_mapping_location(&nam_pairs, &nams1, &nams2, insert_size_distribution) {
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
    let (mut nam_details1, mut nams1) =
        get_nams_by_chaining(&r1.sequence, index, chainer, rescue_distance, mcs_strategy);
    let (mut nam_details2, mut nams2) =
        get_nams_by_chaining(&r2.sequence, index, chainer, rescue_distance, mcs_strategy);

    if nams1.is_empty() && nams2.is_empty() {
        return;
    }

    let nam_pairs = get_nam_pairs(
        &mut nams1,
        &mut nams2,
        insert_size_distribution.mu,
        insert_size_distribution.sigma,
        &nam_details1,
        &nam_details2,
    );

    sort_nams(&mut nams1, &mut nam_details1, rng);
    sort_nams(&mut nams2, &mut nam_details2, rng);

    match get_best_paired_mapping_location(&nam_pairs, &nams1, &nams2, insert_size_distribution) {
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
    }
}

/// Represents either:
/// A proper paired mapping (nam1, nam2, pairing score)
/// Two independent best mappings (one per read)
enum MappedNams<'a> {
    /// A proper paired mapping (nam1, nam2, pairing score)
    Pair(&'a Nam, &'a Nam, f64),
    /// Two independent best mappings (one per read)
    Individual(Option<&'a Nam>, Option<&'a Nam>),
}

/// Choose between:
/// - the best proper pair of mappings
/// - the best individual mappings
///
/// Also updates the insert size distribution using confident pairs.
///
/// For paired mapping and abundance only
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

#[derive(Debug)]
pub struct NamPair {
    pub nam1: Nam,
    pub nam2: Nam,
    pub score: f64,
}

/// Build all plausible forward/revcomp mapping pairings
fn get_nam_pairs(
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

/// Find most forward/revcomp pairs using a linear two-pointer scan.
/// Assumes both slices are sorted by (ref_id, projected_ref_start).
fn find_pairs(fwd: &[Nam], rev: &[Nam], mu: f32, sigma: f32, swap_order: bool) -> Vec<NamPair> {
    let mut out = Vec::new();
    let max_dist = mu + 10.0 * sigma; // distance cutoff from insert size distribution
    let mut rev_ptr = 0;
    let mut last_paired: Option<(usize, f64)> = None;

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
        let mut best: Option<(usize, f32)> = None;
        let mut i = rev_ptr;
        while i < rev.len()
            && rev[i].ref_id == f.ref_id
            && (rev[i].projected_ref_start() - f.projected_ref_start()) as f32 <= max_dist
        {
            let s = rev[i].score;
            if best.is_none_or(|(_, bs)| s > bs) {
                best = Some((i, s));
            }
            i += 1;
        }

        // Best scoring candidate
        let Some((best_id, best_score)) = best else {
            continue;
        };
        let r = &rev[best_id];

        // The pairing score get's a bonus based on the reference distance of the two chosen nam
        // paired and from our current knowledge of the reference distance distribution
        let score = f.score as f64
            + best_score as f64
            + (-20.0f64 + 0.001)
                .max(normal_pdf(f.ref_start.abs_diff(r.ref_start) as f32, mu, sigma).ln() as f64);

        // If the same recomp nam was paired previously, keep only the better scoring pair.
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
