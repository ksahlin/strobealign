use crate::chainer::Chainer;
use crate::details::Details;
use crate::index::StrobemerIndex;
use crate::insertsize::InsertSizeDistribution;
use crate::io::fasta::RefSequence;
use crate::io::paf::PafRecord;
use crate::io::record::{End, SequenceRecord};
use crate::mapper::mapping_quality;
use crate::mcsstrategy::McsStrategy;
use crate::nam::{Nam, get_nams_by_chaining, sort_nams};
use crate::pairing::{PairedNams, get_paired_nams};
use crate::shuffle::shuffle_best;
use bumpalo::Bump;
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
    arena: &Bump,
) -> (Vec<PafRecord>, Details) {
    let (mut nam_details, mut nams) = get_nams_by_chaining(
        &record.sequence,
        index,
        chainer,
        rescue_distance,
        mcs_strategy,
        arena,
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
    arena: &Bump,
) {
    let (_, mut nams) = get_nams_by_chaining(
        &record.sequence,
        index,
        chainer,
        rescue_distance,
        mcs_strategy,
        arena,
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
    arena: &Bump,
) -> (Vec<PafRecord>, Details) {
    let (mut nam_details1, mut nams1) = get_nams_by_chaining(
        &r1.sequence,
        index,
        chainer,
        rescue_distance,
        mcs_strategy,
        arena,
    );
    let (mut nam_details2, mut nams2) = get_nams_by_chaining(
        &r2.sequence,
        index,
        chainer,
        rescue_distance,
        mcs_strategy,
        arena,
    );

    if nams1.is_empty() && nams2.is_empty() {
        nam_details1 += nam_details2;
        return (vec![], nam_details1.into());
    }

    let mut paired_nams = get_paired_nams(
        &mut nams1,
        &mut nams2,
        insert_size_distribution.mu,
        insert_size_distribution.sigma,
        &nam_details1,
        &nam_details2,
    );
    shuffle_best(&mut paired_nams, |p| p.score, rng);

    nam_details1.time_sort_nams = sort_nams(&mut nams1, rng);
    nam_details2.time_sort_nams = sort_nams(&mut nams2, rng);

    let mut records = vec![];
    match get_best_paired_mapping_location(&paired_nams, &nams1, &nams2, insert_size_distribution) {
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
    arena: &Bump,
) {
    let (nam_details1, mut nams1) = get_nams_by_chaining(
        &r1.sequence,
        index,
        chainer,
        rescue_distance,
        mcs_strategy,
        arena,
    );
    let (nam_details2, mut nams2) = get_nams_by_chaining(
        &r2.sequence,
        index,
        chainer,
        rescue_distance,
        mcs_strategy,
        arena,
    );

    if nams1.is_empty() && nams2.is_empty() {
        return;
    }

    let mut paired_nams = get_paired_nams(
        &mut nams1,
        &mut nams2,
        insert_size_distribution.mu,
        insert_size_distribution.sigma,
        &nam_details1,
        &nam_details2,
    );
    shuffle_best(&mut paired_nams, |p| p.score, rng);

    sort_nams(&mut nams1, rng);
    sort_nams(&mut nams2, rng);

    match get_best_paired_mapping_location(&paired_nams, &nams1, &nams2, insert_size_distribution) {
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
            let n_best = paired_nams
                .iter()
                .take_while(|nam_pair| nam_pair.score == joint_score)
                .count();
            let weight_r1 = r1.sequence.len() as f64 / n_best as f64;
            let weight_r2 = r2.sequence.len() as f64 / n_best as f64;
            for PairedNams {
                nam1,
                nam2,
                score: _,
            } in &paired_nams[..n_best]
            {
                abundances[nam1.ref_id] += weight_r1;
                abundances[nam2.ref_id] += weight_r2;
            }
        }
    }
}

enum MappedNams<'a> {
    /// Two independent best NAMs (one per read)
    Individual(Option<&'a Nam<'a>>, Option<&'a Nam<'a>>),
    /// A proper paired NAMs (nam1, nam2, pairing score)
    Pair(&'a Nam<'a>, &'a Nam<'a>, f64),
}

/// Choose between:
/// - the best individual NAMs
/// - the best proper pair of NAMs
///
/// Also updates the insert size distribution using confident pairs.
///
/// For paired-end mapping and abundance estimation modes only
fn get_best_paired_mapping_location<'a>(
    paired_nams: &'a [PairedNams],
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
    if let Some(PairedNams { nam1, nam2, score }) = paired_nams.first()
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
