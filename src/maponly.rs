use fastrand::Rng;
use crate::fasta::RefSequence;
use crate::fastq::SequenceRecord;
use crate::index::StrobemerIndex;
use crate::insertsize::InsertSizeDistribution;
use crate::mapper::get_best_scoring_nam_pairs;
use crate::nam::{get_nams, Nam};
use crate::paf::PafRecord;

/// Map a single-end read to the reference and return PAF records
///
/// This implements "mapping-only" mode in which no base-level alignments are computed
pub fn map_single_end_read(
    record: &SequenceRecord,
    index: &StrobemerIndex,
    references: &[RefSequence],
    rescue_level: usize,
    rng: &mut Rng,
) -> Vec<PafRecord> {
    let (_, nams) = get_nams(&record.sequence, index, rescue_level, rng);

    if nams.is_empty() {
        vec![]
    } else {
        vec![paf_record_from_nam(&nams[0], &record.name, &references, record.sequence.len())]
    }
}

/// Convert Nam into PAF record
fn paf_record_from_nam(nam: &Nam, name: &str, references: &[RefSequence], query_length: usize) -> PafRecord {
    PafRecord {
        query_name: name.into(),
        query_length: query_length as u64,
        query_start: nam.query_start as u64,
        query_end: nam.query_end as u64,
        is_revcomp: nam.is_revcomp,
        target_name: references[nam.ref_id].name.clone(),
        target_length: references[nam.ref_id].sequence.len() as u64,
        target_start: nam.ref_start as u64,
        target_end: nam.ref_end as u64,
        n_matches: nam.n_hits as u64,
        alignment_length: (nam.ref_end - nam.ref_start) as u64,
        mapping_quality: None,
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
    rescue_level: usize,
    insert_size_distribution: &mut InsertSizeDistribution,
    rng: &mut Rng,
) -> Vec<PafRecord> {
    let nams_pair = [
        get_nams(&r1.sequence, index, rescue_level, rng).1,
        get_nams(&r2.sequence, index, rescue_level, rng).1,
    ];

    let best_nam_pair = get_best_paired_map_location(
        &nams_pair[0],
        &nams_pair[1],
        insert_size_distribution,
        // TODO aemb:
        //r1.sequence.len(), r2.sequence.len(),
        //abundances,
        //map_param.output_format == OutputFormat::Abundance
    );

    let mut records = vec![];
    if let Some(nam) = best_nam_pair.0 {
        records.push(paf_record_from_nam(&nam, &r1.name, &references, r1.sequence.len()))
    }
    if let Some(nam) = best_nam_pair.1 {
        records.push(paf_record_from_nam(&nam, &r2.name, &references, r2.sequence.len()))
    }

    records
}

/// Given two lists of NAMs from R1 and R2, find the best location (preferably a proper pair).
/// This is used for mapping-only (PAF) mode and abundances output
fn get_best_paired_map_location(
    nams1: &[Nam],
    nams2: &[Nam],
    insert_size_distribution: &mut InsertSizeDistribution,
) -> (Option<Nam>, Option<Nam>) {
    let nam_pairs = get_best_scoring_nam_pairs(nams1, nams2, insert_size_distribution.mu, insert_size_distribution.sigma);
    if nam_pairs.is_empty() {
        return (None, None);
    }

    // Find first NAM pair that is a proper pair.
    // The first one is also the one with the highest score
    // since nam_pairs is sorted descending by score
    let best_joint_pair = nam_pairs
        .iter()
        .find(|&nam_pair| nam_pair.nam1.is_some() && nam_pair.nam2.is_some());

    let joint_score = if let Some(nam_pair) = best_joint_pair {
        nam_pair.nam1.as_ref().map_or(0, |nam| nam.score) +
        nam_pair.nam2.as_ref().map_or(0, |nam| nam.score)
    } else { 0 };

    // Get individual best scores.
    // nams1 and nams2 are also sorted descending by score.
    let best_individual_nam1 = nams1.first();
    let best_individual_nam2 = nams2.first();

    let individual_score =
        best_individual_nam1.map_or(0, |nam| nam.score) +
        best_individual_nam2.map_or(0, |nam| nam.score);

    // Divisor 2 is penalty for being mapped individually
    if joint_score > individual_score / 2 {
        let best_joint_pair= best_joint_pair.unwrap();
        let best = (best_joint_pair.nam1.clone(), best_joint_pair.nam2.clone());
        if insert_size_distribution.sample_size < 400 {
            insert_size_distribution.update(best.0.as_ref().unwrap().ref_start.abs_diff(best.1.as_ref().unwrap().ref_start));
        }

        best
    } else {
        (best_individual_nam1.cloned(), best_individual_nam2.cloned())
    }

    // TODO abundances
}
