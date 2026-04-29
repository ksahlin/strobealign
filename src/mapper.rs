use crate::aligner::Aligner;
use crate::aligner::{AlignmentInfo, hamming_align, hamming_distance};
use crate::chainer::Chainer;
use crate::cigar::{Cigar, CigarOperation};
use crate::details::Details;
use crate::index::StrobemerIndex;
use crate::insertsize::InsertSizeDistribution;
use crate::io::fasta::RefSequence;
use crate::io::record::SequenceRecord;
use crate::io::sam::{
    MREVERSE, MUNMAP, PAIRED, PROPER_PAIR, READ1, READ2, REVERSE, SECONDARY, SamRecord, UNMAP,
};
use crate::mcsstrategy::McsStrategy;
use crate::nam::{Nam, get_nams_by_chaining, sort_nams};
use crate::pairing::{PairedAlignments, get_paired_alignment, is_proper_pair};
use crate::piecewisealigner::remove_spurious_anchors;
use crate::read::Read;
use crate::revcomp::reverse_complement;
use crate::seeding::SeedingParameters;
use bumpalo::Bump;
use fastrand::Rng;
use std::cmp::{Reverse, min};
use std::time::Instant;

#[derive(Debug)]
pub struct MappingParameters {
    pub max_secondary: usize,
    pub dropoff_threshold: f32,
    pub rescue_distance: usize,
    pub max_tries: usize,
    pub mcs_strategy: McsStrategy,
    pub output_unmapped: bool,
    pub use_ssw: bool,
}

impl Default for MappingParameters {
    fn default() -> Self {
        MappingParameters {
            max_secondary: 0,
            dropoff_threshold: 0.5,
            rescue_distance: 100,
            max_tries: 20,
            mcs_strategy: McsStrategy::default(),
            output_unmapped: true,
            use_ssw: false,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Alignment {
    pub reference_id: usize,
    pub ref_start: usize,
    pub cigar: Cigar,
    pub edit_distance: usize,
    pub soft_clip_left: usize,
    pub soft_clip_right: usize,
    pub score: u32,
    pub length: usize,
    pub is_revcomp: bool,
    /// Whether a gapped alignment function was used to obtain this alignment
    /// (even if true, the alignment can still be without gaps)
    pub gapped: bool,
    pub rescued: bool,
}

impl Alignment {
    fn global_edit_distance(&self) -> usize {
        self.edit_distance + self.soft_clip_left + self.soft_clip_right
    }
}

/// Conversion of an Alignment into a SamRecord
#[derive(Default)]
pub struct SamOutput {
    cigar_eqx: bool,
    details: bool,
    rg_id: Option<String>,
    fastq_comments: bool,
}

impl SamOutput {
    pub fn new(
        details: bool,
        cigar_eqx: bool,
        rg_id: Option<String>,
        fastq_comments: bool,
    ) -> Self {
        SamOutput {
            cigar_eqx,
            details,
            rg_id,
            fastq_comments,
        }
    }

    fn make_record(
        &self,
        alignment: Option<&Alignment>,
        references: &[RefSequence],
        record: &SequenceRecord,
        mapq: u8,
        is_primary: bool,
        details: Details,
    ) -> SamRecord {
        match alignment {
            Some(alignment) => self.make_mapped_record(
                alignment,
                references,
                record,
                mapq,
                is_primary,
                details.clone(),
            ),
            None => self.make_unmapped_record(record, details.clone()),
        }
    }

    /// Convert the alignment into a SamRecord
    fn make_mapped_record(
        &self,
        alignment: &Alignment,
        references: &[RefSequence],
        record: &SequenceRecord,
        mut mapq: u8,
        is_primary: bool,
        details: Details,
    ) -> SamRecord {
        let mut flags = 0;

        if alignment.is_revcomp {
            flags |= REVERSE;
        }
        if !is_primary {
            mapq = 0;
            flags |= SECONDARY;
        }

        let query_sequence = if alignment.is_revcomp {
            reverse_complement(&record.sequence)
        } else {
            record.sequence.clone()
        };
        let query_qualities = if alignment.is_revcomp {
            if let Some(mut qualities) = record.qualities.clone() {
                qualities.reverse();
                Some(qualities)
            } else {
                None
            }
        } else {
            record.qualities.clone()
        };
        let mut cigar = Cigar::new();
        cigar.push(CigarOperation::Softclip, alignment.soft_clip_left);
        cigar.extend(&alignment.cigar);
        cigar.push(CigarOperation::Softclip, alignment.soft_clip_right);
        let reference_name = Some(references[alignment.reference_id].name.clone());
        let details = if self.details { Some(details) } else { None };
        let cigar = if self.cigar_eqx {
            Some(cigar)
        } else {
            Some(cigar.with_m())
        };
        let extra = if self.fastq_comments {
            record.comment.clone()
        } else {
            None
        };
        SamRecord {
            query_name: record.name.clone(),
            flags,
            reference_name,
            pos: Some(alignment.ref_start as u32),
            mapq,
            cigar,
            query_sequence: Some(query_sequence),
            query_qualities,
            edit_distance: Some(alignment.edit_distance as u32),
            alignment_score: Some(alignment.score),
            details,
            rg_id: self.rg_id.clone(),
            extra,
            ..SamRecord::default()
        }
    }

    fn make_unmapped_record(&self, record: &SequenceRecord, details: Details) -> SamRecord {
        let details = if self.details { Some(details) } else { None };
        let extra = if self.fastq_comments {
            record.comment.clone()
        } else {
            None
        };
        SamRecord {
            query_name: record.name.clone(),
            flags: UNMAP,
            query_sequence: Some(record.sequence.clone()),
            query_qualities: record.qualities.clone(),
            details,
            rg_id: self.rg_id.clone(),
            extra,
            ..SamRecord::default()
        }
    }

    pub fn make_unmapped_pair(
        &self,
        record1: &SequenceRecord,
        record2: &SequenceRecord,
        details1: &Details,
        details2: &Details,
    ) -> [SamRecord; 2] {
        let mut sam_records = [
            self.make_unmapped_record(record1, details1.clone()),
            self.make_unmapped_record(record2, details2.clone()),
        ];
        sam_records[0].flags |= PAIRED | READ1 | MUNMAP;
        sam_records[1].flags |= PAIRED | READ2 | MUNMAP;

        sam_records
    }

    // alignment: &Alignment, references: &[RefSequence], record: &SequenceRecord, mut mapq: u8, is_primary: bool, details: Details) -> SamRecord {
    fn make_paired_records(
        &self,
        alignments: [Option<&Alignment>; 2],
        references: &[RefSequence],
        records: [&SequenceRecord; 2],
        mapq: [u8; 2],
        details1: &Details,
        details2: &Details,
        is_primary: bool,
        is_proper: bool,
    ) -> [SamRecord; 2] {
        // Create single-end records
        let mut sam_records = [
            self.make_record(
                alignments[0],
                references,
                records[0],
                mapq[0],
                is_primary,
                details1.clone(),
            ),
            self.make_record(
                alignments[1],
                references,
                records[1],
                mapq[1],
                is_primary,
                details2.clone(),
            ),
        ];
        // Then make them paired
        sam_records[0].flags |= READ1;
        sam_records[1].flags |= READ2;

        for i in 0..2 {
            // Flags
            sam_records[i].flags |= PAIRED;
            if is_proper {
                sam_records[i].flags |= PROPER_PAIR;
            }
            if let Some(mate) = alignments[1 - i] {
                if mate.is_revcomp {
                    sam_records[i].flags |= MREVERSE;
                }
                // PNEXT (position of mate)
                sam_records[i].mate_pos = Some(mate.ref_start as u32);

                // RNEXT (reference name of mate)
                sam_records[i].mate_reference_name =
                    Some(references[mate.reference_id].name.clone());
                if let Some(this) = alignments[i] {
                    // both aligned

                    if this.reference_id == mate.reference_id {
                        // aligned to same reference
                        sam_records[i].mate_reference_name = Some("=".to_string());

                        // TLEN
                        let template_length = if mate.ref_start > this.ref_start {
                            (mate.ref_start - this.ref_start + mate.length) as isize
                        } else if this.ref_start > mate.ref_start {
                            -((this.ref_start - mate.ref_start + this.length) as isize)
                        } else if i == 0 {
                            -(mate.length as isize)
                        } else {
                            this.length as isize
                        };
                        sam_records[i].template_len = Some(template_length as i32);
                    }
                } else {
                    // The SAM specification recommends: "For a[n] unmapped paired-end or
                    // mate-pair read whose mate is mapped, the unmapped read should have
                    // RNAME and POS identical to its mate."
                    sam_records[i].mate_reference_name = Some("=".to_string());
                    sam_records[i].reference_name =
                        Some(references[mate.reference_id].name.clone());
                    sam_records[i].pos = Some(mate.ref_start as u32);
                }
            } else {
                // mate unmapped

                sam_records[i].flags |= MUNMAP;
                // Set RNAME and POS identical to mate
                if let Some(this) = alignments[i] {
                    sam_records[i].mate_reference_name = Some("=".to_string());
                    sam_records[i].mate_pos = Some(this.ref_start as u32);
                }
            }
        }

        sam_records
    }
}

/// Align a single-end read to the reference and return SAM records
pub fn align_single_end_read(
    record: &SequenceRecord,
    index: &StrobemerIndex,
    references: &[RefSequence],
    mapping_parameters: &MappingParameters,
    sam_output: &SamOutput,
    chainer: &Chainer,
    aligner: &Aligner,
    rng: &mut Rng,
    arena: &Bump,
) -> (Vec<SamRecord>, Details) {
    let (mut nam_details, mut nams) = get_nams_by_chaining(
        &record.sequence,
        index,
        chainer,
        mapping_parameters.rescue_distance,
        mapping_parameters.mcs_strategy,
        arena,
    );
    nam_details.time_sort_nams = sort_nams(&mut nams, rng);
    let mut details: Details = nam_details.into();

    let timer = Instant::now();
    if nams.is_empty() {
        return (
            vec![sam_output.make_unmapped_record(record, details.clone())],
            details,
        );
    }
    let mut sam_records = Vec::new();
    let mut alignments = Vec::new();
    let nam_max = nams[0].clone();
    let mut best_edit_distance = usize::MAX;
    let mut best_score = 0;
    let mut best_index = 0;
    let mut second_best_score = 0;
    let mut alignments_with_best_score = 0;
    let mut best_alignment = None;

    let k = index.k();
    let read = Read::new(&record.sequence);
    for (tries, nam) in nams
        .iter_mut()
        .take(mapping_parameters.max_tries)
        .enumerate()
    {
        let score_dropoff = nam.anchors.len() as f32 / nam_max.anchors.len() as f32;

        if (tries > 1 && best_edit_distance == 0)
            || score_dropoff < mapping_parameters.dropoff_threshold
        {
            break;
        }
        let consistent_nam = nam.is_consistent(&read, references, k);
        if !consistent_nam {
            details.inconsistent_nams += 1;
            continue;
        }
        let Some(alignment) = extend_seed(
            aligner,
            nam,
            references,
            &read,
            consistent_nam,
            mapping_parameters.use_ssw,
        ) else {
            continue;
        };

        // outputting Piecewise vs SSW alignments for debugging
        // if log::log_enabled!(log::Level::Trace) {
        //     let (mut ssw, mut pw) = if !mapping_parameters.use_ssw {
        //         (
        //             extend_seed(aligner, nam, references, &read, consistent_nam, true).unwrap(),
        //             alignment.clone(),
        //         )
        //     } else {
        //         (
        //             alignment.clone(),
        //             extend_seed(aligner, nam, references, &read, consistent_nam, false).unwrap(),
        //         )
        //     };
        //     // manually adds the soft clips
        //     let mut cigar = Cigar::new();
        //     cigar.push(CigarOperation::Softclip, pw.soft_clip_left);
        //     cigar.extend(&pw.cigar);
        //     cigar.push(CigarOperation::Softclip, pw.soft_clip_right);
        //     pw.cigar = cigar;
        //
        //     let mut cigar = Cigar::new();
        //     cigar.push(CigarOperation::Softclip, ssw.soft_clip_left);
        //     cigar.extend(&ssw.cigar);
        //     cigar.push(CigarOperation::Softclip, ssw.soft_clip_right);
        //     ssw.cigar = cigar;
        //
        //     trace!("Alignment:[{:?},SSW:{:?},PW:{:?}]", nam.clone(), ssw, pw);
        // }

        details.tried_alignment += 1;
        details.gapped += alignment.gapped as usize;

        if alignment.score >= best_score {
            second_best_score = best_score;
            let mut update_best = false;
            if alignment.score > best_score {
                alignments_with_best_score = 1;
                update_best = true;
            } else {
                assert_eq!(alignment.score, best_score);
                // Two or more alignments have the same best score - count them
                alignments_with_best_score += 1;

                // Pick one randomly using reservoir sampling
                if rng.u32(..alignments_with_best_score) == 0 {
                    update_best = true;
                }
            }
            if update_best {
                best_score = alignment.score;
                best_alignment = Some(alignment.clone());
                best_index = alignments.len();
                if mapping_parameters.max_secondary == 0 {
                    best_edit_distance = alignment.global_edit_distance();
                }
            }
        } else if alignment.score > second_best_score {
            second_best_score = alignment.score;
        }
        if mapping_parameters.max_secondary > 0 {
            alignments.push(alignment);
        }
    }
    if best_alignment.is_none() {
        return (
            vec![sam_output.make_unmapped_record(record, details.clone())],
            details,
        );
    }
    let mapq = (60 * (best_score - second_best_score)).div_ceil(best_score) as u8;

    let best_alignment = best_alignment.unwrap();
    let mut is_primary = true;
    sam_records.push(sam_output.make_mapped_record(
        &best_alignment,
        references,
        record,
        mapq,
        is_primary,
        details.clone(),
    ));

    // Secondary alignments
    if mapping_parameters.max_secondary > 0 {
        // Remove the primary alignment
        alignments.swap_remove(best_index);

        // Highest score first
        alignments.sort_by_key(|k| Reverse(k.score));

        // Output secondary alignments
        //let max_out = min(alignments.len(), mapping_parameters.max_secondary + 1);
        for alignment in alignments.iter().take(mapping_parameters.max_secondary) {
            if alignment.score.saturating_sub(best_score)
                > 2 * aligner.scores.mismatch as u32 + aligner.scores.gap_open as u32
            {
                break;
            }
            is_primary = false;
            // TODO .clone()
            sam_records.push(sam_output.make_mapped_record(
                alignment,
                references,
                record,
                mapq,
                is_primary,
                details.clone(),
            ));
        }
    }
    details.time_extend = timer.elapsed().as_secs_f64();

    (sam_records, details)
}

/// Extend a NAM so that it covers the entire read and return the resulting
/// alignment.
pub fn extend_seed(
    aligner: &Aligner,
    nam: &mut Nam,
    references: &[RefSequence],
    read: &Read,
    consistent_nam: bool,
    use_ssw: bool,
) -> Option<Alignment> {
    let query = if nam.is_revcomp {
        read.rc()
    } else {
        read.seq()
    };
    let refseq = &references[nam.ref_id].sequence;

    let projected_ref_start = nam.ref_start.saturating_sub(nam.query_start);
    let projected_ref_end = min(nam.ref_end + query.len() - nam.query_end, refseq.len());

    // TODO ugly
    let mut info = AlignmentInfo {
        cigar: Default::default(),
        edit_distance: 0,
        ref_start: 0,
        ref_end: 0,
        query_start: 0,
        query_end: 0,
        score: 0,
    };
    let mut result_ref_start = 0;
    let mut gapped = true;
    if projected_ref_start + query.len() == projected_ref_end && consistent_nam {
        let ref_segm_ham = &refseq[projected_ref_start..projected_ref_end];
        if let Some(hamming_dist) = hamming_distance(query, ref_segm_ham) {
            if (hamming_dist as f32 / query.len() as f32) < 0.05 {
                // ungapped worked fine, no need to do gapped alignment
                info = hamming_align(
                    query,
                    ref_segm_ham,
                    aligner.scores.match_,
                    aligner.scores.mismatch,
                    aligner.scores.end_bonus,
                )
                .expect("hamming_dist was successful, this should be as well");
                result_ref_start = projected_ref_start + info.ref_start;
                gapped = false;
            }
        }
    }
    if gapped {
        let padding = read.len() / 10;
        if use_ssw {
            let ref_start = projected_ref_start.saturating_sub(padding);
            let ref_end = min(projected_ref_end + padding, refseq.len());
            let segment = &refseq[ref_start..ref_end];
            info = aligner.align(query, segment)?;
            result_ref_start = ref_start + info.ref_start;
        } else {
            remove_spurious_anchors(&mut nam.anchors);
            info = aligner.align_piecewise(query, refseq, &nam.anchors, padding)?;
            result_ref_start = info.ref_start;
        }
    }
    Some(Alignment {
        cigar: info.cigar.clone(),
        edit_distance: info.edit_distance,
        soft_clip_left: info.query_start,
        soft_clip_right: query.len() - info.query_end,
        score: info.score,
        ref_start: result_ref_start,
        length: info.ref_span(),
        is_revcomp: nam.is_revcomp,
        reference_id: nam.ref_id,
        gapped,
        rescued: false,
    })
}

/// Align a paired-end read pair to the reference and return SAM records
pub fn align_paired_end_read(
    r1: &SequenceRecord,
    r2: &SequenceRecord,
    index: &StrobemerIndex,
    references: &[RefSequence],
    mapping_parameters: &MappingParameters,
    sam_output: &SamOutput,
    seeding_parameters: &SeedingParameters,
    insert_size_distribution: &mut InsertSizeDistribution,
    chainer: &Chainer,
    aligner: &Aligner,
    rng: &mut Rng,
    arena: &Bump,
) -> (Vec<SamRecord>, Details) {
    let (nam_details1, nams1) = get_nams_by_chaining(
        &r1.sequence,
        index,
        chainer,
        mapping_parameters.rescue_distance,
        mapping_parameters.mcs_strategy,
        arena,
    );
    let (nam_details2, nams2) = get_nams_by_chaining(
        &r2.sequence,
        index,
        chainer,
        mapping_parameters.rescue_distance,
        mapping_parameters.mcs_strategy,
        arena,
    );
    let (mut details1, mut details2): (Details, Details) =
        (nam_details1.into(), nam_details2.into());

    if nams1.is_empty() && nams2.is_empty() {
        let mut details_both = details1.clone();
        details_both += details2.clone();
        return (
            sam_output
                .make_unmapped_pair(r1, r2, &details1, &details2)
                .into(),
            details_both,
        );
    }

    let timer = Instant::now();
    let read1 = Read::new(&r1.sequence);
    let read2 = Read::new(&r2.sequence);

    let (paired_alignments, alignments1, alignments2) = get_paired_alignment(
        aligner,
        nams1,
        nams2,
        mapping_parameters.max_tries,
        mapping_parameters.dropoff_threshold,
        references,
        &read1,
        &read2,
        &mut details1,
        &mut details2,
        insert_size_distribution.mu,
        insert_size_distribution.sigma,
        seeding_parameters.syncmer.k,
        rng,
    );

    let mut sam_records = Vec::new();
    match get_best_paired_alignment(
        &paired_alignments,
        &alignments1,
        &alignments2,
        insert_size_distribution,
    ) {
        BestPairedAlignment::Pair(a1, a2, best_score) => {
            let mapq = alignment_quality(&paired_alignments, |p| p.score as f32);
            let proper = is_proper_pair(
                a1,
                a2,
                insert_size_distribution.mu,
                insert_size_distribution.sigma,
            );

            // Primary
            sam_records.extend(sam_output.make_paired_records(
                [Some(a1), Some(a2)],
                references,
                [r1, r2],
                [mapq, mapq],
                &details1,
                &details2,
                true,
                proper,
            ));

            // Secondaries
            let secondary_dropoff = (2 * aligner.scores.mismatch + aligner.scores.gap_open) as f64;
            for pair in paired_alignments
                .iter()
                .skip(1)
                .take(mapping_parameters.max_secondary)
            {
                if (best_score - pair.score) >= secondary_dropoff {
                    break;
                }
                let proper = is_proper_pair(
                    &pair.alignment1,
                    &pair.alignment2,
                    insert_size_distribution.mu,
                    insert_size_distribution.sigma,
                );
                sam_records.extend(sam_output.make_paired_records(
                    [Some(&pair.alignment1), Some(&pair.alignment2)],
                    references,
                    [r1, r2],
                    [0, 0],
                    &details1,
                    &details2,
                    false,
                    proper,
                ));
            }
        }
        BestPairedAlignment::Individual(a1, a2) => {
            let mapq1 = alignment_quality(&alignments1, |a| a.score as f32);
            let mapq2 = alignment_quality(&alignments2, |a| a.score as f32);
            let proper = if let (Some(a1), Some(a2)) = (a1, a2) {
                is_proper_pair(
                    a1,
                    a2,
                    insert_size_distribution.mu,
                    insert_size_distribution.sigma,
                )
            } else {
                false
            };
            sam_records.extend(sam_output.make_paired_records(
                [a1, a2],
                references,
                [r1, r2],
                [mapq1, mapq2],
                &details1,
                &details2,
                true,
                proper,
            ));
        }
    }
    let mut details_both = details1.clone();
    details_both += details2.clone();
    details_both.time_extend = timer.elapsed().as_secs_f64();
    (sam_records, details_both)
}

enum BestPairedAlignment<'a> {
    /// Best proper paired alignment
    Pair(&'a Alignment, &'a Alignment, f64),
    /// Independent best alignments
    Individual(Option<&'a Alignment>, Option<&'a Alignment>),
}

/// Choose between:
/// - the best individual Alignments
/// - the best proper pair of Alignments
///
/// Also updates the insert size distribution using confident pairs.
///
/// For paired-end alignment extension only
fn get_best_paired_alignment<'a>(
    paired_alignments: &'a [PairedAlignments],
    alignments1: &'a [Alignment],
    alignments2: &'a [Alignment],
    insert_size_distribution: &mut InsertSizeDistribution,
) -> BestPairedAlignment<'a> {
    if let Some(PairedAlignments {
        score,
        alignment1,
        alignment2,
    }) = paired_alignments.first()
    {
        // Update insert size
        if insert_size_distribution.sample_size < 400 {
            insert_size_distribution.update(alignment1.ref_start.abs_diff(alignment2.ref_start));
        }

        BestPairedAlignment::Pair(alignment1, alignment2, *score)
    } else {
        BestPairedAlignment::Individual(alignments1.first(), alignments2.first())
    }
}

/// Return mapping quality for the top NAM
pub fn mapping_quality(nams: &[Nam]) -> u8 {
    if nams.len() <= 1 {
        return 60;
    }
    let s1 = nams[0].score;
    let s2 = nams[1].score;
    // from minimap2: MAPQ = 40(1−s2/s1) ·min{1,|M|/10} · log s1
    let min_matches = min(nams[0].anchors.len(), 10) as f32 / 10.0;
    let uncapped_mapq = 40.0 * (1.0 - s2 / s1) * min_matches * s1.ln();

    uncapped_mapq.min(60.0) as u8
}

/// Return mapping quality for the top Alignment or Paired Alignment
fn alignment_quality<T, F>(items: &[T], score: F) -> u8
where
    F: Fn(&T) -> f32,
{
    if items.len() <= 1 {
        return 60;
    }

    let score1 = score(&items[0]);
    let score2 = score(&items[1]);

    if score1 > 0.0 && score2 > 0.0 {
        (score1 - score2).min(60.0) as u8
    } else if score1 > 0.0 && score2 <= 0.0 {
        60
    } else {
        1
    }
}
