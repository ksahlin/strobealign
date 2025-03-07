use std::cmp::{min, Reverse};
use std::collections::{HashMap, HashSet};
use std::mem;
use fastrand::Rng;
use memchr::memmem;
use crate::aligner::{hamming_align, hamming_distance, AlignmentInfo};
use crate::cigar::{Cigar, CigarOperation};
use crate::fasta::RefSequence;
use crate::index::{IndexParameters, StrobemerIndex};
use crate::nam::{reverse_nam_if_needed, Nam};
use crate::revcomp::reverse_complement;
use crate::strobes::RandstrobeIterator;
use crate::syncmers::SyncmerIterator;
use crate::sam::{strip_suffix, SamRecord, MREVERSE, MUNMAP, PAIRED, PROPER_PAIR, READ1, READ2, REVERSE, SECONDARY, UNMAP};
use crate::read::Read;
use crate::aligner::Aligner;
use crate::details::Details;
use crate::fastq::SequenceRecord;
use crate::insertsize::InsertSizeDistribution;
use crate::math::normal_pdf;
use crate::nam;

pub struct MappingParameters {
    pub r: usize,
    pub max_secondary: usize,
    pub dropoff_threshold: f32,
    pub rescue_level: usize,
    pub max_tries: usize,
    pub output_unmapped: bool,
}

impl Default for MappingParameters {
    fn default() -> Self {
        MappingParameters {
            r: 150,
            max_secondary: 0,
            dropoff_threshold: 0.5,
            rescue_level: 2,
            max_tries: 20,
            output_unmapped: true,
        }
    }
}

#[derive(Debug, Clone)]
struct Alignment {
    reference_id: usize,
    ref_start: usize,
    cigar: Cigar,
    edit_distance: usize,
    soft_clip_left: usize,
    soft_clip_right: usize,
    score: u32,
    length: usize,
    is_revcomp: bool,
    /// Whether a gapped alignment function was used to obtain this alignment
    /// (even if true, the alignment can still be without gaps)
    gapped: bool,
}

impl Alignment {
    fn global_edit_distance(&self) -> usize {
        self.edit_distance + self.soft_clip_left + self.soft_clip_right
    }
}

#[derive(Debug)]
pub struct QueryRandstrobe {
    pub hash: u64,
    pub start: usize,
    pub end: usize,
    pub is_revcomp: bool,
}

/// Generate randstrobes for a query sequence and its reverse complement.
/// TODO move to strobes.rs?
pub fn randstrobes_query(seq: &[u8], parameters: &IndexParameters) -> Vec<QueryRandstrobe> {
    let mut randstrobes= Vec::<QueryRandstrobe>::new();
    if seq.len() < parameters.randstrobe.w_max {
        return randstrobes;
    }

    // TODO
    // For the reverse complement, we could re-use the syncmers of the forward
    // sequence because canonical syncmers are invariant under reverse
    // complementing. Only the coordinates need to be adjusted.

    let seq_rc = reverse_complement(seq);
    for (s, is_reverse) in [(seq, false), (&seq_rc, true)] {
        // Generate randstrobes for the forward sequence
        let mut syncmer_iter = SyncmerIterator::new(s, parameters.syncmer.k, parameters.syncmer.s, parameters.syncmer.t);
        let randstrobe_iter = RandstrobeIterator::new(&mut syncmer_iter, &parameters.randstrobe);

        for randstrobe in randstrobe_iter {
            randstrobes.push(
                QueryRandstrobe {
                    hash: randstrobe.hash,
                    start: randstrobe.strobe1_pos,
                    end: randstrobe.strobe2_pos + parameters.syncmer.k,
                    is_revcomp: is_reverse
                }
            );
        }
    }

    randstrobes
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
    pub fn new(details: bool, cigar_eqx: bool, rg_id: Option<String>, fastq_comments: bool) -> Self {
        SamOutput {
            cigar_eqx,
            details,
            rg_id,
            fastq_comments,
        }
    }

    fn make_record(
        &self, alignment: Option<&Alignment>, references: &[RefSequence], record: &SequenceRecord, mapq: u8, is_primary: bool, details: Details
    ) -> SamRecord {
        match alignment {
            Some(alignment) => self.make_mapped_record(alignment, references, record, mapq, is_primary, details.clone()),
            None => self.make_unmapped_record(record, details.clone())
        }
    }

    /// Convert the alignment into a SamRecord
    fn make_mapped_record(
        &self, alignment: &Alignment, references: &[RefSequence], record: &SequenceRecord, mut mapq: u8, is_primary: bool, details: Details
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
            let mut rev = record.qualities.clone();
            rev.reverse();
            rev
        } else {
            record.qualities.clone()
        };
        let mut cigar = Cigar::new();
        cigar.push(CigarOperation::Softclip, alignment.soft_clip_left);
        cigar.extend(&alignment.cigar);
        cigar.push(CigarOperation::Softclip, alignment.soft_clip_right);
        let reference_name = Some(references[alignment.reference_id].name.clone());
        let details = if self.details { Some(details) } else { None };
        let cigar = if self.cigar_eqx { Some(cigar) } else { Some(cigar.with_m()) };
        let extra = if self.fastq_comments { record.comment.clone() } else { None };
        SamRecord {
            query_name: strip_suffix(&record.name).into(),
            flags,
            reference_name,
            pos: Some(alignment.ref_start as u32),
            mapq,
            cigar,
            query_sequence: Some(query_sequence),
            query_qualities: Some(query_qualities),
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
        let extra = if self.fastq_comments { record.comment.clone() } else { None };
        SamRecord {
            query_name: strip_suffix(&record.name).into(),
            flags: UNMAP,
            query_sequence: Some(record.sequence.clone()),
            query_qualities: Some(record.qualities.clone()),
            details,
            rg_id: self.rg_id.clone(),
            extra,
            ..SamRecord::default()
        }
    }

    fn make_unmapped_pair(&self, records: [&SequenceRecord; 2], details: &[Details; 2]) -> [SamRecord; 2] {
        let mut sam_records = [
            self.make_unmapped_record(records[0], details[0].clone()),
            self.make_unmapped_record(records[1], details[1].clone()),
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
        details: &[Details; 2],
        is_primary: bool,
        is_proper: bool,
    ) -> [SamRecord; 2] {
        // Create single-end records
        let mut sam_records = [
            self.make_record(alignments[0], references, records[0], mapq[0], is_primary, details[0].clone()),
            self.make_record(alignments[1], references, records[1], mapq[1], is_primary, details[1].clone()),
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
            if let Some(mate) = alignments[1-i] {
                if mate.is_revcomp {
                    sam_records[i].flags |= MREVERSE;
                }
                // RNEXT (reference name of mate)
                if let Some(this) = alignments[i] {
                    // both aligned
                    if this.reference_id == mate.reference_id {
                        sam_records[i].mate_reference_name = Some("=".to_string());

                        // TLEN
                        let template_length =
                            if mate.ref_start > this.ref_start {
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
                    sam_records[i].reference_name = Some(references[mate.reference_id].name.clone());
                    sam_records[i].pos = Some(mate.ref_start as u32);
                }
                // PNEXT (position of mate)
                sam_records[i].mate_pos = Some(mate.ref_start as u32);
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
    aligner: &Aligner,
    rng: &mut Rng,
) -> Vec<SamRecord> {
    let mut details = Details::default();

    let (nam_details, mut nams) = nam::get_nams(&record.sequence, index, mapping_parameters.rescue_level, rng);
    details.nam_rescue = nam_details.nam_rescue;
    details.nams = nams.len();
    details.n_hits = nam_details.n_hits;
    details.n_rescue_hits = nam_details.n_rescue_hits;

    // Timer extend_timer;
    if nams.is_empty() {
        return vec![sam_output.make_unmapped_record(record, details)];
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

    let k = index.parameters.syncmer.k;
    let read = Read::new(&record.sequence);
    for (tries, nam) in nams.iter_mut().enumerate() {
        let score_dropoff = nam.n_hits as f32 / nam_max.n_hits as f32;

        // TODO iterate over slice of nams instead of tracking tries
        if tries >= mapping_parameters.max_tries || (tries > 1 && best_edit_distance == 0) || score_dropoff < mapping_parameters.dropoff_threshold {
            break;
        }
        let consistent_nam = reverse_nam_if_needed(nam, &read, references, k);
        details.inconsistent_nams += (!consistent_nam) as usize;
        let alignment = extend_seed(aligner, nam, references, &read, consistent_nam);
        if alignment.is_none() {
            continue;
        }
        let alignment = alignment.unwrap();
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
                best_index = tries;
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
        return vec![sam_output.make_unmapped_record(record, details)];
    }
    let mapq = ((60 * (best_score - second_best_score) + best_score - 1) / best_score) as u8;

    let best_alignment = best_alignment.unwrap();
    let mut is_primary = true;
    sam_records.push(
        sam_output.make_mapped_record(&best_alignment, references, record, mapq, is_primary, details.clone())
    );

    // Secondary alignments
    if mapping_parameters.max_secondary > 0 {
        // Remove the primary alignment
        alignments.swap_remove(best_index);

        // Highest score first
        alignments.sort_by_key(|k| Reverse(k.score));

        // Output secondary alignments
        //let max_out = min(alignments.len(), mapping_parameters.max_secondary + 1);
        for (i, alignment) in alignments.iter().enumerate() {

            if i >= mapping_parameters.max_secondary || alignment.score - best_score > 2 * aligner.scores.mismatch as u32 + aligner.scores.gap_open as u32 {
                break;
            }
            is_primary = false;
            // TODO .clone()
            sam_records.push(sam_output.make_mapped_record(alignment, references, record, mapq, is_primary, details.clone()));
        }
    }
    // statistics.tot_extend += extend_timer.duration();
    // statistics += details;
    sam_records
}

/// Extend a NAM so that it covers the entire read and return the resulting
/// alignment.
fn extend_seed(
    aligner: &Aligner,
    nam: &Nam,
    references: &[RefSequence],
    read: &Read,
    consistent_nam: bool,
) -> Option<Alignment> {
    let query = if nam.is_revcomp { read.rc() } else { read.seq() };
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
                info = hamming_align(query, ref_segm_ham, aligner.scores.match_, aligner.scores.mismatch, aligner.scores.end_bonus).expect(
                    "hamming_dist was successful, this should be as well"
                );
                result_ref_start = projected_ref_start + info.ref_start;
                gapped = false;
            }
        }
    }
    if gapped {
        let ref_start = projected_ref_start.saturating_sub(50);
        let ref_end = min(projected_ref_end + 50, refseq.len());
        let segment = &refseq[ref_start..ref_end];
        info = aligner.align(query, segment)?;
        result_ref_start = ref_start + info.ref_start;
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
    })
}

// TODO alignment statistics
/// Align a paired-end read pair to the reference and return SAM records
pub fn align_paired_end_read(
    r1: &SequenceRecord,
    r2: &SequenceRecord,
    index: &StrobemerIndex,
    references: &[RefSequence],
    mapping_parameters: &MappingParameters,
    sam_output: &SamOutput,
    index_parameters: &IndexParameters,
    insert_size_distribution: &mut InsertSizeDistribution,
    aligner: &Aligner,
    rng: &mut Rng,
) -> Vec<SamRecord> {
    let mut details = [Details::default(), Details::default()];
    let mut nams_pair = [vec![], vec![]];

    for is_revcomp in [0, 1] {
        let record = if is_revcomp == 0 { r1 } else { r2 };
        let (nam_details, nams) = nam::get_nams(&record.sequence, index, mapping_parameters.rescue_level, rng);
        details[is_revcomp].nam_rescue = nam_details.nam_rescue;
        details[is_revcomp].nams = nams.len();
        details[is_revcomp].n_hits = nam_details.n_hits;
        details[is_revcomp].n_rescue_hits = nam_details.n_rescue_hits;
        nams_pair[is_revcomp] = nams;
    }

    // Timer extend_timer;
    let read1 = Read::new(&r1.sequence); // TODO pass r1, r2 to extend_paired_seeds instead
    let read2 = Read::new(&r2.sequence);
    let alignment_pairs = extend_paired_seeds(
        aligner, &mut nams_pair, &read1, &read2,
        index_parameters.syncmer.k, references, &mut details,
        mapping_parameters.dropoff_threshold, insert_size_distribution,
        mapping_parameters.max_tries
    );

    let mut sam_records = Vec::new();

    match alignment_pairs {
        // Typical case: both reads map uniquely and form a proper pair.
        // Then the mapping quality is computed based on the NAMs.
        AlignedPairs::Proper((alignment1, alignment2)) => {
            let is_proper = is_proper_pair(Some(&alignment1), Some(&alignment2), insert_size_distribution.mu, insert_size_distribution.sigma);
            if is_proper
                && insert_size_distribution.sample_size < 400
                && alignment1.edit_distance + alignment2.edit_distance < 3
            {
                insert_size_distribution.update(alignment1.ref_start.abs_diff(alignment2.ref_start));
            }

            let mapq1 = proper_pair_mapq(&nams_pair[0]);
            let mapq2 = proper_pair_mapq(&nams_pair[1]);

            details[0].best_alignments = 1;
            details[1].best_alignments = 1;
            let is_primary = true;

            sam_records.extend(
                sam_output.make_paired_records([Some(&alignment1), Some(&alignment2)], references, [r1, r2], [mapq1, mapq2], &details, is_primary, is_proper)
            );
        },
        AlignedPairs::WithScores(mut alignment_pairs) => {
            alignment_pairs.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());
            deduplicate_scored_pairs(&mut alignment_pairs);

            // If there are multiple top-scoring alignments (all with the same score),
            // pick one randomly and move it to the front.
            let i = count_best_alignment_pairs(&alignment_pairs);
            details[0].best_alignments = i;
            details[1].best_alignments = i;
            if i > 1 {
                let random_index = rng.usize(..i);
                alignment_pairs.swap(0, random_index);
            }

            let secondary_dropoff = 2 * aligner.scores.mismatch + aligner.scores.gap_open;
            sam_records.extend(aligned_pairs_to_sam(
                &alignment_pairs,
                sam_output,
                references,
                mapping_parameters.max_secondary,
                secondary_dropoff as f64,
                r1,
                r2,
                insert_size_distribution.mu,
                insert_size_distribution.sigma,
                &details
            ));
        }
    }
    // TODO
    // statistics.tot_extend += extend_timer.duration();
    // statistics += details[0];
    // statistics += details[1];

    sam_records
}

#[derive(Debug)]
enum AlignedPairs {
    Proper((Alignment, Alignment)),
    WithScores(Vec<ScoredAlignmentPair>),
}

/// Given two lists of NAMs for the two reads in a pair, pair them up and
/// compute base-level alignments
fn extend_paired_seeds(
    aligner: &Aligner,
    nams: &mut [Vec<Nam>; 2],
    read1: &Read,
    read2: &Read,
    k: usize,
    references: &[RefSequence],
    details: &mut [Details; 2],
    dropoff: f32,
    insert_size_distribution: &InsertSizeDistribution,
    max_tries: usize,
) -> AlignedPairs {
    let mu = insert_size_distribution.mu;
    let sigma = insert_size_distribution.sigma;

    if nams[0].is_empty() && nams[1].is_empty() {
         // None of the reads have any NAMs
        return AlignedPairs::WithScores(vec![]);
    }

    if !nams[0].is_empty() && nams[1].is_empty() {
        // Only read 1 has NAMS: attempt to rescue read 2
        return AlignedPairs::WithScores(rescue_read(
            read2,
            read1,
            aligner,
            references,
            &mut nams[0],
            max_tries,
            dropoff,
            details,
            k,
            mu,
            sigma
        ));
    }

    if nams[0].is_empty() && !nams[1].is_empty() {
        // Only read 2 has NAMS: attempt to rescue read 1
        let mut swapped_details = [details[0].clone(), details[1].clone()];
        let mut pairs = rescue_read(
            read1,
            read2,
            aligner,
            references,
            &mut nams[1],
            max_tries,
            dropoff,
            &mut swapped_details,
            k,
            mu,
            sigma
        );
        details[0] += swapped_details[1].clone();
        details[1] += swapped_details[0].clone();
        for pair in &mut pairs {
            mem::swap(&mut pair.alignment1, &mut pair.alignment2);
        }

        return AlignedPairs::WithScores(pairs);
    }

    // Both reads have NAMs
    assert!(!nams[0].is_empty() && !nams[1].is_empty());

    // Deal with the typical case that both reads map uniquely and form a proper pair
    if top_dropoff(&nams[0]) < dropoff && top_dropoff(&nams[1]) < dropoff && is_proper_nam_pair(&nams[0][0], &nams[1][0], mu, sigma) {
        let mut n_max1 = nams[0][0].clone();
        let mut n_max2 = nams[1][0].clone();

        let consistent_nam1 = reverse_nam_if_needed(&mut n_max1, read1, references, k);
        details[0].inconsistent_nams += !consistent_nam1 as usize;
        let consistent_nam2 = reverse_nam_if_needed(&mut n_max2, read2, references, k);
        details[1].inconsistent_nams += !consistent_nam2 as usize;

        let alignment1 = extend_seed(aligner, &n_max1, references, read1, consistent_nam1);
        let alignment2 = extend_seed(aligner, &n_max2, references, read2, consistent_nam2);
        if let (Some(alignment1), Some(alignment2)) = (alignment1, alignment2) {
            details[0].tried_alignment += 1;
            details[0].gapped += alignment1.gapped as usize;
            details[1].tried_alignment += 1;
            details[1].gapped += alignment2.gapped as usize;

            return AlignedPairs::Proper((alignment1, alignment2));
        }
        // TODO what if one of the alignments is None?
    }

    // Do a full search for highest-scoring pair
    // Get top hit counts for all locations.
    // The joint hit count is the sum of hits of the two mates.
    // Then align as long as score dropoff or cnt < 20

    let reads = [read1, read2];
    // Cache for already computed alignments. Maps NAM ids to alignments.
    // TODO rename
    let mut alignment_cache = [HashMap::new(), HashMap::new()];

    // These keep track of the alignments that would be best if we treated
    // the paired-end read as two single-end reads.
    let mut a_indv_max= [None, None];
    for i in 0..2 {
        let consistent_nam = reverse_nam_if_needed(&mut nams[i][0], reads[i], references, k);
        details[i].inconsistent_nams += !consistent_nam as usize;
        a_indv_max[i] = extend_seed(aligner, &nams[i][0], references, reads[i], consistent_nam);
        details[i].tried_alignment += 1;
        details[i].gapped += a_indv_max[i].as_ref().map_or(0, |a| a.score) as usize;
        alignment_cache[i].insert(nams[i][0].nam_id, a_indv_max[i].clone());
    }

    // Turn pairs of high-scoring NAMs into pairs of alignments
    let nam_pairs = get_best_scoring_nam_pairs(&nams[0], &nams[1], mu, sigma);
    let mut alignment_pairs = vec![];
    let max_score = nam_pairs[0].n_hits;
    for nam_pair in nam_pairs {
        let score_ = nam_pair.n_hits;
        let namsp = [nam_pair.nam1, nam_pair.nam2];
        let score_dropoff = score_ as f32 / max_score as f32;

        if alignment_pairs.len() >= max_tries || score_dropoff < dropoff {
            break;
        }

        // Get alignments for the two NAMs, either by computing the alignment,
        // retrieving it from the cache or by attempting a rescue (if the NAM
        // actually is a dummy, that is, only the partner is available)

        let mut alignments = [None, None];
        for i in 0..2 {
            let alignment;
            if let Some(mut this_nam) = namsp[i].clone() {
                if alignment_cache[i].contains_key(&this_nam.nam_id) {
                    alignment = alignment_cache[i].get(&this_nam.nam_id).unwrap().clone();
                } else {
                    let consistent_nam = reverse_nam_if_needed(&mut this_nam, reads[i], references, k);
                    details[i].inconsistent_nams += !consistent_nam as usize;
                    alignment = extend_seed(aligner, &this_nam, references, reads[i], consistent_nam);
                    details[i].tried_alignment += 1;
                    details[i].gapped += alignment.as_ref().map_or(false, |a| a.gapped) as usize;
                    alignment_cache[i].insert(this_nam.nam_id, alignment.clone());
                }
            } else {
                let mut other_nam = namsp[1-i].clone().unwrap();
                details[1-i].inconsistent_nams += !reverse_nam_if_needed(&mut other_nam, reads[1-i], references, k) as usize;
                alignment = rescue_align(aligner, &other_nam, references, reads[i], mu, sigma, k);
                if alignment.is_some() {
                    details[i].mate_rescue += 1;
                    details[i].tried_alignment += 1;
                }
            }
            if alignment.as_ref().map_or(0, |a| a.score) > a_indv_max[i].as_ref().map_or(0, |a| a.score) {
                a_indv_max[i] = alignment.clone();
            }
            alignments[i] = alignment.clone();
        }
        
        if alignments[0].is_none() || alignments[1].is_none() {
            continue;
        }
        let a1 = alignments[0].as_ref().unwrap();
        let a2 = alignments[1].as_ref().unwrap();
        let r1_r2 = a2.is_revcomp && !a1.is_revcomp && (a1.ref_start <= a2.ref_start) && ((a2.ref_start - a1.ref_start) < (mu + 10.0 * sigma) as usize); // r1 ---> <---- r2
        let r2_r1 = a1.is_revcomp && !a2.is_revcomp && (a2.ref_start <= a1.ref_start) && ((a1.ref_start - a2.ref_start) < (mu + 10.0 * sigma) as usize); // r2 ---> <---- r1

        let combined_score=
            if r1_r2 || r2_r1 {
                // Treat as a pair
                let x = a1.ref_start.abs_diff(a2.ref_start);
                a1.score as f64 + a2.score as f64 + (-20.0f64 + 0.001).max(normal_pdf(x as f32, mu, sigma).ln() as f64)
            } else {
                // Treat as two single-end reads
                // 20 corresponds to a value of log(normal_pdf(x, mu, sigma)) of more than 5 stddevs away (for most reasonable values of stddev)
                a1.score as f64 + a2.score as f64 - 20.0
            };

        let aln_pair = ScoredAlignmentPair { score: combined_score, alignment1: Some(a1.clone()), alignment2: Some(a2.clone()) };
        alignment_pairs.push(aln_pair);
    }

    // Finally, add highest scores of both mates as individually mapped
    // 20 corresponds to  a value of log( normal_pdf(x, mu, sigma ) ) of more than 5 stddevs away (for most reasonable values of stddev)
    if let (Some(a1), Some(a2)) = (&a_indv_max[0], &a_indv_max[1]) {
        let combined_score = a1.score as f64 + a2.score as f64 - 20.0;
        alignment_pairs.push(ScoredAlignmentPair {score: combined_score, alignment1: Some(a1.clone()), alignment2: Some(a2.clone())});
    }

    AlignedPairs::WithScores(alignment_pairs)
}

/// Align a pair of reads for which only one has NAMs. For the other, rescue
/// is attempted by aligning it locally.
fn rescue_read(
    read2: &Read,  // read to be rescued
    read1: &Read,  // read that has NAMs
    aligner: &Aligner,
    references: &[RefSequence],
    nams1: &mut [Nam],
    max_tries: usize,
    dropoff: f32,
    details: &mut [Details; 2],
    k: usize,
    mu: f32,
    sigma: f32
) -> Vec<ScoredAlignmentPair> {
    let n_max1_hits = nams1[0].n_hits;

    let mut alignments1 = vec![];
    let mut alignments2 = vec![];
    for nam in nams1.iter_mut().take(max_tries) {
        let score_dropoff1 = nam.n_hits as f32 / n_max1_hits as f32;
        // only consider top hits (as minimap2 does) and break if below dropoff cutoff.
        if score_dropoff1 < dropoff {
            break;
        }
        let consistent_nam = reverse_nam_if_needed(nam, read1, references, k);
        details[0].inconsistent_nams += !consistent_nam as usize;
        if let Some(alignment) = extend_seed(aligner, nam, references, read1, consistent_nam) {
            details[0].gapped += alignment.gapped as usize;
            alignments1.push(alignment);
            details[0].tried_alignment += 1;

            let a2 = rescue_align(aligner, nam, references, read2, mu, sigma, k);
            if a2.is_some() {
                details[1].mate_rescue += 1;
            }
            alignments2.push(a2);
        }
    }
    /*
    TODO This should not be necessary because we later sort the pairs by score
    alignments1.sort_by_key(|a| Reverse(a.score));
    alignments2.sort_by_key(|a| Reverse(a.score));
    */
    let mut pairs = vec![];
    for a1 in alignments1 {
        for a2 in &alignments2 {
            if let Some(a2) = a2 {
                let dist = a1.ref_start.abs_diff(a2.ref_start);
                let mut score = (a1.score + a2.score) as f32;
                if (a1.is_revcomp ^ a2.is_revcomp) && (dist as f32) < mu + 4.0 * sigma {
                    score += normal_pdf(dist as f32, mu, sigma).ln();
                } else {
                    // 10 corresponds to a value of log(normal_pdf(dist, mu, sigma)) of more than 4 stddevs away
                    score -= 10.0;
                }
                pairs.push(ScoredAlignmentPair {score: score as f64, alignment1: Some(a1.clone()), alignment2: Some(a2.clone())});
            } else {
                let score = (a1.score as f64) - 10.0;
                pairs.push(ScoredAlignmentPair {score, alignment1: Some(a1.clone()), alignment2: None});
            }
        }
    }

    pairs
}

/// Align a read to the reference given the mapping location of its mate.
fn rescue_align(
    aligner: &Aligner,
    mate_nam: &Nam,
    references: &[RefSequence],
    read: &Read,
    mu: f32,
    sigma: f32,
    k: usize,
) -> Option<Alignment> {
    let read_len = read.len();

    let (r_tmp, ref_start, ref_end) =
        if mate_nam.is_revcomp {
            (
                read.seq(),
                mate_nam.projected_ref_start().saturating_sub((mu + 5.0 * sigma) as usize),
                mate_nam.projected_ref_start() + read_len/2  // at most half read overlap
            )
        } else {
            (
                read.rc(), // mate is rc since fr orientation
                (mate_nam.ref_end + read_len - mate_nam.query_end).saturating_sub(read_len / 2),  // at most half read overlap
                mate_nam.ref_end + read_len - mate_nam.query_end + (mu + 5.0 * sigma) as usize
            )
        };

    let ref_len = references[mate_nam.ref_id].sequence.len();
    let ref_start = ref_start.min(ref_len);
    let ref_end = ref_end.min(ref_len);

    if ref_end < ref_start + k {
//        std::cerr << "RESCUE: Caught Bug3! ref start: " << ref_start << " ref end: " << ref_end << " ref len:  " << ref_len << std::endl;
        return None
    }
    let ref_segm = &references[mate_nam.ref_id].sequence[ref_start..ref_end];

    if !has_shared_substring(r_tmp, ref_segm, k) {
        return None
    }
    let info = aligner.align(r_tmp, ref_segm).unwrap();
    Some(Alignment {
        reference_id: mate_nam.ref_id,
        ref_start: ref_start + info.ref_start,
        edit_distance: info.edit_distance,
        soft_clip_left: info.query_start,
        soft_clip_right: read_len - info.query_end,
        score: info.score,
        length: info.ref_span(),
        cigar: info.cigar,
        is_revcomp: !mate_nam.is_revcomp,
        gapped: true,
    })
}

/// Determine (roughly) whether the read sequence has some l-mer (with l = k*2/3)
/// in common with the reference sequence
fn has_shared_substring(
    read_seq: &[u8],
    ref_seq: &[u8],
    k: usize,
) -> bool {
    let sub_size = 2 * k / 3;
    let step_size = k / 3;
    for i in (0..read_seq.len().saturating_sub(sub_size)).step_by(step_size) {
        let submer = &read_seq[i..i+sub_size];
        if memmem::find(ref_seq, submer).is_some() {
            return true;
        }
    }

    false
}

fn is_proper_pair(alignment1: Option<&Alignment>, alignment2: Option<&Alignment>, mu: f32, sigma: f32) -> bool {
    match (alignment1, alignment2) {
        (None, None) => false,
        (Some(_), None) => false,
        (None, Some(_)) => false,
        (Some(a1), Some(a2)) => {
            let dist = a2.ref_start as isize - a1.ref_start as isize;
            let same_reference = a1.reference_id == a2.reference_id;
            let r1_r2 = !a1.is_revcomp && a2.is_revcomp && dist >= 0; // r1 ---> <---- r2
            let r2_r1 = !a2.is_revcomp && a1.is_revcomp && dist <= 0; // r2 ---> <---- r1
            let rel_orientation_good = r1_r2 || r2_r1;
            let insert_good = dist.unsigned_abs() <= (mu + sigma * 6.0) as usize;

            same_reference && insert_good && rel_orientation_good
        }
    }
}

fn is_proper_nam_pair(nam1: &Nam, nam2: &Nam, mu: f32, sigma: f32) -> bool {
    if nam1.ref_id != nam2.ref_id || nam1.is_revcomp == nam2.is_revcomp {
        return false;
    }
    let r1_ref_start = nam1.projected_ref_start();
    let r2_ref_start = nam2.projected_ref_start();

    // r1 ---> <---- r2
    let r1_r2 = nam2.is_revcomp && (r1_ref_start <= r2_ref_start) && ((r2_ref_start - r1_ref_start) as f32) < mu + 10.0 * sigma;

     // r2 ---> <---- r1
    let r2_r1 = nam1.is_revcomp && (r2_ref_start <= r1_ref_start) && ((r1_ref_start - r2_ref_start) as f32) < mu + 10.0 * sigma;

    r1_r2 || r2_r1
}

/// Find high-scoring NAMs and NAM pairs. Proper pairs are preferred, but also
/// high-scoring NAMs that could not be paired up are returned (these are paired
/// with None in the returned vector).
pub fn get_best_scoring_nam_pairs(
    nams1: &[Nam],
    nams2: &[Nam],
    mu: f32,
    sigma: f32,
) -> Vec<NamPair> {
    let mut nam_pairs = vec![];
    if nams1.is_empty() && nams2.is_empty() {
        return nam_pairs;
    }

    // Find NAM pairs that appear to be proper pairs
    let mut added_n1 = HashSet::new();
    let mut added_n2 = HashSet::new();
    let mut best_joint_hits = 0;
    for nam1 in nams1 {
        for nam2 in nams2 {
            let joint_hits = nam1.n_hits + nam2.n_hits;
            if joint_hits < best_joint_hits / 2 {
                break;
            }
            if is_proper_nam_pair(nam1, nam2, mu, sigma) {
                nam_pairs.push(NamPair{n_hits: joint_hits, nam1: Some(nam1.clone()), nam2: Some(nam2.clone())});
                added_n1.insert(nam1.nam_id);
                added_n2.insert(nam2.nam_id);
                best_joint_hits = joint_hits.max(best_joint_hits);
            }
        }
    }

    // Find high-scoring R1 NAMs that are not part of a proper pair
    if !nams1.is_empty() {
        let best_joint_hits1 = if best_joint_hits > 0 { best_joint_hits } else { nams1[0].n_hits };
        for nam1 in nams1 {
            if nam1.n_hits < best_joint_hits1 / 2 {
                break;
            }
            if added_n1.contains(&nam1.nam_id) {
                continue;
            }
            nam_pairs.push(NamPair{n_hits: nam1.n_hits, nam1: Some(nam1.clone()), nam2: None});
        }
    }

    // Find high-scoring R2 NAMs that are not part of a proper pair
    if !nams2.is_empty() {
        let best_joint_hits2 = if best_joint_hits > 0 { best_joint_hits } else { nams2[0].n_hits };
        for nam2 in nams2 {
            if nam2.n_hits < best_joint_hits2 / 2 {
                break;
            }
            if added_n2.contains(&nam2.nam_id) {
                continue;
            }
            nam_pairs.push(NamPair{n_hits: nam2.n_hits, nam1: None, nam2: Some(nam2.clone())});
        }
    }
    nam_pairs.sort_by_key(|nam_pair| Reverse(nam_pair.n_hits));

    nam_pairs
}

fn proper_pair_mapq(nams: &[Nam]) -> u8 {
    if nams.len() <= 1 {
        return 60;
    }
    let s1 = nams[0].score;
    let s2 = nams[1].score;
    // from minimap2: MAPQ = 40(1−s2/s1) ·min{1,|M|/10} · log s1
    let min_matches = min(nams[0].n_hits, 10) as f32 / 10.0;
    let uncapped_mapq = 40.0 * (1 - s2 / s1) as f32 * min_matches * (s1 as f32).ln();

    uncapped_mapq.min(60.0) as u8
}

#[derive(Debug)]
pub struct NamPair {
    pub n_hits: usize,
    pub nam1: Option<Nam>,
    pub nam2: Option<Nam>,
}

#[derive(Debug, Clone)]
pub struct ScoredAlignmentPair {
    score: f64,
    alignment1: Option<Alignment>,
    alignment2: Option<Alignment>,
}

/// Remove consecutive identical alignment pairs and leave only the first.
fn deduplicate_scored_pairs(pairs: &mut Vec<ScoredAlignmentPair>) {
    // TODO use Vec::dedup(_by...)
    if pairs.len() < 2 {
        return;
    }

    fn is_same(a1: Option<&Alignment>, a2: Option<&Alignment>) -> bool {
        match (a1, a2) {
            (None, Some(_)) => false,
            (Some(_), None) => false,
            (None, None) => true,
            (Some(a1), Some(a2)) => {
                a1.ref_start == a2.ref_start && a1.reference_id == a2.reference_id
            },
        }
    }
    let mut k = 0;
    let mut j = 1;
    for i in 1..pairs.len() {
        if !is_same(pairs[k].alignment1.as_ref(), pairs[i].alignment1.as_ref()) || !is_same(pairs[k].alignment2.as_ref(), pairs[i].alignment2.as_ref()) {
            pairs[j] = pairs[i].clone();
            j += 1;
            k = i;
        }
    }
    pairs.truncate(j);
}

fn count_best_alignment_pairs(pairs: &[ScoredAlignmentPair]) -> usize {
    if pairs.is_empty() {
        0
    } else {
        pairs.iter().take_while(|x| x.score == pairs[0].score).count()
    }
}

fn aligned_pairs_to_sam(
    high_scores: &[ScoredAlignmentPair],
    sam_output: &SamOutput,
    references: &[RefSequence],
    max_secondary: usize,
    secondary_dropoff: f64,
    record1: &SequenceRecord,
    record2: &SequenceRecord,
    mu: f32,
    sigma: f32,
    details: &[Details; 2],
) -> Vec<SamRecord> {

    let mut records = vec![];
    if high_scores.is_empty() {
        records.extend(sam_output.make_unmapped_pair([record1, record2], details));
        return records;
    }

    let mapq = joint_mapq_from_high_scores(high_scores);
    let best_aln_pair = &high_scores[0];

    if max_secondary == 0 {
        let alignment1 = &best_aln_pair.alignment1;
        let alignment2 = &best_aln_pair.alignment2;

        let is_proper = is_proper_pair(alignment1.as_ref(), alignment2.as_ref(), mu, sigma);
        let is_primary = true;
        records.extend(sam_output.make_paired_records(
            [alignment1.as_ref(), alignment2.as_ref()], references, [record1, record2], [mapq, mapq], details, is_primary, is_proper
        ));
    } else {
        let mut is_primary = true;
        let s_max = best_aln_pair.score;
        for aln_pair in high_scores.iter().take(max_secondary) {
            let alignment1 = &aln_pair.alignment1;
            let alignment2 = &aln_pair.alignment2;
            let s_score = aln_pair.score;
            if s_max - s_score < secondary_dropoff {
                let is_proper = is_proper_pair(alignment1.as_ref(), alignment2.as_ref(), mu, sigma);
                let mapq = if is_primary { mapq } else { 0 };
                records.extend(sam_output.make_paired_records(
                    [alignment1.as_ref(), alignment2.as_ref()], references, [record1, record2], [mapq, mapq], details, is_proper, is_primary
                ));
            } else {
                break;
            }
            is_primary = false;
        }
    }

    records
}

fn joint_mapq_from_high_scores(pairs: &[ScoredAlignmentPair]) -> u8 {
    if pairs.len() <= 1 {
        return 60;
    }
    let score1 = pairs[0].score;
    let score2 = pairs[1].score;
    if score1 == score2 {
        return 0;
    }
    if score1 > 0.0 && score2 > 0.0 {
        (score1 - score2).min(60.0) as u8
    } else if score1 > 0.0 && score2 <= 0.0 {
        60
    } else {
        1
    }
}


/// compute dropoff of the first (top) NAM
fn top_dropoff(nams: &[Nam]) -> f32 {
    let n_max = &nams[0];
    if n_max.n_hits <= 2 {
        1.0
    } else if nams.len() > 1 {
        nams[1].n_hits as f32 / n_max.n_hits as f32
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use crate::cigar::Cigar;
    use crate::mapper::{count_best_alignment_pairs, deduplicate_scored_pairs, Alignment, ScoredAlignmentPair};

    fn dummy_alignment() -> Alignment {
        Alignment {
            reference_id: 0,
            ref_start: 0,
            cigar: Cigar::default(),
            edit_distance: 0,
            soft_clip_left: 0,
            soft_clip_right: 0,
            score: 0,
            length: 0,
            is_revcomp: false,
            gapped: false,
        }
    }

    #[test]
    fn test_count_best_alignment_pairs() {
        let mut pairs = vec![];
        fn add_alignment(pairs: &mut Vec<ScoredAlignmentPair>, score: f64) {
            pairs.push(ScoredAlignmentPair { score, alignment1: Some(dummy_alignment()), alignment2: Some(dummy_alignment()) });
        }

        assert_eq!(count_best_alignment_pairs(&pairs), 0);
        add_alignment(&mut pairs, 10.0);
        assert_eq!(count_best_alignment_pairs(&pairs), 1);

        add_alignment(&mut pairs, 10.0);
        assert_eq!(count_best_alignment_pairs(&pairs), 2);

        add_alignment(&mut pairs, 5.0);
        assert_eq!(count_best_alignment_pairs(&pairs), 2);

        pairs[1].score = 5.0;
        assert_eq!(count_best_alignment_pairs(&pairs), 1);
    }

    #[test]
    fn test_deduplicate_scored_pairs() {
        let a1 = Some(
            Alignment {
                reference_id: 0,
                ref_start: 1906,
                .. dummy_alignment()
            },
        );
        let a2 = Some(
            Alignment {
                reference_id: 0,
                ref_start: 123,
                .. dummy_alignment()
            },
        );
        let mut alignment_pairs = vec![
            ScoredAlignmentPair {
                score: 733.0,
                alignment1: a1.clone(),
                alignment2: a2.clone(),
            },
            ScoredAlignmentPair {
                score: 724.0,
                alignment1: a1.clone(),
                alignment2: a2.clone(),
            },
        ];
        deduplicate_scored_pairs(&mut alignment_pairs);
        assert_eq!(alignment_pairs.len(), 1);
    }
}
